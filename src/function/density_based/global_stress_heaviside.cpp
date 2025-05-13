/*
 *   Copyright (C) 2023 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   SolidPrep is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   SolidPrep is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with SolidPrep.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <mpich-x86_64/mpi.h>
#include "function/density_based/global_stress_heaviside.hpp"
#include "logger.hpp"
#include "math/matrix.hpp"
#include "optimizer.hpp"
#include "project_specification/data_map.hpp"
#include "project_data.hpp"

namespace function::density_based{

GlobalStressHeaviside::GlobalStressHeaviside(const projspec::DataMap& data):
    mesh(data.proj->topopt_mesher.get()),
    fem(data.proj->topopt_fea.get()),
    max_stress(data.get_double("max_stress")),
    C(data.get_double("C")),
    pc(data.proj->topopt_penalization),
    pt(data.get_double("pt")),
    psiK(data.proj->topopt_psi),
    psiS(-data.get_double("psiS")){}

double GlobalStressHeaviside::calculate(const DensityBasedOptimizer* const op, const std::vector<double>& u, const std::vector<double>& x){
    (void)op;
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    double result = 0;

    const size_t s_size = this->mesh->elem_info->get_D_dimension();
    auto x_it = x.cbegin();
    if(mpi_id == 0){
        for(const auto& g:this->mesh->geometries){
            if(g->do_topopt){
                const size_t num_den = g->number_of_densities_needed();
                math::Matrix D(s_size, s_size);
                if(g->with_void){
                    for(const auto& e:g->mesh){
                        const auto c = e->get_centroid();
                        g->materials.get_D(x_it, this->psiS, e.get(), c, D);
                        const double S = e->get_stress_at(D, c, u, this->vm_eps);
                        const double rho = this->relaxed_rho(*x_it);
                        const double Se = rho*S;
                        const double H = this->heaviside(Se);

                        result += H*Se;

                        x_it += num_den;
                    }
                } else {
                    for(const auto& e:g->mesh){
                        const auto c = e->get_centroid();
                        g->materials.get_D(x_it, this->psiS, e.get(), c, D);
                        const double Se = e->get_stress_at(D, c, u, this->vm_eps);
                        const double H = this->heaviside(Se);

                        result += H*Se;
                        x_it += num_den;
                    }
                }
            } else {
                for(const auto& e:g->mesh){
                    const auto c = e->get_centroid();
                    const auto D = g->materials.get_D(e.get(), c);
                    const double Se = e->get_stress_at(D, c, u, this->vm_eps);
                    const double H = this->heaviside(Se);
                    result += H*Se;
                }
            }
        }
        logger::quick_log("Calculating adjoint problem...{");
    }

    return result;
}

double GlobalStressHeaviside::calculate_with_gradient(const DensityBasedOptimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
    (void) op;
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
    std::vector<std::vector<double>> fl(this->mesh->sub_problems->size());
    std::vector<double> l(mesh->max_dofs,0);
    for(size_t i = 0; i < fl.size(); ++i){
        fl[i].resize(mesh->max_dofs);
    }

    double result = 0;

    const size_t s_size = this->mesh->elem_info->get_D_dimension();
    auto x_it = x.cbegin();
    if(mpi_id == 0){
        for(const auto& g:this->mesh->geometries){
            if(g->do_topopt){
                const size_t num_den = g->number_of_densities_needed();
                math::Matrix D(s_size, s_size);
                if(g->with_void){
                    for(const auto& e:g->mesh){
                        const auto c = e->get_centroid();
                        g->materials.get_D(x_it, this->psiS, e.get(), c, D);
                        const double S = e->get_stress_at(D, c, u, this->vm_eps);
                        const double rho = this->relaxed_rho(*x_it);
                        const double Se = rho*S;
                        const double H = this->heaviside(Se);
                        const double dH = this->heaviside_grad(Se);
                        for(size_t i = 0; i < fl.size(); ++i){
                            const auto& ui = fem->sub_u[i];
                            e->get_virtual_load(D, (dH*rho*Se + H)*rho/S, e->get_centroid(), ui, fl[i]);
                        }
                        result += H*Se;

                        x_it += num_den;
                    }
                } else {
                    for(const auto& e:g->mesh){
                        const auto c = e->get_centroid();
                        g->materials.get_D(x_it, this->psiS, e.get(), c, D);
                        const double Se = e->get_stress_at(D, c, u, this->vm_eps);
                        const double H = this->heaviside(Se);
                        const double dH = this->heaviside_grad(Se);
                        for(size_t i = 0; i < fl.size(); ++i){
                            const auto& ui = fem->sub_u[i];
                            e->get_virtual_load(D, (dH*Se + H)/Se, e->get_centroid(), ui, fl[i]);
                        }
                        result += H*Se;

                        x_it += num_den;
                    }
                }
            } else {
                for(const auto& e:g->mesh){
                    const auto c = e->get_centroid();
                    const auto D = g->materials.get_D(e.get(), c);
                    const double Se = e->get_stress_at(D, c, u, this->vm_eps);
                    const double H = this->heaviside(Se);
                    const double dH = this->heaviside_grad(Se);
                    for(size_t i = 0; i < fl.size(); ++i){
                        const auto& ui = fem->sub_u[i];
                        e->get_virtual_load(D, (dH*Se + H)/Se, e->get_centroid(), ui, fl[i]);
                    }
                    result += H*Se;
                }
            }
        }
        logger::quick_log("Calculating adjoint problem...{");
    }
    this->fem->calculate_displacements_adjoint(this->mesh, fl, l);
    if(mpi_id == 0){
        logger::quick_log("} Done.");

        logger::quick_log("Calculating stress gradient...");

        x_it = x.cbegin();
        auto grad_it = grad.begin();
        for(const auto& g:this->mesh->geometries){
            if(g->do_topopt){
                const size_t num_den = g->number_of_densities_needed();
                std::vector<math::Matrix> gradD_K(num_den, math::Matrix(s_size, s_size));
                std::vector<math::Matrix> gradD_S(num_den, math::Matrix(s_size, s_size));
                math::Matrix D_S(math::Matrix(s_size, s_size));
                if(g->with_void){
                    for(const auto& e:g->mesh){
                        const auto c = e->get_centroid();
                        g->materials.get_gradD(x_it, psiK, e.get(), c, gradD_K);
                        g->materials.get_gradD(x_it, psiS, e.get(), c, gradD_S);
                        g->materials.get_D(x_it, psiS, e.get(), c, D_S);

                        double lKu = pc*std::pow(*x_it, pc-1)*e->get_compliance(gradD_K[0], this->mesh->thickness, u, l);
                        const double S = e->get_stress_at(D_S, c, u, this->vm_eps);
                        const double rho = this->relaxed_rho(*x_it);
                        const double drho = this->relaxed_rho_grad(*x_it);
                        const double Se = rho*S;
                        const double H = this->heaviside(Se);
                        const double dH = this->heaviside_grad(Se);

                        *grad_it = (dH*rho*Se + H)*drho*Se - lKu;

                        const double mult = rho*this->heaviside_grad(S)/S;
                        const double rho_lKu = std::pow(*x_it, pc);

                        ++x_it;
                        ++grad_it;
                        for(size_t i = 1; i < num_den; ++i){
                            double lKu = rho_lKu*e->get_compliance(gradD_K[i], this->mesh->thickness, u, l);
                            const auto c = e->get_centroid();
                            const double Se = e->von_Mises_derivative(D_S, gradD_S[i], mult, c, u);
                            const double H = this->heaviside(Se);
                            const double dH = this->heaviside_grad(Se);

                            *grad_it = (dH*S + H)*Se - lKu;

                            ++x_it;
                            ++grad_it;
                        }
                    }
                } else {
                    for(const auto& e:g->mesh){
                        const auto c = e->get_centroid();
                        g->materials.get_gradD(x_it, psiK, e.get(), c, gradD_K);
                        g->materials.get_gradD(x_it, psiS, e.get(), c, gradD_S);
                        g->materials.get_D(x_it, psiS, e.get(), c, D_S);

                        const double S = e->get_stress_at(D_S, c, u, this->vm_eps);
                        const double mult = this->heaviside_grad(S)/S;
                        for(size_t i = 0; i < num_den; ++i){
                            double lKu = e->get_compliance(gradD_K[i], this->mesh->thickness, u, l);
                            const double Se = e->von_Mises_derivative(D_S, gradD_S[i], mult, c, u);
                            const double H = this->heaviside(Se);
                            const double dH = this->heaviside_grad(Se);

                            *grad_it = (dH*S + H)*Se - lKu;

                            ++x_it;
                            ++grad_it;
                        }
                    }
                }
            }
        }
        logger::quick_log("Done.");
    }

    return result;
}

using namespace projspec;
const bool GlobalStressHeaviside::reg = Factory<DensityBasedFunction>::add(
    [](const DataMap& data){
        return std::make_unique<GlobalStressHeaviside>(data);
    },
    ObjectRequirements{
        "global_stress_heaviside",
        {
            DataEntry{.name = "C", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "max_stress", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "pt", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "psi", .type = TYPE_DOUBLE, .required = true},
        }
    }
);

}

