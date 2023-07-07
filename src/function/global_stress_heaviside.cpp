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

#include "function/global_stress_heaviside.hpp"
#include <algorithm>
#include <mpich-x86_64/mpi.h>

namespace function{

GlobalStressHeaviside::GlobalStressHeaviside(const Meshing* const mesh, FiniteElement* fem, double max_stress, double C, double pc, double pt, double psiK, double psiS):
    mesh(mesh), fem(fem), max_stress(max_stress), C(C), pc(pc), pt(pt), psiK(psiK), psiS(psiS){}

double GlobalStressHeaviside::calculate(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x){
    (void)x;

    return 0;
}

double GlobalStressHeaviside::calculate_with_gradient(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
    std::vector<double> fl(u.size(), 0);

    double result = 0;

    auto x_it = x.cbegin();
    if(mpi_id == 0){
        for(const auto& g:this->mesh->geometries){
            const size_t num_mat = g->number_of_materials();
            if(g->do_topopt){
                if(num_mat == 1){
                    const auto D = g->materials.get_D();
                    for(const auto& e:g->mesh){
                        const auto c = e->get_centroid();
                        const double S = e->get_stress_at(D, c, u, this->vm_eps);
                        const double rho = this->relaxed_rho(*x_it);
                        const double Se = rho*S;
                        const double H = this->heaviside(Se);
                        const double dH = this->heaviside_grad(Se);
                        e->get_virtual_load(D, (dH*rho*Se + H)*rho/S, e->get_centroid(), u, fl);

                        result += H*Se;

                        ++x_it;
                    }
                } else {
                    auto D = g->materials.get_D();
                    for(const auto& e:g->mesh){
                        g->materials.get_D(x_it, g->with_void, pt, this->K_MIN, this->psiS, D);
                        const auto c = e->get_centroid();
                        const double Se = e->get_stress_at(D, c, u, this->vm_eps);
                        const double H = this->heaviside(Se);
                        const double dH = this->heaviside_grad(Se);
                        e->get_virtual_load(D, (dH*Se + H)/Se, e->get_centroid(), u, fl);

                        result += H*Se;
                    }
                }
            } else {
                const auto D = g->materials.get_D();
                for(const auto& e:g->mesh){
                    const auto c = e->get_centroid();
                    const double Se = e->get_stress_at(D, c, u, this->vm_eps);
                    const double H = this->heaviside(Se);
                    const double dH = this->heaviside_grad(Se);
                    e->get_virtual_load(D, (dH*Se + H)/Se, e->get_centroid(), u, fl);
                    result += this->heaviside(Se)*Se;
                }
            }
        }
        logger::quick_log("Calculating adjoint problem...{");
    }
    auto l = this->fem->calculate_displacements(this->mesh, fl, x, pc);
    if(mpi_id == 0){
        logger::quick_log("} Done.");

        logger::quick_log("Calculating stress gradient...");

        x_it = x.cbegin();
        auto grad_it = grad.begin();
        for(const auto& g:this->mesh->geometries){
            const size_t num_den = g->number_of_densities_needed();
            const size_t num_mat = g->number_of_materials();
            if(g->do_topopt){
                if(num_mat == 1){
                    const auto D = g->materials.get_D();
                    for(const auto& e:g->mesh){
                        double lKu = pc*std::pow(*x_it, pc-1)*e->get_compliance(D, this->mesh->thickness, u, l);
                        const auto c = e->get_centroid();
                        const double S = e->get_stress_at(D, c, u, this->vm_eps);
                        const double rho = this->relaxed_rho(*x_it);
                        const double drho = this->relaxed_rho_grad(*x_it);
                        const double Se = rho*S;
                        const double H = this->heaviside(Se);
                        const double dH = this->heaviside_grad(Se);

                        *grad_it = (dH*rho*Se + H)*drho*Se - lKu;

                        ++x_it;
                        ++grad_it;
                    }
                } else {
                    std::vector<std::vector<double>> gradD_K(num_den, g->materials.get_D());
                    std::vector<std::vector<double>> gradD_S(num_den, g->materials.get_D());
                    std::vector<double> D_S(g->materials.get_D());
                    for(const auto& e:g->mesh){
                        g->materials.get_gradD(x_it, g->with_void, pc, this->K_MIN, psiK, gradD_K);
                        x_it -= num_den;
                        g->materials.get_gradD(x_it, g->with_void, pt, this->K_MIN, psiS, gradD_S);
                        x_it -= num_den;
                        g->materials.get_D(x_it, g->with_void, pt, this->K_MIN, this->psiS, D_S);
                        x_it -= num_den;
                        for(size_t i = 0; i < num_den; ++i){
                            double lKu = e->get_compliance(gradD_K[i], this->mesh->thickness, u, l);
                            const auto c = e->get_centroid();
                            const double S = e->get_stress_at(D_S, c, u, this->vm_eps);
                            const double Se = e->von_Mises_derivative(D_S, gradD_S[i], this->heaviside_grad(S)/S, c, u);
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

}

