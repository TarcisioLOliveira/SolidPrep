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

#include <limits>
#include <mpich-x86_64/mpi.h>
#include "function/density_based/mechanostat.hpp"
#include "logger.hpp"
#include "math/matrix.hpp"
#include "utils.hpp"
#include "optimizer.hpp"
#include "project_data.hpp"

namespace function::density_based{

Mechanostat::Mechanostat(const projspec::DataMap& data):
    mesh(data.proj->topopt_mesher.get()),
    fem(data.proj->topopt_fea.get()),
    beta(data.get_double("beta")),
    pc(data.proj->topopt_penalization),
    psiK(data.proj->topopt_psi),
    t({
        data.get_array("traction")->get_double(0),
        data.get_array("traction")->get_double(1)
    }),
    c({
        data.get_array("compression")->get_double(0),
        data.get_array("compression")->get_double(1)
    }),
    s({
        data.get_array("shear")->get_double(0),
        data.get_array("shear")->get_double(1)
    }),
    K_e1({0.5*std::pow((t[0]+c[0])/(2*t[0]*c[0]), 2), 
          0.5*std::pow((t[1]+c[1])/(2*t[1]*c[1]), 2)}),
    K_g({1.0/(2*s[0]*s[0]),
         1.0/(2*s[1]*s[1])}),
    K_e2({-(t[0]-c[0])/(2*t[0]*c[0]), 
          -(t[1]-c[1])/(2*t[1]*c[1])}),
    problem_type(data.proj->type){

}

void Mechanostat::initialize_views(Visualization* viz){
    this->shadow_view = viz->add_view("Normalized Strain", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::OTHER);
    this->gradient_view = viz->add_view("Normalized Strain Gradient", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::OTHER);
}

void Mechanostat::initialize(const DensityBasedOptimizer* const op){
    (void)op;
    size_t elem_num = 0;
    for(const auto& g:mesh->geometries){
        elem_num += g->mesh.size();
    }
    this->He.resize(elem_num, std::numeric_limits<double>::quiet_NaN());
    this->gradHe.resize(elem_num, std::numeric_limits<double>::quiet_NaN());
}

double Mechanostat::calculate(const DensityBasedOptimizer* const op, const std::vector<double>& u, const std::vector<double>& x){
    (void)op;
    (void)u;
    (void)x;
    return 0;
}
double Mechanostat::calculate_with_gradient(const DensityBasedOptimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
    (void)op;
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
    std::vector<std::vector<double>> fl(this->mesh->sub_problems->size());
    std::vector<double> l(mesh->max_dofs,0);
    for(size_t i = 0; i < fl.size(); ++i){
        fl[i].resize(mesh->max_dofs,0);
    }

    double result = 0;

    const size_t s_size = this->mesh->elem_info->get_D_dimension();
    const size_t dof = this->mesh->elem_info->get_dof_per_node();
    const size_t num_nodes = this->mesh->elem_info->get_nodes_per_element();

    if(mpi_id == 0){
        auto x_it = x.cbegin();
        auto He_it = this->He.begin();
        for(const auto& g:this->mesh->geometries){
            if(g->do_topopt){
                const size_t num_den = g->number_of_densities_needed();
                for(const auto& e:g->mesh){
                    const auto c = e->get_centroid();

                    const auto B = e->get_B(c);
                    const double rho = this->relaxed_rho(*x_it);

                    double H_e = 0;
                    for(size_t it = 0; it < fl.size(); ++it){
                        const auto& ui = fem->sub_u[it];
                        const math::Vector eps(1e6*e->get_strain_vector(c, ui));
                        if(this->problem_type == utils::PROBLEM_TYPE_2D){
                            const double eps_lhs1 = rho*this->LHS_2D(0, eps);
                            const double eps_lhs2 = rho*this->LHS_2D(1, eps);
                            const StrainVector2D deps_lhs1 = rho*this->dLHS_2D(0, eps);
                            const StrainVector2D deps_lhs2 = rho*this->dLHS_2D(1, eps);
                            const double H1 = Hm(eps_lhs1);
                            const double H2 = Hp(eps_lhs2);
                            const double Hr = H1*eps_lhs1 + H2*eps_lhs2;
                            const StrainVector2D dH = (dHm(eps_lhs1)*eps_lhs1 + H1)*deps_lhs1 +
                                                      (dHp(eps_lhs2)*eps_lhs2 + H2)*deps_lhs2;

                            H_e = Hr;
                            math::Vector dHB(dH.T()*B);
                            for(size_t i = 0; i < num_nodes; ++i){
                                for(size_t j = 0; j < dof; ++j){
                                    const long pos = e->nodes[i]->u_pos[j];
                                    if(pos > -1){
                                        fl[it][pos] += dHB[i*dof + j];
                                    }
                                }
                            }
                        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                            const double eps_lhs1 = rho*this->LHS_3D(0, eps);
                            const double eps_lhs2 = rho*this->LHS_3D(1, eps);
                            const StrainVector3D deps_lhs1 = rho*this->dLHS_3D(0, eps);
                            const StrainVector3D deps_lhs2 = rho*this->dLHS_3D(1, eps);
                            const double H1 = Hm(eps_lhs1);
                            const double H2 = Hp(eps_lhs2);
                            const double Hr = H1*eps_lhs1 + H2*eps_lhs2;
                            const StrainVector3D dH = (dHm(eps_lhs1)*eps_lhs1 + H1)*deps_lhs1 +
                                                      (dHp(eps_lhs2)*eps_lhs2 + H2)*deps_lhs2;

                            H_e = Hr;
                            math::Vector dHB(dH.T()*B);

                            for(size_t i = 0; i < num_nodes; ++i){
                                for(size_t j = 0; j < dof; ++j){
                                    const long pos = e->nodes[i]->u_pos[j];
                                    if(pos > -1){
                                        fl[it][pos] += dHB[i*dof + j];
                                    }
                                }
                            }
                        }
                    }

                    *He_it = H_e;
                    result += H_e;

                    x_it += num_den;
                    ++He_it;
                }
            } else {
                He_it += g->mesh.size();
            }
        }
        logger::quick_log("Calculating adjoint problem...{");
    }
    this->shadow_view->update_view(this->He);
    this->fem->calculate_displacements_adjoint(this->mesh, fl, l);
    if(mpi_id == 0){
        logger::quick_log("} Done.");

        logger::quick_log("Calculating strain gradient...");

        size_t xi = 0;
        size_t gi = 0;
        size_t ghi = 0;
        for(const auto& g:this->mesh->geometries){
            if(g->do_topopt){
                const size_t num_den = g->number_of_densities_needed();
                std::vector<math::Matrix> gradD_K(num_den, math::Matrix(s_size, s_size));
                if(g->with_void){
                    #pragma omp parallel for
                    for(size_t i = 0; i < g->mesh.size(); ++i){
                        const auto& e = g->mesh[i];
                        const auto c = e->get_centroid();

                        auto x_it = x.cbegin() + xi + i*num_den;
                        auto grad_it = grad.begin() + gi + i*num_den;
                        auto gradHe_it = this->gradHe.begin() + ghi + i;

                        g->materials.get_gradD(x_it, psiK, e.get(), c, gradD_K);
                        double lKu = pc*std::pow(*x_it, pc-1)*e->get_compliance(gradD_K[0], this->mesh->thickness, u, l);

                        const math::Vector eps(1e6*e->get_strain_vector(c, u));

                        const double rho = this->relaxed_rho(*x_it);
                        const double drho = this->relaxed_rho_grad(*x_it);

                        double dH_e = 0;
                        if(this->problem_type == utils::PROBLEM_TYPE_2D){
                            const double eps_lhs1 = this->LHS_2D(0, eps);
                            const double eps_lhs2 = this->LHS_2D(1, eps);
                            const double H1 = Hm(rho*eps_lhs1);
                            const double H2 = Hp(rho*eps_lhs2);
                            const auto dH = drho*((dHm(rho*eps_lhs1)*rho*eps_lhs1 + H1)*eps_lhs1 +
                                                  (dHp(rho*eps_lhs2)*rho*eps_lhs2 + H2)*eps_lhs2);

                            dH_e = dH;
                        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                            const double eps_lhs1 = this->LHS_3D(0, eps);
                            const double eps_lhs2 = this->LHS_3D(1, eps);
                            const double H1 = Hm(rho*eps_lhs1);
                            const double H2 = Hp(rho*eps_lhs2);
                            const auto dH = drho*((dHm(rho*eps_lhs1)*rho*eps_lhs1 + H1)*eps_lhs1 +
                                                  (dHp(rho*eps_lhs2)*rho*eps_lhs2 + H2)*eps_lhs2);

                            dH_e = dH;
                        }

                        *grad_it = dH_e - lKu;
                        *gradHe_it = *grad_it;

                        const double rho_lKu = std::pow(*x_it, pc);

                        ++x_it;
                        ++grad_it;
                        ++gradHe_it;
                        for(size_t j = 1; j < num_den; ++j){
                            double lKu = rho_lKu*e->get_compliance(gradD_K[j], this->mesh->thickness, u, l);

                            *grad_it = -lKu;

                            ++x_it;
                            ++grad_it;
                        }
                    }
                } else {
                    #pragma omp parallel for
                    for(size_t i = 0; i < g->mesh.size(); ++i){
                        const auto& e = g->mesh[i];
                        const auto c = e->get_centroid();

                        auto x_it = x.cbegin() + xi + i*num_den;
                        auto grad_it = grad.begin() + gi + i*num_den;
                        auto gradHe_it = this->gradHe.begin() + ghi + i;

                        const math::Vector eps(1e6*e->get_strain_vector(c, u));

                        const double rho = this->relaxed_rho(*x_it);
                        const double drho = this->relaxed_rho_grad(*x_it);

                        double dH_e = 0;
                        if(this->problem_type == utils::PROBLEM_TYPE_2D){
                            const double eps_lhs1 = this->LHS_2D(0, eps);
                            const double eps_lhs2 = this->LHS_2D(1, eps);
                            const double H1 = Hm(rho*eps_lhs1);
                            const double H2 = Hp(rho*eps_lhs2);
                            const auto dH = drho*((dHm(rho*eps_lhs1)*rho*eps_lhs1 + H1)*eps_lhs1 +
                                                  (dHp(rho*eps_lhs2)*rho*eps_lhs2 + H2)*eps_lhs2);

                            dH_e = dH;
                        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                            const double eps_lhs1 = this->LHS_3D(0, eps);
                            const double eps_lhs2 = this->LHS_3D(1, eps);
                            const double H1 = Hm(rho*eps_lhs1);
                            const double H2 = Hp(rho*eps_lhs2);
                            const auto dH = drho*((dHm(rho*eps_lhs1)*rho*eps_lhs1 + H1)*eps_lhs1 +
                                                  (dHp(rho*eps_lhs2)*rho*eps_lhs2 + H2)*eps_lhs2);

                            dH_e = dH;
                        }

                        g->materials.get_gradD(x_it, psiK, e.get(), c, gradD_K);
                        double lKu = e->get_compliance(gradD_K[0], this->mesh->thickness, u, l);

                        *grad_it = dH_e - lKu;

                        ++x_it;
                        ++grad_it;
                        for(size_t j = 1; j < num_den; ++j){
                            double lKu = e->get_compliance(gradD_K[j], this->mesh->thickness, u, l);

                            *grad_it = -lKu;

                            ++x_it;
                            ++grad_it;
                        }
                        *gradHe_it = *(grad_it - num_den);
                    }
                }
                xi += g->mesh.size()*num_den;
                gi += g->mesh.size()*num_den;
                ghi += g->mesh.size();
            }
        }
        logger::quick_log("Done.");
    }
    this->gradient_view->update_view(this->gradHe);

    return result;
}

using namespace projspec;
const bool Mechanostat::reg = Factory<DensityBasedFunction>::add(
    [](const DataMap& data){
        return std::make_unique<Mechanostat>(data);
    },
    ObjectRequirements{
        "mechanostat",
        {
            DataEntry{.name = "beta", .type = TYPE_STRING, .required = true},
            DataEntry{.name = "traction", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 2,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
            DataEntry{.name = "compression", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 2,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
            DataEntry{.name = "shear", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 2,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
        }
    }
);

}
