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

#include "function/mechanostat.hpp"
#include "utils.hpp"
#include <Eigen/src/Core/Map.h>
#include <mpich-x86_64/mpi.h>

namespace function{

Mechanostat::Mechanostat(const Meshing* const mesh, FiniteElement* fem, double pc, double psiK, double beta, Range traction, Range compression, Range shear, utils::ProblemType type):
    mesh(mesh), fem(fem), beta(beta), pc(pc), psiK(psiK),
    t(traction), c(compression), s(shear),
    K_e1({0.5*std::pow((t[0]+c[0])/(2*t[0]*c[0]), 2), 
          0.5*std::pow((t[1]+c[1])/(2*t[1]*c[1]), 2)}),
    K_g({1.0/(2*s[0]*s[0]),
         1.0/(2*s[1]*s[1])}),
    K_e2({(t[0]-c[0])/(2*t[0]*c[0]), 
          (t[1]-c[1])/(2*t[1]*c[1])}),
   problem_type(type){

}

double Mechanostat::calculate(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x){

}
double Mechanostat::calculate_with_gradient(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
    std::vector<double> fl(u.size(), 0);

    double result = 0;

    const size_t s_size = this->mesh->elem_info->get_D_dimension();
    const size_t k_size = this->mesh->elem_info->get_k_dimension();
    const size_t dof = this->mesh->elem_info->get_dof_per_node();
    const size_t num_nodes = this->mesh->elem_info->get_nodes_per_element();
    auto x_it = x.cbegin();

    if(mpi_id == 0){
        for(const auto& g:this->mesh->geometries){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                const auto c = e->get_centroid();
                const auto eps_vec = e->get_strain_vector(c, u);

                double rho = 1.0;
                if(g->do_topopt && g->with_void){
                    rho = this->relaxed_rho(*x_it);
                }

                const auto B = e->get_B(c);

                double H_e = 0;
                if(this->problem_type == utils::PROBLEM_TYPE_2D){
                    const StrainVector2D eps = Eigen::Map<const StrainVector2D>(eps_vec.data(), 3);
                    const double eps_lhs1 = rho*this->LHS_2D(0, eps);
                    const double eps_lhs2 = rho*this->LHS_2D(1, eps);
                    const auto deps_lhs1 = rho*this->dLHS_2D(0, eps);
                    const auto deps_lhs2 = rho*this->dLHS_2D(1, eps);
                    const double H1 = Hm(eps_lhs1);
                    const double H2 = Hp(eps_lhs2);
                    const double Hr = H1*eps_lhs1 + H2*eps_lhs2;
                    const auto dH = (dHm(eps_lhs1)*eps_lhs1 + H1)*deps_lhs1 +
                                    (dHp(eps_lhs2)*eps_lhs2 + H2)*deps_lhs2;

                    H_e = Hr;
                    StrainVector2D dHB{0,0,0};
                    for(size_t i = 0; i < s_size; ++i){
                        for(size_t j = 0; j < k_size; ++j){
                            dHB[i] += dH[i]*B[i*k_size + j];
                        }
                    }
                    for(size_t i = 0; i < num_nodes; ++i){
                        for(size_t j = 0; j < dof; ++j){
                            const long pos = e->nodes[i]->u_pos[j];
                            if(pos > -1){
                                fl[pos] += dHB[i*dof + j];
                            }
                        }
                    }
                } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                    const StrainVector3D eps = Eigen::Map<const StrainVector3D>(eps_vec.data(), 6);
                    const double eps_lhs1 = rho*this->LHS_3D(0, eps);
                    const double eps_lhs2 = rho*this->LHS_3D(1, eps);
                    const auto deps_lhs1 = rho*this->dLHS_3D(0, eps);
                    const auto deps_lhs2 = rho*this->dLHS_3D(1, eps);
                    const double H1 = Hm(eps_lhs1);
                    const double H2 = Hp(eps_lhs2);
                    const double Hr = H1*eps_lhs1 + H2*eps_lhs2;
                    const auto dH = (dHm(eps_lhs1)*eps_lhs1 + H1)*deps_lhs1 +
                                    (dHp(eps_lhs2)*eps_lhs2 + H2)*deps_lhs2;

                    H_e = Hr;
                    StrainVector3D dHB{0,0,0,0,0,0};
                    for(size_t i = 0; i < s_size; ++i){
                        for(size_t j = 0; j < k_size; ++j){
                            dHB[i] += dH[i]*B[i*k_size + j];
                        }
                    }
                    for(size_t i = 0; i < num_nodes; ++i){
                        for(size_t j = 0; j < dof; ++j){
                            const long pos = e->nodes[i]->u_pos[j];
                            if(pos > -1){
                                fl[pos] += dHB[i*dof + j];
                            }
                        }
                    }
                }

                result += H_e;

                x_it += num_den;
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
            if(g->do_topopt){
                const size_t num_den = g->number_of_densities_needed();
                std::vector<std::vector<double>> gradD_K(num_den, std::vector<double>(s_size*s_size, 0));
                std::vector<double> D_S(std::vector<double>(s_size*s_size, 0));
                for(const auto& e:g->mesh){
                    const auto c = e->get_centroid();
                    g->materials.get_gradD(x_it, psiK, c, gradD_K);

                    if(g->with_void){
                        double lKu = pc*std::pow(*x_it, pc-1)*e->get_compliance(gradD_K[0], this->mesh->thickness, u, l);
                        const auto eps_vec = e->get_strain_vector(c, u);

                        const double rho = this->relaxed_rho(*x_it);
                        const double drho = this->relaxed_rho_grad(*x_it);

                        double dH_e = 0;
                        if(this->problem_type == utils::PROBLEM_TYPE_2D){
                            const StrainVector2D eps = Eigen::Map<const StrainVector2D>(eps_vec.data(), 3);
                            const double eps_lhs1 = this->LHS_2D(0, eps);
                            const double eps_lhs2 = this->LHS_2D(1, eps);
                            const double H1 = Hm(rho*eps_lhs1);
                            const double H2 = Hp(rho*eps_lhs2);
                            const auto dH = drho*((dHm(rho*eps_lhs1)*rho*eps_lhs1 + H1)*eps_lhs1 +
                                                  (dHp(rho*eps_lhs2)*rho*eps_lhs2 + H2)*eps_lhs2);

                            dH_e = dH;
                        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                            const StrainVector3D eps = Eigen::Map<const StrainVector3D>(eps_vec.data(), 6);
                            const double eps_lhs1 = this->LHS_3D(0, eps);
                            const double eps_lhs2 = this->LHS_3D(1, eps);
                            const double H1 = Hm(rho*eps_lhs1);
                            const double H2 = Hp(rho*eps_lhs2);
                            const auto dH = drho*((dHm(rho*eps_lhs1)*rho*eps_lhs1 + H1)*eps_lhs1 +
                                                  (dHp(rho*eps_lhs2)*rho*eps_lhs2 + H2)*eps_lhs2);

                            dH_e = dH;
                        }

                        *grad_it = dH_e - lKu;
                        //*grad_it = - lKu;

                        const double rho_lKu = std::pow(*x_it, pc);

                        ++x_it;
                        ++grad_it;
                        for(size_t i = 1; i < num_den; ++i){
                            double lKu = rho_lKu*e->get_compliance(gradD_K[i], this->mesh->thickness, u, l);

                            *grad_it = -lKu;

                            ++x_it;
                            ++grad_it;
                        }
                    } else {
                        for(size_t i = 0; i < num_den; ++i){
                            double lKu = e->get_compliance(gradD_K[i], this->mesh->thickness, u, l);

                            *grad_it = -lKu;

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
