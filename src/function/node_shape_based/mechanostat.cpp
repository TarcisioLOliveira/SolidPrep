/*
 *   Copyright (C) 2024 Tarcísio Ladeia de Oliveira.
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

#include <mpi.h>
#include "function/node_shape_based/mechanostat.hpp"
#include "logger.hpp"
#include "math/matrix.hpp"
#include "optimizer.hpp"

namespace function::node_shape_based{

Mechanostat::Mechanostat(const Meshing* const mesh, SolverManager* fem, double beta, Range traction, Range compression, Range shear, utils::ProblemType type):
    mesh(mesh), fem(fem), beta(beta),
    t(traction), c(compression), s(shear),
    K_e1({0.5*std::pow((t[0]+c[0])/(2*t[0]*c[0]), 2), 
          0.5*std::pow((t[1]+c[1])/(2*t[1]*c[1]), 2)}),
    K_g({1.0/(2*s[0]*s[0]),
         1.0/(2*s[1]*s[1])}),
    K_e2({-(t[0]-c[0])/(2*t[0]*c[0]), 
          -(t[1]-c[1])/(2*t[1]*c[1])}),
   problem_type(type){

}

void Mechanostat::initialize_views(Visualization* viz){
    this->shadow_view = viz->add_view("Normalized Strain", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::OTHER);
    this->gradient_view = viz->add_view("Normalized Strain Gradient", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::OTHER);
}

void Mechanostat::initialize(const NodeShapeBasedOptimizer* const op){
    (void)op;
    size_t elem_num = 0;
    for(const auto& g:mesh->geometries){
        elem_num += g->mesh.size();
    }
    this->He.resize(elem_num, std::numeric_limits<double>::quiet_NaN());
    this->gradHe.resize(elem_num, std::numeric_limits<double>::quiet_NaN());
}

double Mechanostat::calculate(const NodeShapeBasedOptimizer* const op, const std::vector<double>& u){
    (void)op;
    (void)u;
    return 0;
}

double Mechanostat::calculate_with_gradient(const NodeShapeBasedOptimizer* const op, const std::vector<double>& u, std::vector<double>& grad){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    std::vector<std::vector<double>> fl(this->mesh->sub_problems->size());
    std::vector<double> l(mesh->max_dofs,0);
    for(size_t i = 0; i < fl.size(); ++i){
        fl[i].resize(mesh->max_dofs, 0);
    }

    auto nodes = op->shape_handler.get_nodes();
    const size_t node_num = this->mesh->elem_info->get_nodes_per_element();
    const size_t kw = this->mesh->elem_info->get_k_dimension();
    const size_t dof = this->mesh->elem_info->get_dof_per_node();
    const size_t num_nodes = this->mesh->elem_info->get_nodes_per_element();
    std::fill(grad.begin(), grad.end(), 0);

    double result = 0;

    if(mpi_id == 0){
        auto He_it = this->He.begin();
        for(const auto& g:this->mesh->geometries){
            for(const auto& e:g->mesh){
                const auto c = e->get_centroid();

                const auto B = e->get_B(c);

                double H_e = 0;
                for(size_t it = 0; it < fl.size(); ++it){
                    const auto& ui = fem->sub_u[it];
                    const math::Vector eps(1e6*e->get_strain_vector(c, ui));
                    if(this->problem_type == utils::PROBLEM_TYPE_2D){
                        const double eps_lhs1 = this->LHS_2D(0, eps);
                        const double eps_lhs2 = this->LHS_2D(1, eps);
                        const StrainVector2D deps_lhs1 = this->dLHS_2D(0, eps);
                        const StrainVector2D deps_lhs2 = this->dLHS_2D(1, eps);
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
                        const double eps_lhs1 = this->LHS_3D(0, eps);
                        const double eps_lhs2 = this->LHS_3D(1, eps);
                        const StrainVector3D deps_lhs1 = this->dLHS_3D(0, eps);
                        const StrainVector3D deps_lhs2 = this->dLHS_3D(1, eps);
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

                ++He_it;
            }
        }
        logger::quick_log("Calculating adjoint problem...{");
    }
    this->shadow_view->update_view(this->He);
    this->fem->calculate_displacements_adjoint(this->mesh, fl, l);
    if(mpi_id == 0){
        logger::quick_log("} Done.");

        logger::quick_log("Calculating stress gradient...");

        math::Vector ue(kw);
        math::Vector le(kw);
        math::Matrix D;

        for(size_t i = 0; i < nodes.size(); ++i){
            const auto& shn = nodes[i];
            for(const auto& e:shn.elements){
                const auto& g = this->mesh->elem_geom_mapping.at(e.e);
                const auto c = e.e->get_centroid();
                D = g->materials.get_D(e.e, c);

                const math::Vector eps(1e6*e.e->get_strain_vector(c, u));

                for(size_t n = 0; n < node_num; ++n){
                    for(size_t j = 0; j < dof; ++j){
                        ue[n*dof + j] = u[e.e->nodes[n]->u_pos[j]];
                        le[n*dof + j] = l[e.e->nodes[n]->u_pos[j]];
                    }
                }

                math::Vector dH_e;
                if(this->problem_type == utils::PROBLEM_TYPE_2D){
                    const double eps_lhs1 = this->LHS_2D(0, eps);
                    const double eps_lhs2 = this->LHS_2D(1, eps);
                    const StrainVector2D deps_lhs1 = this->dLHS_2D(0, eps);
                    const StrainVector2D deps_lhs2 = this->dLHS_2D(1, eps);
                    const double H1 = Hm(eps_lhs1);
                    const double H2 = Hp(eps_lhs2);
                    StrainVector2D dH = (dHm(eps_lhs1)*eps_lhs1 + H1)*deps_lhs1 +
                                        (dHp(eps_lhs2)*eps_lhs2 + H2)*deps_lhs2;

                    dH_e = std::move(dH);
                } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                    const double eps_lhs1 = this->LHS_3D(0, eps);
                    const double eps_lhs2 = this->LHS_3D(1, eps);
                    const StrainVector3D deps_lhs1 = this->dLHS_3D(0, eps);
                    const StrainVector3D deps_lhs2 = this->dLHS_3D(1, eps);
                    const double H1 = Hm(eps_lhs1);
                    const double H2 = Hp(eps_lhs2);
                    StrainVector3D dH = (dHm(eps_lhs1)*eps_lhs1 + H1)*deps_lhs1 +
                                        (dHp(eps_lhs2)*eps_lhs2 + H2)*deps_lhs2;

                    dH_e = std::move(dH);
                }

                for(size_t j = 0; j < dof; ++j){
                    const double ldKu = le.T()*e.e->get_dk_sh(D, this->mesh->thickness, e.node_num, j)*ue;

                    auto dB(e.e->get_dB_sh(c, e.node_num, j));

                    const double mult = dH_e.T()*dB*ue;

                    const double dS = e.e->von_Mises_derivative_sh(D, mult, c, u, e.node_num, j);

                    grad[i*dof + j] += dS - ldKu;
                }
            }
        }
        logger::quick_log("Done.");
    }

    return result;
}

}
