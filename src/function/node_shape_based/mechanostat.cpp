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
#include "project_data.hpp"

namespace function::node_shape_based{

Mechanostat::Mechanostat(const projspec::DataMap& data):
    mesh(data.proj->topopt_mesher.get()),
    fem(data.proj->topopt_fea.get()),
    beta(data.get_double("beta")),
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
    //this->gradient_view = viz->add_view("Normalized Strain Gradient", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::OTHER);
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

    #define CONSIDER_SPREAD 1

    #ifdef CONSIDER_SPREAD
        const auto& clusters = op->shape_handler.get_clusters();
        std::vector<std::vector<double>> shape_fl;
        shape_fl.resize(clusters.size());
        for(size_t i = 0; i < clusters.size(); ++i){
            shape_fl[i].resize(clusters[i].matrix_width, 0);
        }
    #endif

    if(mpi_id == 0){
        auto He_it = this->He.begin();
        //math::Vector ue(kw);
        for(const auto& g:this->mesh->geometries){
            //const auto& cluster = op->shape_handler.get_cluster(g);
            //auto& sfl = shape_fl[cluster->id];
            for(const auto& e:g->mesh){
                const auto c = e->get_centroid();

                const auto B = e->get_B(c);

                //const auto D = g->materials.get_D(e.get(), c);
                //const auto D0 = g->materials.get_D_orig(e.get(), c);

                double H_e = 0;
                math::Vector dH_e;
                const size_t it = 0;

                //const auto& ui = fem->get_sub_u(it);
                const math::Vector eps(1e6*e->get_strain_vector(c, u));
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

                    *He_it = this->c[1]*eps_lhs2;
                    H_e = Hr;
                    dH_e = std::move(dH);
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

                    *He_it = this->c[1]*eps_lhs2;
                    H_e = Hr;
                    dH_e = std::move(dH);
                }
                math::Vector dHB(1e6*(dH_e.T()*B));
                for(size_t i = 0; i < num_nodes; ++i){
                    for(size_t j = 0; j < dof; ++j){
                        const long pos = e->nodes[i]->u_pos[j];
                        fl[it][pos] += dHB[i*dof + j];
                    }
                }

                //for(size_t n = 0; n < node_num; ++n){
                //    for(size_t j = 0; j < dof; ++j){
                //        ue[n*dof + j] = u[e->nodes[n]->u_pos[j]];
                //    }
                //}
                //for(size_t i = 0; i < node_num; ++i){
                //    const auto id = cluster->id_mapping.at(e->nodes[i]->id);
                //    if(id >= 0){
                //        for(size_t j = 0; j < dof; ++j){
                //            const auto dB(e->get_dB_sh(c, i, j));
                //            const auto dHxi(1e6*(dH_e.T()*dB*ue));
                //            sfl[id*dof + j] += dHxi;
                //        }
                //    }
                //}

                //*He_it = H_e;
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

        #pragma omp parallel
        {
            math::Vector ue(kw);
            math::Vector le(kw);
            //math::Matrix D;
            //math::Matrix D0;
            //math::Matrix dD;
            //math::Matrix dD0;
            #pragma omp for
            for(size_t i = 0; i < nodes.size(); ++i){
                const auto& shn = nodes[i];
                for(const auto& e:shn.elements){
                    const auto& g = this->mesh->elem_geom_mapping.at(e.e);
                    const auto c = e.e->get_centroid();
                    //D = g->materials.get_D(e.e, c);
                    //D0 = g->materials.get_D_orig(e.e, c);

                    const math::Vector eps(1e6*e.e->get_strain_vector(c, u));

                    for(size_t n = 0; n < node_num; ++n){
                        for(size_t j = 0; j < dof; ++j){
                            ue[n*dof + j] = u[e.e->nodes[n]->u_pos[j]];
                            le[n*dof + j] = l[e.e->nodes[n]->u_pos[j]];
                        }
                    }

                    for(size_t j = 0; j < dof; ++j){
                        //dD = g->materials.get_D_dsh(e.e, c, e.node_num, j);
                        //dD0 = g->materials.get_D_orig_dsh(e.e, c, e.node_num, j);

                        const double ldKu = le.T()*e.e->get_dk_sh(g->materials, this->mesh->thickness, e.node_num, j)*ue;

                        const auto dB(e.e->get_dB_sh(c, e.node_num, j));
                        const auto deps(1e6*dB*ue);

                        double dH_e = 0;
                        if(this->problem_type == utils::PROBLEM_TYPE_2D){
                            const double eps_lhs1 = this->LHS_2D(0, eps);
                            const double eps_lhs2 = this->LHS_2D(1, eps);
                            const auto deps_lhs1 = this->dLHS_2D(0, eps).T()*deps;
                            const auto deps_lhs2 = this->dLHS_2D(1, eps).T()*deps;
                            const double H1 = Hm(eps_lhs1);
                            const double H2 = Hp(eps_lhs2);
                            const auto dH = (dHm(eps_lhs1)*eps_lhs1 + H1)*deps_lhs1 +
                                            (dHp(eps_lhs2)*eps_lhs2 + H2)*deps_lhs2;

                            dH_e = dH;
                        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                            const double eps_lhs1 = this->LHS_3D(0, eps);
                            const double eps_lhs2 = this->LHS_3D(1, eps);
                            //const auto deps_lhs1 = this->dUd_3D_full(0, eps, deps, D0, D, dD0, dD);
                            //const auto deps_lhs2 = this->dUd_3D_full(1, eps, deps, D0, D, dD0, dD);
                            const auto deps_lhs1 = this->dLHS_3D(0, eps).T()*deps;
                            const auto deps_lhs2 = this->dLHS_3D(1, eps).T()*deps;
                            const double H1 = Hm(eps_lhs1);
                            const double H2 = Hp(eps_lhs2);
                            const auto dH = (dHm(eps_lhs1)*eps_lhs1 + H1)*deps_lhs1 +
                                            (dHp(eps_lhs2)*eps_lhs2 + H2)*deps_lhs2;

                            dH_e = dH;
                        }

                        grad[i*dof + j] += dH_e - ldKu;
                    }
                }
                // Remove normal component from external boundary-adjacent nodes
                // double dot = 0;
                // for(size_t j = 0; j < dof; ++j){
                //     dot += grad[i*dof + j]*shn.bound_norm[j];
                // }
                // for(size_t j = 0; j < dof; ++j){
                //     grad[i*dof + j] -= dot*shn.bound_norm[j];
                // }
            }
        }
        #ifdef CONSIDER_SPREAD
            // dK
            for(const auto& g:this->mesh->geometries){
                const auto& cluster = op->shape_handler.get_cluster(g);
                auto& sfl = shape_fl[cluster->id];
                std::vector<double> ldku(g->mesh.size()*dof*node_num);
                logger::quick_log("geom", g->id);
                #pragma omp parallel
                {
                    math::Vector ue(kw);
                    math::Vector le(kw);
                    #pragma omp for
                    for(size_t ei = 0; ei < g->mesh.size(); ++ei){
                        const auto& e = g->mesh[ei];
                        for(size_t ni = 0; ni < node_num; ++ni){
                            for(size_t nj = 0; nj < dof; ++nj){
                                ue[ni*dof + nj] = u[e->nodes[ni]->u_pos[nj]];
                                le[ni*dof + nj] = l[e->nodes[ni]->u_pos[nj]];
                            }
                        }
                        const auto dk = e->get_dk_sh_full(g->materials, this->mesh->thickness);
                        for(size_t i = 0; i < node_num; ++i){
                            for(size_t j = 0; j < dof; ++j){
                                ldku[ei*kw + i*dof + j] = le.T()*dk[i*dof + j]*ue;
                            }
                        }
                    }
                }
                math::Vector ue(kw);
                for(size_t ei = 0; ei < g->mesh.size(); ++ei){
                    const auto& e = g->mesh[ei];
                    for(size_t ni = 0; ni < node_num; ++ni){
                        for(size_t nj = 0; nj < dof; ++nj){
                            ue[ni*dof + nj] = u[e->nodes[ni]->u_pos[nj]];
                        }
                    }
                    for(size_t i = 0; i < node_num; ++i){
                        const auto id = cluster->id_mapping.at(e->nodes[i]->id);
                        if(id >= 0){
                            for(size_t j = 0; j < dof; ++j){
                                const auto c = e->get_centroid();
                                const math::Vector eps(1e6*e->get_strain_vector(c, u));


                                const auto dB(e->get_dB_sh(c, i, j));
                                const auto deps(1e6*dB*ue);

                                double dH_e = 0;
                                if(this->problem_type == utils::PROBLEM_TYPE_2D){
                                    const double eps_lhs1 = this->LHS_2D(0, eps);
                                    const double eps_lhs2 = this->LHS_2D(1, eps);
                                    const auto deps_lhs1 = this->dLHS_2D(0, eps).T()*deps;
                                    const auto deps_lhs2 = this->dLHS_2D(1, eps).T()*deps;
                                    const double H1 = Hm(eps_lhs1);
                                    const double H2 = Hp(eps_lhs2);
                                    const auto dH = (dHm(eps_lhs1)*eps_lhs1 + H1)*deps_lhs1 +
                                                    (dHp(eps_lhs2)*eps_lhs2 + H2)*deps_lhs2;

                                    dH_e = dH;
                                } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                                    const double eps_lhs1 = this->LHS_3D(0, eps);
                                    const double eps_lhs2 = this->LHS_3D(1, eps);
                                    const auto deps_lhs1 = this->dLHS_3D(0, eps).T()*deps;
                                    const auto deps_lhs2 = this->dLHS_3D(1, eps).T()*deps;
                                    const double H1 = Hm(eps_lhs1);
                                    const double H2 = Hp(eps_lhs2);
                                    const auto dH = (dHm(eps_lhs1)*eps_lhs1 + H1)*deps_lhs1 +
                                                    (dHp(eps_lhs2)*eps_lhs2 + H2)*deps_lhs2;

                                    dH_e = dH;
                                }

                                sfl[id*dof + j] += dH_e - ldku[ei*kw + i*dof + j];
                            }
                        }
                    }
                }
            }
        #endif
        op->contact_gradient(mesh, fem, u, l, grad);
        #ifdef CONSIDER_SPREAD
            for(size_t i = 0; i < clusters.size(); ++i){
                clusters[i].solver->solve_adjoint(shape_fl[i]);
            }
            op->shape_handler.adjoint_grad_append(grad, shape_fl);
        #endif
        logger::quick_log("Done.");
    }

    return result;
}

using namespace projspec;
const bool Mechanostat::reg = Factory<NodeShapeBasedFunction>::add(
    [](const DataMap& data){
        return std::make_unique<Mechanostat>(data);
    },
    ObjectRequirements{
        "mechanostat",
        {
            DataEntry{.name = "beta", .type = TYPE_DOUBLE, .required = true},
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
