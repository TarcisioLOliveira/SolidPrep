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

#include "global_stiffness_matrix.hpp"
#include "logger.hpp"
#include "math/matrix.hpp"
#include "math/slicer.hpp"
#include "utils.hpp"
#include <limits>

GlobalStiffnessMatrix::GlobalStiffnessMatrix(const FiniteElement::ContactData& data):
    HK(data.HK),
    HMU(data.HMU),
    EPS_PENALTY(data.EPS_DISPL_SIMPLE), EPS_DISPL(data.EPS_DISPL_SIMPLE), LAG_DISPL_SIMPLE(data.EPS_DISPL_SIMPLE)
    {

    }

void GlobalStiffnessMatrix::generate_base(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type){
    size_t D_offset = 0;
    (void)lambda;
    (void)l_num;

    /*
    size_t L = u_size;
    if(type != FiniteElement::ContactType::RIGID){
        L += l_num;
        //if(type == FiniteElement::ContactType::FRICTIONLESS_DISPL_SIMPLE){
        //    L += 1;
        //}
    }
    */
    for(auto& g : mesh->geometries){
        if((topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous()){
            this->add_geometry(mesh, node_positions, g, D_offset, D_cache);
            D_offset += g->mesh.size();
        } else {
            this->add_geometry(mesh, node_positions, g);
        }
    }
    if(mesh->springs->size() > 0){
        this->add_springs(mesh, node_positions);
    }
    // Reserve space but do not apply matrices
    // Only if using coordinate field?
    if(this->first_time){
        this->first_time = false;

        const size_t dof      = mesh->elem_info->get_dof_per_node();
        const size_t node_num = mesh->elem_info->get_nodes_per_element();
        const size_t bnum     = mesh->elem_info->get_boundary_nodes_per_element();
        std::vector<long> u_pos(dof*node_num);

        const math::Matrix k(u_pos.size(), u_pos.size(), 1);

        for(auto& g : mesh->geometries){
            for(auto& e : g->mesh){
                for(size_t i = 0; i < node_num; ++i){
                    const auto& n = e->nodes[i];
                    for(size_t j = 0; j < dof; ++j){
                        const size_t p = n->id*dof + j;
                        u_pos[i*dof + j] = node_positions[p];
                    }
                }
                this->reserve_element_matrix(k, u_pos);
            }
        }
        if(type == FiniteElement::ContactType::FRICTIONLESS_DISPL_LOG){

            this->lag_pos_cache.resize(mesh->paired_boundary.size()*node_num, -1);
            size_t cur_id = 0;
            for(const auto& e:mesh->paired_boundary){
                for(size_t i = 0; i < node_num; ++i){
                    for(size_t j = 0; j < bnum; ++j){
                        if(e.b1->parent->nodes[i] == e.b1->nodes[j]){
                            const auto id = mesh->lag_node_map.at(e.b1->parent->nodes[i]->id);
                            this->lag_pos_cache[cur_id*node_num + i] = id;
                            break;
                        }
                    }
                }
                ++cur_id;
            }

            //this->add_frictionless_log(mesh, node_positions, u_ext, lambda, true);
        } else if(type == FiniteElement::ContactType::FRICTIONLESS_DISPL_SIMPLE){
            const std::vector<double> empty_lambda(lambda.size(), 0);
            this->add_frictionless_simple(mesh, node_positions, u_ext, empty_lambda, true);
        }
        // Internal loads 2
        if(mesh->internal_loads2 != nullptr && mesh->internal_loads2->size() > 0){
            const size_t dof      = mesh->elem_info->get_dof_per_node();
            const size_t node_num = mesh->elem_info->get_nodes_per_element();
            const size_t s_size   = mesh->elem_info->get_D_dimension();
            //const size_t il2_num  = mesh->internal_loads2->size();

            std::vector<long> u_pos(dof*node_num);
            std::vector<long> f_pos(dof*node_num);
            std::vector<long> s_pos(s_size);
            std::vector<long> s_pos2(s_size);
            math::Matrix M(s_size, u_pos.size(), 1);
            math::Matrix F(u_pos.size(), u_pos.size(), 1);
            for(const auto& il:*(mesh->internal_loads2)){
                for(size_t i = 0; i < 6; ++i){
                    s_pos[i] = il.constraint_offset + i;
                    s_pos2[i] = il.constraint_offset + s_size + i;
                }
                for(const auto& b:il.submesh){
                    for(size_t i = 0; i < node_num; ++i){
                        const auto& n = b->parent->nodes[i];
                        const size_t f_id = il.F_id_map.at(n->id);
                        for(size_t j = 0; j < dof; ++j){
                            const size_t p = n->id*dof + j;
                            u_pos[i*dof + j] = node_positions[p];
                        }
                        for(size_t j = 0; j < dof; ++j){
                            const size_t fp = f_id*dof + j;
                            f_pos[i*dof + j] = il.il2_offset + fp;
                        }
                    }
                    this->reserve_element_matrix(F, f_pos);
                    //this->reserve_block_symmetric(F, u_pos, f_pos);
                    //this->reserve_block_symmetric(M, s_pos, u_pos);
                    this->reserve_block_symmetric(M, s_pos2, u_pos);
                }
            }
        }
    }
    // Internal loads 2
    if(mesh->internal_loads2 != nullptr && mesh->internal_loads2->size() > 0){
        const size_t dof      = mesh->elem_info->get_dof_per_node();
        const size_t node_num = mesh->elem_info->get_nodes_per_element();
        const size_t s_size   = mesh->elem_info->get_D_dimension();
        //const auto prob_type = mesh->elem_info->get_problem_type();
        //const size_t dim = (prob_type == utils::PROBLEM_TYPE_2D) ? 2 : 3;

        std::vector<long> u_pos(dof*node_num);
        std::vector<long> f_pos(dof*node_num);
        std::vector<long> s_pos(s_size);
        std::vector<long> s_pos2(s_size);
        for(const auto& il:*(mesh->internal_loads2)){
            for(size_t i = 0; i < 6; ++i){
                s_pos[i] = il.constraint_offset + i;
                s_pos2[i] = il.constraint_offset + s_size + i;
            }
            for(const auto& b:il.submesh){
                for(size_t i = 0; i < node_num; ++i){
                    const auto& n = b->parent->nodes[i];
                    const size_t f_id = il.F_id_map.at(n->id);
                    for(size_t j = 0; j < dof; ++j){
                        const size_t p = n->id*dof + j;
                        u_pos[i*dof + j] = node_positions[p];
                    }
                    if(f_id >= 0){
                        for(size_t j = 0; j < dof; ++j){
                            const size_t fp = f_id*dof + j;
                            f_pos[i*dof + j] = il.il2_offset + fp;
                        }
                    } else {
                        for(size_t j = 0; j < dof; ++j){
                            f_pos[i*dof + j] = -1;
                        }
                    }
                }
                //const auto H = il.get_hess(b);
                //this->insert_element_matrix(H, u_pos);
                const auto K = il.get_K_mat(b);
                const auto F = il.get_F_mat(b);
                const auto M = il.get_constr(b);
                const auto M2 = il.get_constr2(b);
                this->insert_element_matrix(K, f_pos);
                this->insert_block_symmetric(F, u_pos, f_pos);
                this->insert_block_symmetric(M, s_pos, f_pos);
                this->insert_block_symmetric(M2, s_pos2, u_pos);
            }
        }
    }
    if(topopt){
        for(size_t i = 0; i < u_size; ++i){
            this->add_to_matrix(i, i, this->K_MIN);
        }
    }
}

void GlobalStiffnessMatrix::add_geometry(const Meshing* const mesh, const std::vector<long>& node_positions, const Geometry* const g){
    const double t = mesh->thickness;
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    std::vector<long> u_pos(dof*node_num);
    const auto D = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
    for(auto& e : g->mesh){
        for(size_t i = 0; i < node_num; ++i){
            const auto& n = e->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                const size_t p = n->id*dof + j;
                u_pos[i*dof + j] = node_positions[p];
            }
        }
        math::Matrix k(e->get_k(D, t));
        this->insert_element_matrix(k, u_pos);
    }
}

void GlobalStiffnessMatrix::add_geometry(const Meshing* const mesh, const std::vector<long>& node_positions, const Geometry* const g, const size_t D_offset, const std::vector<math::Matrix>& D_cache){
    const double t = mesh->thickness;
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    std::vector<long> u_pos(dof*node_num);

    size_t it = 0;
    for(const auto& e:g->mesh){
        for(size_t i = 0; i < node_num; ++i){
            const auto& n = e->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                const size_t p = n->id*dof + j;
                u_pos[i*dof + j] = node_positions[p];
            }
        }
        const auto& D = D_cache[D_offset + it];
        const math::Matrix k(e->get_k(D, t));
        this->insert_element_matrix(k, u_pos);
        ++it;
    }
}

void GlobalStiffnessMatrix::add_springs(const Meshing * const mesh, const std::vector<long>& node_positions){
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnode_num = mesh->elem_info->get_boundary_nodes_per_element();
    std::vector<long> u_pos(dof*node_num);

    std::vector<gp_Pnt> points(bnode_num);
    for(size_t it = 0; it < mesh->springs->size(); ++it){
        for(const auto& b : mesh->springs->at(it).submesh){
            const auto c = b->get_centroid(bnode_num);
            const auto& K = mesh->springs->at(it).get_K(b->parent, c);
            const auto& e = b->parent;
            for(size_t i = 0; i < bnode_num; ++i){
                points[i] = b->nodes[i]->point;
            }
            for(size_t i = 0; i < node_num; ++i){
                const auto& n = e->nodes[i];
                for(size_t j = 0; j < dof; ++j){
                    const size_t p = n->id*dof + j;
                    u_pos[i*dof + j] = node_positions[p];
                }
            }
            const math::Matrix R(e->get_R(K, mesh->thickness, points));
            this->insert_element_matrix(R, u_pos);
        }
    }
}

void GlobalStiffnessMatrix::add_nonlinear_bcs(const Meshing * const mesh, const std::vector<long>& node_positions, const std::vector<double>& u_ext){
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnode_num = mesh->elem_info->get_boundary_nodes_per_element();
    std::vector<long> u_pos(dof*node_num);

    std::vector<gp_Pnt> points(bnode_num);

    // Contact springs
    for(size_t it = 0; it < mesh->contact_springs->size(); ++it){
        const auto& cs = mesh->contact_springs->at(it);
        const math::Vector normal({
            cs.normal.X(),
            cs.normal.Y(),
            cs.normal.Z()
        });

        for(const auto& b : mesh->contact_springs->at(it).submesh){
            const auto c = b->get_centroid(bnode_num);
            const auto& K = mesh->contact_springs->at(it).get_K(b->parent, c);
            const auto& e = b->parent;
            for(size_t i = 0; i < bnode_num; ++i){
                points[i] = b->nodes[i]->point;
            }
            for(size_t i = 0; i < node_num; ++i){
                const auto& n = e->nodes[i];
                for(size_t j = 0; j < dof; ++j){
                    const size_t p = n->id*dof + j;
                    u_pos[i*dof + j] = node_positions[p];
                }
            }
            const math::Matrix R(e->get_R_contact(K, mesh->thickness, points, normal, u_ext));
            this->insert_element_matrix(R, u_pos);
        }
    }

    // // Internal loads 2
    // if(mesh->internal_loads2->size() > 0){
    //     const size_t dof      = mesh->elem_info->get_dof_per_node();
    //     const size_t node_num = mesh->elem_info->get_nodes_per_element();

    //     std::vector<long> u_pos(dof*node_num + 6);
    //     const size_t il_offset = mesh->load_vector[0].size() - mesh->internal_loads2->size()*6;
    //     size_t il_id = 0;
    //     for(const auto& il:*(mesh->internal_loads2)){
    //         for(size_t i = 0; i < 6; ++i){
    //             u_pos[dof*node_num + i] = il_offset + 6*il_id + i;
    //         }
    //         for(const auto& b:il.submesh){
    //             for(size_t i = 0; i < node_num; ++i){
    //                 const auto& n = b->parent->nodes[i];
    //                 for(size_t j = 0; j < dof; ++j){
    //                     const size_t p = n->id*dof + j;
    //                     u_pos[i*dof + j] = node_positions[p];
    //                 }
    //             }
    //             const auto M = il.get_hess(b);
    //             this->insert_element_matrix(M, u_pos);
    //         }
    //         ++il_id;
    //     }
    // }
}

void GlobalStiffnessMatrix::append_Ku_nl_bcs(const Meshing* const mesh, const std::vector<long>& node_positions, const std::vector<double>& u_ext, std::vector<double>& Ku) const{
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnode_num = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();

    std::vector<gp_Pnt> points(node_num);

    // Contact springs
    for(size_t it = 0; it < mesh->contact_springs->size(); ++it){
        const auto& cs = mesh->contact_springs->at(it);
        const math::Vector normal({
            cs.normal.X(),
            cs.normal.Y(),
            cs.normal.Z()
        });

        for(const auto& b : mesh->contact_springs->at(it).submesh){
            const auto c = b->get_centroid(bnode_num);
            const auto& K = mesh->contact_springs->at(it).get_K(b->parent, c);
            const auto& e = b->parent;
            for(size_t i = 0; i < bnode_num; ++i){
                points[i] = b->nodes[i]->point;
            }
            const math::Vector Kue(e->get_Ku_cs(K, mesh->thickness, points, normal, u_ext, u_ext));
            for(size_t i = 0; i < node_num; ++i){
                for(size_t j = 0; j < dof; ++j){
                    const size_t p = e->nodes[i]->id*dof + j;
                    const auto pos = node_positions[p];
                    if(pos >= 0){
                        Ku[pos] += Kue[i*dof + j];
                    }
                }
            }
        }
    }

    // // Internal loads 2
    // if(mesh->internal_loads2->size() > 0){
    //     //const size_t dof      = mesh->elem_info->get_dof_per_node();
    //     //const size_t node_num = mesh->elem_info->get_nodes_per_element();

    //     const size_t il_offset = mesh->load_vector[0].size() - mesh->internal_loads2->size()*6;
    //     size_t il_id = 0;
    //     for(const auto& il:*(mesh->internal_loads2)){
    //         const auto mult = il.calculate_multipliers(u_ext);
    //         //for(const auto& b:il.submesh){
    //         //    const math::Vector Kue(il.get_grad(b, mult));
    //         //    for(size_t i = 0; i < node_num; ++i){
    //         //        for(size_t j = 0; j < dof; ++j){
    //         //            const size_t p = b->parent->nodes[i]->id*dof + j;
    //         //            const auto pos = node_positions[p];
    //         //            if(pos >= 0){
    //         //                Ku[pos] += Kue[i*dof + j];
    //         //            }
    //         //        }
    //         //    }
    //         //    for(size_t i = 0; i < 6; ++i){
    //         //        Ku[il_offset + 6*il_id + i] += Kue[dof*node_num + i];
    //         //    }
    //         //}
    //         for(size_t i = 0; i < 6; ++i){
    //             Ku[il_offset + 6*il_id + i] += mult[i];
    //         }
    //         ++il_id;
    //     }
    // }
}

void GlobalStiffnessMatrix::append_dKu_nl_bcs(const Meshing* const mesh, const std::vector<long>& node_positions, const std::vector<double>& u_ext, const std::vector<double>& du_ext, const std::vector<double>& du, std::vector<double>& Ku) const{
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnode_num = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();

    std::vector<gp_Pnt> points(node_num);

    // Contact springs
    for(size_t it = 0; it < mesh->contact_springs->size(); ++it){
        const auto& cs = mesh->contact_springs->at(it);
        const math::Vector normal({
            cs.normal.X(),
            cs.normal.Y(),
            cs.normal.Z()
        });

        for(const auto& b : mesh->contact_springs->at(it).submesh){
            const auto c = b->get_centroid(bnode_num);
            const auto& K = mesh->contact_springs->at(it).get_K(b->parent, c);
            const auto& e = b->parent;
            for(size_t i = 0; i < bnode_num; ++i){
                points[i] = b->nodes[i]->point;
            }
            const math::Vector Kue(e->get_Ku_cs(K, mesh->thickness, points, normal, u_ext, du_ext));
            for(size_t i = 0; i < node_num; ++i){
                for(size_t j = 0; j < dof; ++j){
                    const size_t p = e->nodes[i]->id*dof + j;
                    const auto pos = node_positions[p];
                    if(pos >= 0){
                        Ku[pos] += Kue[i*dof + j];
                    }
                }
            }
        }
    }

    // // Internal loads 2
    // if(mesh->internal_loads2->size() > 0){
    //     const size_t dof      = mesh->elem_info->get_dof_per_node();
    //     const size_t node_num = mesh->elem_info->get_nodes_per_element();

    //     const size_t il_offset = mesh->load_vector[0].size() - mesh->internal_loads2->size()*6;
    //     size_t il_id = 0;
    //     math::Vector duv_ext(dof*node_num + 6);
    //     for(const auto& il:*(mesh->internal_loads2)){
    //         // for(size_t i = 0; i < 6; ++i){
    //         //     duv_ext[dof*node_num + i] = du[il_offset + il_id*6 + i];
    //         // }
    //         // for(const auto& b:il.submesh){
    //         //     math::VectorSliceGeneralView duv(du_ext, math::slicer::from_node_upos(b->parent->nodes, node_num, dof));
    //         //     for(size_t i = 0; i < dof*node_num; ++i){
    //         //         duv_ext[i] = duv[i];
    //         //     }
    //         //     const math::Vector Kue(il.get_hess(b)*duv_ext);
    //         //     for(size_t i = 0; i < node_num; ++i){
    //         //         for(size_t j = 0; j < dof; ++j){
    //         //             const size_t p = b->parent->nodes[i]->id*dof + j;
    //         //             const auto pos = node_positions[p];
    //         //             if(pos >= 0){
    //         //                 Ku[pos] += Kue[i*dof + j];
    //         //             }
    //         //         }
    //         //     }
    //         //     for(size_t i = 0; i < 6; ++i){
    //         //         Ku[il_offset + 6*il_id + i] += Kue[dof*node_num + i];
    //         //     }
    //         // }
    //         (void) il;
    //         for(size_t i = 0; i < 6; ++i){
    //             Ku[il_offset + 6*il_id + i] += du[il_offset + il_id*6 + i];
    //         }
    //         ++il_id;
    //     }
    // }
}

void GlobalStiffnessMatrix::append_il2_f(const Meshing* const mesh, std::vector<double>& f) const{
    // Internal loads 2
    //if(mesh->internal_loads2->size() > 0){
    //    //const size_t dof      = mesh->elem_info->get_dof_per_node();
    //    //const size_t node_num = mesh->elem_info->get_nodes_per_element();

    //    size_t il_id = 0;
    //    for(const auto& il:*(mesh->internal_loads2)){
    //        const auto mult = il.get_loads();
    //        for(size_t i = 0; i < 6; ++i){
    //            f[mesh->il2_offset+ 6*il_id + i] += mult[i];
    //        }
    //        ++il_id;
    //    }
    //}
}

void GlobalStiffnessMatrix::calculate_dimensions(const Meshing* const mesh, const std::vector<long>& node_positions, const size_t matrix_width){
    const size_t k_dim    = mesh->elem_info->get_k_dimension();
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();

    W = matrix_width;
    N = k_dim;

    for(auto& g:mesh->geometries){
        for(auto& e:g->mesh){
            size_t min_i = 0;
            size_t max_i = 0;
            long min = std::numeric_limits<long>::max();
            long max = -1;
            std::vector<long> pos;
            pos.reserve(k_dim);
            for(size_t i = 0; i < node_num; ++i){
                const auto& n = e->nodes[i];
                for(size_t j = 0; j < dof; ++j){
                    const size_t p = n->id*dof + j;
                    pos.push_back(node_positions[p]);
                }
            }
            for(size_t i = 0; i < pos.size(); ++i){
                if(pos[i] > -1){
                    if(pos[i] < min){
                        min = pos[i];
                        min_i = i;
                    }
                }
                if(pos[i] > max){
                    max = pos[i];
                    max_i = i;
                }
            }
            size_t N_candidate = pos[max_i] - pos[min_i] + 1;
            if(N_candidate > N){
                N = N_candidate;
            }
        }
    }
}

void GlobalStiffnessMatrix::add_contacts(const Meshing * const mesh, const std::vector<long>& node_positions, const std::vector<double>& u_ext){
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();

    const size_t kw = mesh->elem_info->get_k_dimension();
    std::vector<gp_Pnt> points(node_num);
    std::vector<long> u_pos(2*kw);
    size_t added = 0;
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < bnum; ++i){
            points[i] = e.b1->nodes[i]->point;
        }
        for(size_t i = 0; i < node_num; ++i){
            const auto n1 = e.b1->parent->nodes[i];
            const auto n2 = e.b2->parent->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                u_pos[dof*i + j] = node_positions[n1->u_pos[j]];
                u_pos[dof*(node_num + i) + j] = node_positions[n2->u_pos[j]];
            }
        }
        const auto MM(EPS_PENALTY*e.b1->parent->get_MnMn(e.b2->parent, u_ext, points, -e.b1->normal));
        for(size_t i = 0; i < 2*kw; ++i){
            if(std::abs(MM(i, i)) > 1e-14){
                ++added;
                break;
            }
        }
        this->insert_element_matrix(MM, u_pos);
    }
    this->final_flush_matrix();
    logger::quick_log("added", added);
}
void GlobalStiffnessMatrix::append_Ku_penalty(const Meshing* const mesh, const std::vector<long>& node_positions, const std::vector<double>& u_ext, const std::vector<double>& lambda, std::vector<double>& Ku) const{
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();

    const size_t kw = mesh->elem_info->get_k_dimension();
    std::vector<gp_Pnt> points(node_num);
    std::vector<long> u_pos(2*kw);
    std::vector<long> l_pos(bnum);
    math::Vector u1(kw);
    math::Vector u2(kw);
    math::Vector uv(2*kw);
    math::Vector l_e(bnum);
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < bnum; ++i){
            points[i] = e.b1->nodes[i]->point;
        }
        for(size_t i = 0; i < bnum; ++i){
            l_pos[i] = mesh->lag_node_map.at(e.b1->nodes[i]->id);
            l_e[i] = lambda[l_pos[i]];
        }
        for(size_t i = 0; i < node_num; ++i){
            const auto n1 = e.b1->parent->nodes[i];
            const auto n2 = e.b2->parent->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                u_pos[dof*i + j] = node_positions[n1->u_pos[j]];
                u_pos[dof*(node_num + i) + j] = node_positions[n2->u_pos[j]];
                uv[dof*i + j] = u_ext[n1->u_pos[j]];
                uv[dof*(node_num + i) + j] = u_ext[n2->u_pos[j]];
                u1[dof*i + j] = u_ext[n1->u_pos[j]];
                u2[dof*i + j] = u_ext[n2->u_pos[j]];
            }
        }
        const auto MM(EPS_PENALTY*e.b1->parent->get_MnMn(e.b2->parent, u_ext, points, -e.b1->normal));
        const auto uL(e.elem->fl2_uL(u1, u2));
        const auto Kue = MM*uv;
        for(size_t i = 0; i < 2*kw; ++i){
            if(u_pos[i] >= 0){
                Ku[u_pos[i]] += Kue[i];
            }
        }
    }
}
void GlobalStiffnessMatrix::append_dKu_penalty(const Meshing* const mesh, const std::vector<long>& node_positions, const std::vector<double>& u_ext, const std::vector<double>& lambda, const std::vector<double>& du, std::vector<double>& dKu) const{
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();

    const size_t kw = mesh->elem_info->get_k_dimension();
    std::vector<gp_Pnt> points(node_num);
    std::vector<long> u_pos(2*kw);
    std::vector<long> l_pos(bnum);
    math::Vector u1(kw);
    math::Vector u2(kw);
    math::Vector duv(2*kw);
    math::Vector l_e(bnum);
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < bnum; ++i){
            points[i] = e.b1->nodes[i]->point;
        }
        for(size_t i = 0; i < bnum; ++i){
            l_pos[i] = mesh->lag_node_map.at(e.b1->nodes[i]->id);
            l_e[i] = lambda[l_pos[i]];
        }
        for(size_t i = 0; i < node_num; ++i){
            const auto n1 = e.b1->parent->nodes[i];
            const auto n2 = e.b2->parent->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                const auto p1 = node_positions[n1->u_pos[j]];
                const auto p2 = node_positions[n2->u_pos[j]];
                u_pos[dof*i + j] = p1;
                u_pos[dof*(node_num + i) + j] = p2;
                u1[dof*i + j] = u_ext[n1->u_pos[j]];
                u2[dof*i + j] = u_ext[n2->u_pos[j]];
                if(p1 >= 0){
                    duv[dof*i + j] = du[p1];
                }
                if(p2 >= 0){
                    duv[dof*(node_num + i) + j] = du[p2];
                }
            }
        }
        const auto MM(EPS_PENALTY*e.b1->parent->get_MnMn(e.b2->parent, u_ext, points, -e.b1->normal));
        const auto uL(e.elem->fl2_uL(u1, u2));
        const auto dKue = MM*duv;
        for(size_t i = 0; i < 2*kw; ++i){
            if(u_pos[i] >= 0){
                dKu[u_pos[i]] += dKue[i];
            }
        }
    }
}
void GlobalStiffnessMatrix::add_frictionless_log(const Meshing * const mesh, const std::vector<long>& node_positions, const std::vector<double>& u_ext, const std::vector<double>& lambda, bool stub){
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();
    const size_t u_size = mesh->load_vector[0].size();

    std::vector<double> u_ext_zero;
    if(stub){
        u_ext_zero.resize(u_ext.size(), 0);
    }

    const size_t kw = mesh->elem_info->get_k_dimension();
    std::vector<gp_Pnt> points(node_num);
    math::Vector L(node_num);
    std::vector<long> u_pos(2*kw);
    std::vector<long> Mn_pos(1, u_size);
    size_t cur_id = 0;
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < bnum; ++i){
            points[i] = e.b1->nodes[i]->point;
        }
        for(size_t i = 0; i < node_num; ++i){
            const auto id = this->lag_pos_cache[cur_id*node_num + i];
            if(id > -1){
                L[i] = lambda[id];
            } else {
                L[i] = 0;
            }
        }
        ++cur_id;
        for(size_t i = 0; i < node_num; ++i){
            const auto n1 = e.b1->parent->nodes[i];
            const auto n2 = e.b2->parent->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                u_pos[dof*i + j] = node_positions[n1->u_pos[j]];
                u_pos[dof*(node_num + i) + j] = node_positions[n2->u_pos[j]];
            }
        }
        if(!stub){
            const auto MM(e.b1->parent->get_MnMn_log(e.b2->parent, u_ext, points, -e.b1->normal, HMU, HK, mesh->thickness, L));
            this->insert_element_matrix(MM, u_pos);
        } else {
            const auto MM(e.b1->parent->get_MnMn_log(e.b2->parent, u_ext_zero, points, -e.b1->normal, HMU, HK, mesh->thickness, L));
            this->reserve_element_matrix(MM, u_pos);
        }
    }
    if(!stub){
        this->final_flush_matrix();
    }
}

void GlobalStiffnessMatrix::add_frictionless_simple(const Meshing * const mesh, const std::vector<long>& node_positions, const std::vector<double>& u_ext, const std::vector<double>& lambda, const bool stub){
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();
    const size_t u_size = mesh->load_vector[0].size();

    const size_t kw = mesh->elem_info->get_k_dimension();

    std::vector<gp_Pnt> points(bnum);
    std::vector<long> u_pos(2*kw);
    std::vector<long> l_pos(bnum);
    std::vector<double> u1(kw), u2(kw);
    std::vector<double> l_e(bnum);
    std::vector<double> l_e_stub(bnum, 10);

    math::Matrix uu(2*kw, 2*kw, 0); 
    math::Matrix uL(2*kw, bnum, 0); 
    math::Matrix LL(bnum, bnum, 0); 
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < bnum; ++i){
            points[i] = e.b1->nodes[i]->point;
        }
        for(size_t i = 0; i < node_num; ++i){
            const auto n1 = e.b1->parent->nodes[i];
            const auto n2 = e.b2->parent->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                u_pos[dof*i + j] = node_positions[n1->u_pos[j]];
                u_pos[dof*(node_num + i) + j] = node_positions[n2->u_pos[j]];
            }
        }
        for(size_t i = 0; i < bnum; ++i){
            l_pos[i] = mesh->lag_node_map.at(e.b1->nodes[i]->id);
            l_e[i] = lambda[l_pos[i]];
            l_pos[i] += u_size;
        }
        for(size_t i = 0; i < node_num; ++i){
            const auto n1 = e.b1->parent->nodes[i];
            const auto n2 = e.b2->parent->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                if(!stub){
                    u1[dof*i + j] = u_ext[n1->u_pos[j]];
                    u2[dof*i + j] = u_ext[n2->u_pos[j]];
                } else {
                    u2[dof*i + j] = -e.b1->normal.Coord(1+j);
                }
            }
        }

        if(!stub){
            uu = this->LAG_DISPL_SIMPLE*e.elem->fl2_uu(l_e, u1, u2);
            uL = this->LAG_DISPL_SIMPLE*e.elem->fl2_uL(u1, u2);
            this->insert_element_matrix(uu, u_pos);
            this->insert_block_symmetric(uL, u_pos, l_pos);
        } else {
            uu = this->LAG_DISPL_SIMPLE*e.elem->fl2_uu(l_e_stub, u1, u2);
            uL = this->LAG_DISPL_SIMPLE*e.elem->fl2_uL(u1, u2);
            this->reserve_element_matrix(uu, u_pos);
            this->reserve_block_symmetric(uL, u_pos, l_pos);
        }
    }
    if(!stub){
        this->final_flush_matrix();
    }
}

void GlobalStiffnessMatrix::append_Ku_frictionless_simple(const Meshing* const mesh, const std::vector<double>& u, std::vector<double>& Ku) const{
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();
    const size_t u_size = mesh->load_vector[0].size();

    const size_t kw = mesh->elem_info->get_k_dimension();
    std::vector<long> u1_pos(kw);
    std::vector<long> u2_pos(kw);
    std::vector<long> l_pos(bnum);
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < bnum; ++i){
            l_pos[i] = mesh->lag_node_map.at(e.b1->nodes[i]->id) + u_size;
        }
        for(size_t i = 0; i < node_num; ++i){
            const auto n1 = e.b1->parent->nodes[i];
            const auto n2 = e.b2->parent->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                u1_pos[dof*i + j] = mesh->node_positions[0][n1->u_pos[j]];
                u2_pos[dof*i + j] = mesh->node_positions[0][n2->u_pos[j]];
            }
        }
        e.elem->fl2_Ku_lambda(this->LAG_DISPL_SIMPLE, u1_pos, u2_pos, l_pos, u, Ku);
    }
}

void GlobalStiffnessMatrix::append_dKu_frictionless_simple(const Meshing* const mesh, const std::vector<double>& u, const std::vector<double>& du, std::vector<double>& Ku) const{
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();
    const size_t u_size = mesh->load_vector[0].size();

    const size_t kw = mesh->elem_info->get_k_dimension();
    std::vector<long> u1_pos(kw);
    std::vector<long> u2_pos(kw);
    std::vector<long> l_pos(bnum);
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < bnum; ++i){
            l_pos[i] = mesh->lag_node_map.at(e.b1->nodes[i]->id) + u_size;
        }
        for(size_t i = 0; i < node_num; ++i){
            const auto n1 = e.b1->parent->nodes[i];
            const auto n2 = e.b2->parent->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                u1_pos[dof*i + j] = mesh->node_positions[0][n1->u_pos[j]];
                u2_pos[dof*i + j] = mesh->node_positions[0][n2->u_pos[j]];
            }
        }
        e.elem->fl2_dKu_lambda(this->LAG_DISPL_SIMPLE, u1_pos, u2_pos, l_pos, u, du, Ku);
    }
}

double GlobalStiffnessMatrix::constraint_frictionless_log(const Meshing* const mesh, const std::vector<double>& u_ext, const std::vector<double>& lambda) const{
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();

    std::vector<gp_Pnt> points(node_num);
    math::Vector L(node_num);
    size_t cur_id = 0;
    double constr = -this->LOG_TOL;
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < bnum; ++i){
            points[i] = e.b1->nodes[i]->point;
        }
        for(size_t i = 0; i < node_num; ++i){
            const auto id = this->lag_pos_cache[cur_id*node_num + i];
            if(id > -1){
                L[i] = lambda[id];
            } else {
                L[i] = 0;
            }
        }
        ++cur_id;
        constr += e.b1->parent->get_log_integ(e.b2->parent, u_ext, points, -e.b1->normal, HMU, HK, mesh->thickness, L);
    }

    return constr;
}

void GlobalStiffnessMatrix::append_Ku_frictionless_log(const Meshing* const mesh, const std::vector<double>& u_ext, const std::vector<double>& lambda, std::vector<double>& Ku) const{
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();

    std::vector<gp_Pnt> points(node_num);
    math::Vector L(node_num);
    size_t cur_id = 0;
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < bnum; ++i){
            points[i] = e.b1->nodes[i]->point;
        }
        for(size_t i = 0; i < node_num; ++i){
            const auto id = this->lag_pos_cache[cur_id*node_num + i];
            if(id > -1){
                L[i] = lambda[id];
            } else {
                L[i] = 0;
            }
        }
        ++cur_id;
        e.b1->parent->Ku_log(e.b2->parent, mesh->node_positions[0], u_ext, points, -e.b1->normal, Ku, HMU, HK, mesh->thickness, L);
    }
}

void GlobalStiffnessMatrix::append_dKu_frictionless_log(const Meshing* const mesh, const std::vector<double>& u, const std::vector<double>& du, const std::vector<double>& lambda, std::vector<double>& Ku) const{
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();

    std::vector<gp_Pnt> points(node_num);
    math::Vector L(node_num);
    size_t cur_id = 0;
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < bnum; ++i){
            points[i] = e.b1->nodes[i]->point;
        }
        for(size_t i = 0; i < node_num; ++i){
            const auto id = this->lag_pos_cache[cur_id*node_num + i];
            if(id > -1){
                L[i] = lambda[id];
            } else {
                L[i] = 0;
            }
        }
        ++cur_id;
        e.b1->parent->dKu_log(e.b2->parent, mesh->node_positions[0], u, du, points, -e.b1->normal, Ku, HMU, HK, mesh->thickness, L);
    }
}
