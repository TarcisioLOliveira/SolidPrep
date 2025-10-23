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
#include "cblas.h"
#include "logger.hpp"
#include "math/matrix.hpp"
#include <limits>

GlobalStiffnessMatrix::GlobalStiffnessMatrix(double EPS_DISPL_SIMPLE):
    EPS_PENALTY(EPS_DISPL_SIMPLE), EPS_DISPL(EPS_DISPL_SIMPLE), LAG_DISPL_SIMPLE(EPS_DISPL_SIMPLE){

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
    if(this->first_time){
        this->first_time = false;
        if(type == FiniteElement::ContactType::FRICTIONLESS_DISPL_LOG){
            this->add_frictionless_log(mesh, node_positions, u_ext, true);
        } else if(type == FiniteElement::ContactType::FRICTIONLESS_DISPL_SIMPLE){
            const std::vector<double> empty_lambda(lambda.size(), 0);
            this->add_frictionless_simple(mesh, node_positions, u_ext, empty_lambda, true);
        } else if(type == FiniteElement::ContactType::FRICTIONLESS_DISPL_CONSTR){
            this->add_frictionless_part2(mesh, node_positions, u_ext, lambda, D_cache, topopt, true);
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
        const math::Matrix k(e->get_k(D, t));
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
    const double t = mesh->thickness;
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnode_num = mesh->elem_info->get_boundary_nodes_per_element();

    std::vector<gp_Pnt> points(bnode_num);

    std::vector<long> u_pos(dof*node_num);
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
            const math::Matrix R(e->get_R(K, t, points));
            this->insert_element_matrix(R, u_pos);
        }
    }
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
        //logger::quick_log(e.b1->geom_id, e.b2->geom_id);
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
        //logger::quick_log(e.b1->geom_id, e.b2->geom_id);
        const auto MM(EPS_PENALTY*e.b1->parent->get_MnMn(e.b2->parent, u_ext, points, -e.b1->normal));
        const auto uL(e.elem->fl2_uL(l_e, u1, u2));
        //const auto Kue = uL*l_e + MM*uv;
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
        //logger::quick_log(e.b1->geom_id, e.b2->geom_id);
        const auto MM(EPS_PENALTY*e.b1->parent->get_MnMn(e.b2->parent, u_ext, points, -e.b1->normal));
        const auto uL(e.elem->fl2_uL(l_e, u1, u2));
        //const auto dKue = uL*l_e + MM*duv;
        const auto dKue = MM*duv;
        for(size_t i = 0; i < 2*kw; ++i){
            if(u_pos[i] >= 0){
                dKu[u_pos[i]] += dKue[i];
            }
        }
    }
}
void GlobalStiffnessMatrix::add_frictionless_part2(const Meshing * const mesh, const std::vector<long>& node_positions, const std::vector<double>& u_ext, const std::vector<double>& lambda, const std::vector<math::Matrix>& D_cache, bool topopt, bool stub){
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();
    const size_t u_size = mesh->load_vector[0].size();
    const size_t max_size = u_size + lambda.size();

    const size_t kw = mesh->elem_info->get_k_dimension();
    const size_t l_num = mesh->lag_node_map.size();
    std::vector<long> u_pos(kw);
    std::vector<long> l_pos(3*bnum);
    std::vector<double> u_e(kw);
    std::vector<double> l_n(bnum);
    std::vector<double> l_p1(bnum);
    std::vector<double> l_p2(bnum);

    const double MULT = stub ? 0 : 1;

    const auto insert_elements = [&](const PairedBoundaryElements& e, const math::Matrix& D)->void{
        for(size_t i = 0; i < bnum; ++i){
            const size_t j = i + bnum;
            const size_t k = j + bnum;
            l_pos[i] = mesh->lag_node_map.at(e.b1->nodes[i]->id);
            l_pos[j] = l_pos[i] + l_num;
            l_pos[k] = l_pos[j] + l_num;
            l_n[i] = lambda[l_pos[i]];
            l_p1[i] = lambda[l_pos[j]];
            l_p2[i] = lambda[l_pos[k]];
            l_pos[i] += u_size;
            l_pos[j] += u_size;
            l_pos[k] += u_size;
        }
        for(size_t i = 0; i < node_num; ++i){
            const auto n1 = e.b1->parent->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                u_pos[dof*i + j] = node_positions[n1->u_pos[j]];
                //u_e[dof*i + j] = u_ext[n1->u_pos[j]];
            }
        }
        const auto uL(MULT*e.elem->fl3_uL(D, l_n));
        const auto LL(MULT*e.elem->fl3_LL(D, l_n, l_p1, l_p2, u_ext));

        this->insert_block_symmetric(uL, u_pos, l_pos);
        this->insert_element_matrix(LL, l_pos);
    };

    size_t current_geom = 0;
    size_t D_offset = 0;
    math::Matrix D_const;
    for(const auto& e:mesh->paired_boundary){
        const auto data = mesh->contact_data.at(e.elem->e1);
        const auto g = mesh->geometries[data.geom_id];
        while(data.geom_id != current_geom){
            if((topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous()){
                D_offset += g->mesh.size();
            } else {
                D_const = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
            }
            ++current_geom;
        }
        if((topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous()){
            const auto& D = D_cache[D_offset + data.mesh_pos];
            insert_elements(e, D);
        } else {
            insert_elements(e, D_const);
        }
    }
    if(!stub){
        for(size_t i = u_size; i < max_size; ++i){
            this->add_to_matrix(i, i, this->K_MIN);
        }
        this->final_flush_matrix();
    }
}
void GlobalStiffnessMatrix::add_frictionless_log(const Meshing * const mesh, const std::vector<long>& node_positions, const std::vector<double>& u_ext, bool stub){
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();
    const size_t u_size = mesh->load_vector[0].size();

    const size_t kw = mesh->elem_info->get_k_dimension();
    std::vector<gp_Pnt> points(node_num);
    std::vector<long> u_pos(2*kw);
    std::vector<long> Mn_pos(1, u_size);
    if(!stub){
        double constr = -this->LOG_TOL;
        for(const auto& e:mesh->paired_boundary){
            for(size_t i = 0; i < bnum; ++i){
                points[i] = e.b1->nodes[i]->point;
            }
            constr += e.b1->parent->get_log_integ(e.b2->parent, u_ext, points, -e.b1->normal, HC, HK);
        }
        constr = std::max(constr, 0.0);
        const double mult = this->LAG_DISPL_LOG + this->MU_LOG*constr;
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
            //logger::quick_log(e.b1->geom_id, e.b2->geom_id);
            const auto Kue = e.b1->parent->Kue_log(e.b2->parent, u_ext, points, -e.b1->normal, HC, HK);
            const auto MM(mult*e.b1->parent->get_MnMn_log(e.b2->parent, u_ext, points, -e.b1->normal, HC, HK)
                    + MU_LOG*Kue*Kue.T());
            this->insert_element_matrix(MM, u_pos);
            this->insert_block_symmetric(Kue, u_pos, Mn_pos);
        }
        this->final_flush_matrix();
    } else {
        const math::Matrix Mn(2*kw, 1, 0);
        for(const auto& e:mesh->paired_boundary){
            for(size_t i = 0; i < node_num; ++i){
                const auto n1 = e.b1->parent->nodes[i];
                const auto n2 = e.b2->parent->nodes[i];
                for(size_t j = 0; j < dof; ++j){
                    u_pos[dof*i + j] = node_positions[n1->u_pos[j]];
                    u_pos[dof*(node_num + i) + j] = node_positions[n2->u_pos[j]];
                }
            }
            this->insert_block_symmetric(Mn, u_pos, Mn_pos);
        }
    }
}

void GlobalStiffnessMatrix::append_Ku_frictionless(const Meshing* const mesh, const std::vector<double>& u, std::vector<double>& Ku, const std::vector<math::Matrix>& D_cache, bool topopt) const{
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t l_num = mesh->lag_node_map.size();
    const size_t dof = mesh->elem_info->get_dof_per_node();
    const size_t u_size = mesh->load_vector[0].size();

    const size_t kw = mesh->elem_info->get_k_dimension();
    std::vector<long> u_pos(kw);
    std::vector<long> l_pos(3*bnum);

    const auto insert_elements = [&](const PairedBoundaryElements& e, const math::Matrix& D)->void{
        for(size_t i = 0; i < bnum; ++i){
            const size_t j = i + bnum;
            const size_t k = j + bnum;
            l_pos[i] = mesh->lag_node_map.at(e.b1->nodes[i]->id) + u_size;
            l_pos[j] = l_pos[i] + l_num;
            l_pos[k] = l_pos[j] + l_num;
        }
        for(size_t i = 0; i < node_num; ++i){
            const auto n1 = e.b1->parent->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                u_pos[dof*i + j] = mesh->node_positions[0][n1->u_pos[j]];
            }
        }
        e.elem->fl3_Ku(D, u_pos, l_pos, u, Ku);
    };

    size_t current_geom = 0;
    size_t D_offset = 0;
    math::Matrix D_const;
    for(const auto& e:mesh->paired_boundary){
        const auto data = mesh->contact_data.at(e.elem->e1);
        const auto g = mesh->geometries[data.geom_id];
        while(data.geom_id != current_geom){
            if((topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous()){
                D_offset += g->mesh.size();
            } else {
                D_const = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
            }
            ++current_geom;
        }
        if((topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous()){
            const auto& D = D_cache[D_offset + data.mesh_pos];
            insert_elements(e, D);
        } else {
            insert_elements(e, D_const);
        }
    }
}

void GlobalStiffnessMatrix::append_dKu_frictionless(const Meshing* const mesh, const std::vector<double>& u, const std::vector<double>& du, const double eta, std::vector<double>& Ku, const std::vector<math::Matrix>& D_cache, bool topopt) const{
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t l_num = mesh->lag_node_map.size();
    const size_t dof = mesh->elem_info->get_dof_per_node();
    const size_t u_size = mesh->load_vector[0].size();

    const size_t kw = mesh->elem_info->get_k_dimension();
    std::vector<long> u_pos(kw);
    std::vector<long> l_pos(3*bnum);

    const auto insert_elements = [&](const PairedBoundaryElements& e, const math::Matrix& D)->void{
        for(size_t i = 0; i < bnum; ++i){
            const size_t j = i + bnum;
            const size_t k = j + bnum;
            l_pos[i] = mesh->lag_node_map.at(e.b1->nodes[i]->id) + u_size;
            l_pos[j] = l_pos[i] + l_num;
            l_pos[k] = l_pos[j] + l_num;
        }
        for(size_t i = 0; i < node_num; ++i){
            const auto n1 = e.b1->parent->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                u_pos[dof*i + j] = mesh->node_positions[0][n1->u_pos[j]];
            }
        }
        e.elem->fl3_dKu(D, eta, u_pos, l_pos, u, du, Ku);
    };

    size_t current_geom = 0;
    size_t D_offset = 0;
    math::Matrix D_const;
    for(const auto& e:mesh->paired_boundary){
        const auto data = mesh->contact_data.at(e.elem->e1);
        const auto g = mesh->geometries[data.geom_id];
        while(data.geom_id != current_geom){
            if((topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous()){
                D_offset += g->mesh.size();
            } else {
                D_const = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
            }
            ++current_geom;
        }
        if((topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous()){
            const auto& D = D_cache[D_offset + data.mesh_pos];
            insert_elements(e, D);
        } else {
            insert_elements(e, D_const);
        }
    }
}

void GlobalStiffnessMatrix::add_frictionless_simple(const Meshing * const mesh, const std::vector<long>& node_positions, const std::vector<double>& u_ext, const std::vector<double>& lambda, bool stub){
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();
    const size_t u_size = mesh->load_vector[0].size();
    const size_t max_size = u_size + lambda.size();

    const size_t kw = mesh->elem_info->get_k_dimension();

    std::vector<gp_Pnt> points(bnum);
    std::vector<long> u_pos(2*kw);
    std::vector<long> l_pos(bnum);
    std::vector<double> u1(kw), u2(kw);
    std::vector<double> l_e(bnum);

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
                u1[dof*i + j] = u_ext[n1->u_pos[j]];
                u2[dof*i + j] = u_ext[n2->u_pos[j]];
            }
        }
        if(!stub){ 
            uu = this->LAG_DISPL_SIMPLE*e.elem->fl2_uu(l_e, u1, u2);
            uL = this->LAG_DISPL_SIMPLE*e.elem->fl2_uL(l_e, u1, u2);
            //LL = this->LAG_DISPL_SIMPLE*e.elem->fl2_LL(l_e, u1, u2);
        }

        this->insert_element_matrix(uu, u_pos);
        this->insert_block_symmetric(uL, u_pos, l_pos);
        //this->insert_element_matrix(LL, l_pos);
    }
    if(!stub){
        //for(size_t i = u_size; i < max_size; ++i){
        //    //this->add_to_matrix(i, i, 0);
        //    this->add_to_matrix(i, i, 0);
        //    //this->add_to_matrix(i, i, 5);
        //}
        this->final_flush_matrix();
    //} else {
    //    for(size_t i = u_size; i < max_size; ++i){
    //        this->add_to_matrix(i, i, 0);
    //    }
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

void GlobalStiffnessMatrix::append_Ku_frictionless_log(const Meshing* const mesh, const std::vector<double>& u_ext, std::vector<double>& Ku) const{
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();

    std::vector<gp_Pnt> points(node_num);
    double constr = -this->LOG_TOL;
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < bnum; ++i){
            points[i] = e.b1->nodes[i]->point;
        }
        constr += e.b1->parent->get_log_integ(e.b2->parent, u_ext, points, -e.b1->normal, HC, HK);
    }
    Ku.back() = constr;
    const double mult = this->LAG_DISPL_LOG + this->MU_LOG*constr;
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < bnum; ++i){
            points[i] = e.b1->nodes[i]->point;
        }
        e.b1->parent->Ku_log(mult, e.b2->parent, mesh->node_positions[0], u_ext, points, -e.b1->normal, Ku, HC, HK);
    }
}

void GlobalStiffnessMatrix::append_dKu_frictionless_log(const Meshing* const mesh, const std::vector<double>& u_ext, const std::vector<double>& u, const std::vector<double>& du, const double eta, std::vector<double>& Ku) const{
    (void)eta;
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();

    std::vector<gp_Pnt> points(node_num);
    double constr = -LOG_TOL;
    double constr_deriv = 0;
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < bnum; ++i){
            points[i] = e.b1->nodes[i]->point;
        }
        constr += e.b1->parent->get_log_integ(e.b2->parent, u_ext, points, -e.b1->normal, HC, HK);
        constr_deriv += e.b1->parent->get_log_integ_deriv(e.b2->parent, mesh->node_positions[0], u, du, points, -e.b1->normal, HC, HK);
    }
    constr = std::max(constr, 0.0);
    Ku.back() = constr_deriv;
    const double mult = this->LAG_DISPL_LOG + this->MU_LOG*constr;
    const double mult_deriv = du.back() + this->MU_LOG*constr_deriv;
    //const double mult_deriv = this->MU_LOG*constr_deriv;
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < bnum; ++i){
            points[i] = e.b1->nodes[i]->point;
        }
        e.b1->parent->Ku_log(mult_deriv, e.b2->parent, mesh->node_positions[0], u_ext, points, -e.b1->normal, Ku, HC, HK);
        e.b1->parent->dKu_log(mult, e.b2->parent, mesh->node_positions[0], u, du, points, -e.b1->normal, Ku, HC, HK);
    }
}
