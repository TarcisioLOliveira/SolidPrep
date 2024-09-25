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
#include "project_data.hpp"
#include <limits>

void GlobalStiffnessMatrix::generate_base(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache, const std::vector<double>& u_ext, const FiniteElement::MatrixType type){
    size_t D_offset = 0;
    size_t L = u_size;
    if(type == FiniteElement::MatrixType::FRICTIONLESS){
        L += l_num;
    }
    const std::vector<double> lambda(l_num, 0);
    //if(type == FiniteElement::MatrixType::LAMBDA_SLIDING){
    //    L += 2*l_num;
    //} else if(type == FiniteElement::MatrixType::LAMBDA_HESSIAN){
    //    L += 3*l_num;
    //}
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
    if(mesh->paired_boundary.size() > 0){
        this->add_contacts(mesh, node_positions, u_ext);
        //this->add_frictionless_part1(mesh, node_positions, u_ext);
        //this->add_frictionless_part2(mesh, node_positions, u_ext, lambda);
        //for(size_t i = u_size; i < L; ++i){
        //    this->add_to_matrix(i, i, this->K_MIN);
        //}
    }
    //if(type == FiniteElement::MatrixType::LAMBDA_SLIDING){
    //    this->generate_expansion(mesh, u_size, l_num, node_positions, topopt, D_cache, Meshing::LambdaType::PARALLEL);
    //} else if(type == FiniteElement::MatrixType::LAMBDA_HESSIAN){
    //    this->generate_expansion(mesh, u_size, l_num, node_positions, topopt, D_cache, Meshing::LambdaType::ALL);
    //}
    if(topopt){
        for(size_t i = 0; i < u_size; ++i){
            this->add_to_matrix(i, i, this->K_MIN);
        }
    }
}

void GlobalStiffnessMatrix::generate_expansion(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache, const Meshing::LambdaType type){
    const double t = mesh->thickness;
    const size_t bnode_num = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t kw = mesh->elem_info->get_k_dimension();
    const size_t lw = mesh->elem_info->get_dof_per_node();
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    // TODO: Springs
    // TODO: DENSITY

    std::vector<long> u_pos(kw);
    size_t geom = 0;
    const Geometry* g = mesh->geometries[geom];
    bool use_D_cache = (topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous();
    size_t D_offset = 0;
    std::vector<double> D;
    if(!use_D_cache){
        D = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
    }
    std::vector<double> k;
    std::vector<size_t> lv_it(node_num);
    std::vector<const std::vector<size_t>*> lv_vec(node_num);
    for(const auto& l_e:mesh->lambda_affected_elements){
        if(l_e.geom_id != geom){
            if(!use_D_cache){
                D_offset += g->mesh.size();
            }
            geom = l_e.geom_id;
            const Geometry* g = mesh->geometries[geom];
            use_D_cache = (topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous();
            if(!use_D_cache){
                D = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
            }
        }
        if(use_D_cache){
            k = l_e.e->get_k(D_cache[l_e.e->id - D_offset], t);
        } else {
            k = l_e.e->get_k(D, t);
        }
        const size_t ln = l_e.lambdas.size();
        std::vector<double> R(ln*(kw*lw), 0);

        for(size_t i = 0; i < node_num; ++i){
            const auto& n = l_e.e->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                const size_t p = n->id*dof + j;
                u_pos[i*dof + j] = node_positions[p];
            }
        }

        std::fill(lv_it.begin(), lv_it.end(), 0);
        for(size_t i = 0; i < node_num; ++i){
            const auto& vi = mesh->node_lambda_map.find(l_e.e->nodes[i]->id);
            if(vi != mesh->node_lambda_map.end()){
                lv_vec[i] = &vi->second;
            } else {
                lv_vec[i] = nullptr;
            }
        }
        for(size_t l_it = 0; l_it < ln; ++l_it){
            const auto l_i = l_e.lambdas[l_it];
            const auto& lambda = mesh->lambda_elements[l_i];
            std::vector<double> R0 = 
                {lambda.n.X(), lambda.p1.X(), lambda.p2.X(),
                 lambda.n.Y(), lambda.p1.Y(), lambda.p2.Y(),
                 lambda.n.Z(), lambda.p1.Z(), lambda.p2.Z()};

            for(size_t i = 0; i < node_num; ++i){
                if(lv_vec[i] != nullptr && lv_it[i] < lv_vec[i]->size() && lv_vec[i]->at(lv_it[i]) == l_i){
                    for(size_t ii = 0; ii < lw; ++ii){
                        for(size_t jj = 0; jj < lw; ++jj){
                            R[(i*lw + ii)*(lw*ln) + (jj + l_it*lw)] = -R0[ii*lw + jj];
                        }
                    }
                    ++lv_it[i];
                }
            }
        }
        this->insert_expansion_matrices(k, R, u_pos, l_e.lambdas, lw, l_num, u_size, type);
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
        const std::vector<double> k = e->get_k(D, t);
        this->insert_element_matrix(k, u_pos);
    }
}

void GlobalStiffnessMatrix::add_geometry(const Meshing* const mesh, const std::vector<long>& node_positions, const Geometry* const g, const size_t D_offset, const std::vector<std::vector<double>>& D_cache){
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
        std::vector<double> k = e->get_k(D, t);
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
            const std::vector<double> R = e->get_R(K, t, points);
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


void GlobalStiffnessMatrix::insert_expansion_matrices(const std::vector<double>& k, const std::vector<double>& R, const std::vector<long>& u_pos, const std::vector<size_t>& l_i, const size_t dof, const size_t l_num, const size_t u_num, const Meshing::LambdaType type){
    const size_t kw = u_pos.size();
    const size_t lw = dof;
    const size_t ln = l_i.size();

    std::vector<long> l_pos(lw*ln);
    for(size_t i = 0; i < l_i.size(); ++i){
        l_pos[3*i + 0] = (type != Meshing::LambdaType::PARALLEL) ? static_cast<long>(u_num + l_i[i] + 2*l_num) : -1;
        l_pos[3*i + 1] = (type != Meshing::LambdaType::NORMAL) ? static_cast<long>(u_num + l_i[i]) : -1;
        l_pos[3*i + 2] = (type != Meshing::LambdaType::NORMAL) ? static_cast<long>(u_num + l_i[i] + l_num) : -1;
    }

    std::vector<double> KR(kw*lw*l_i.size(), 0);
    // No idea why the BLAS multiplication is not working
    // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, kw, lw, kw, 1, k.data(), kw, R.data(), kw, 0, KR.data(), kw);
    for(size_t i = 0; i < kw; ++i){
        for(size_t j = 0; j < lw*ln; ++j){
            for(size_t l = 0; l < kw; ++l){
                KR[i*lw*ln + j] += k[i*kw + l]*R[l*lw*ln + j];
            }
        }
    }
    
    this->insert_block_symmetric(KR, u_pos, l_pos);

    std::vector<double> RKR(lw*lw*ln*ln, 0);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, lw*ln, lw*ln, kw, 1, R.data(), lw*ln, KR.data(), lw*ln, 0, RKR.data(), lw*ln);

    this->insert_element_matrix(RKR, l_pos);
}


void GlobalStiffnessMatrix::add_contacts(const Meshing * const mesh, const std::vector<long>& node_positions, const std::vector<double>& u_ext){
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();

    const size_t kw = mesh->elem_info->get_k_dimension();
    std::vector<gp_Pnt> points(node_num);
    std::vector<long> u_pos(2*kw);
    size_t added = 0;
    size_t negative = 0;
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
        auto MM = e.b1->parent->get_MnMn(e.b2->parent, u_ext, points, e.b1->normal);
        cblas_dscal(MM.size(), EPS, MM.data(), 1);
        for(size_t i = 0; i < 2*kw; ++i){
            if(std::abs(MM[i*2*kw + i]) > 1e-7){
                ++added;
                break;
            }
        }
        for(size_t i = 0; i < 2*kw; ++i){
            if(std::abs(MM[i*2*kw + i]) > 1e-7 && MM[i*2*kw + i] < 0){
                ++negative;
                break;
            }
        }
        this->insert_element_matrix(MM, u_pos);
    }
    logger::quick_log("added", added);
    logger::quick_log("negative", negative);
}

void GlobalStiffnessMatrix::add_frictionless_part1(const Meshing * const mesh, const std::vector<long>& node_positions, const std::vector<double>& u_ext){
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();

    const size_t kw = mesh->elem_info->get_k_dimension();
    std::vector<gp_Pnt> points(bnum);
    std::vector<long> u_pos(2*kw);
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
        auto uu = e.b1->parent->get_MnMn(e.b2->parent, u_ext, points, e.b1->normal);
        cblas_dscal(uu.size(), this->EPS, uu.data(), 1);
        this->insert_element_matrix(uu, u_pos);
    }
}

void GlobalStiffnessMatrix::add_frictionless_part2(const Meshing * const mesh, const std::vector<long>& node_positions, const std::vector<double>& u_ext, const std::vector<double>& lambda){
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();
    const size_t u_size = mesh->load_vector[0].size();

    const size_t kw = mesh->elem_info->get_k_dimension();
    std::vector<long> u_pos(2*kw);
    std::vector<long> l_pos(bnum);
    std::vector<double> u1(kw), u2(kw);
    std::vector<double> l_e(bnum);
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < bnum; ++i){
            l_pos[i] = mesh->lag_node_map.at(e.b1->nodes[i]->id);
            l_e[i] = lambda[l_pos[i]];
            l_pos[i] += u_size;
        }
        for(size_t i = 0; i < node_num; ++i){
            const auto n1 = e.b1->parent->nodes[i];
            const auto n2 = e.b2->parent->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                u_pos[dof*i + j] = node_positions[n1->u_pos[j]];
                u_pos[dof*(node_num + i) + j] = node_positions[n2->u_pos[j]];
                u1[dof*i + j] = u_ext[n1->u_pos[j]];
                u2[dof*i + j] = u_ext[n2->u_pos[j]];
            }
        }
        //logger::quick_log(e.b1->geom_id, e.b2->geom_id);
        auto LL = e.elem->fl2_LL(l_e, u1, u2);
        auto uL = e.elem->fl2_uL(l_e);
        cblas_dscal(uL.size(), this->EPS, uL.data(), 1);
        cblas_dscal(LL.size(), this->EPS, LL.data(), 1);
        //for(auto& m:MM){
        //    m *= EPS;
        //}
        this->insert_block_symmetric(uL, u_pos, l_pos);
        this->insert_element_matrix(LL, l_pos);
    }
}

void GlobalStiffnessMatrix::append_Ku_frictionless(const Meshing* const mesh, const std::vector<double>& u, std::vector<double>& Ku) const{
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
        e.elem->fl2_Ku_lambda(this->EPS, u1_pos, u2_pos, l_pos, u, Ku);
    }
}
