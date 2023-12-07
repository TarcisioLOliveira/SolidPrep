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

void GlobalStiffnessMatrix::generate_base(const Meshing * const mesh, const std::vector<double>& density, const double pc, const double psi){
    if(density.size() == 0){
        for(auto& g : mesh->geometries){
            this->add_geometry(mesh, g);
        }
    } else {
        auto rho = density.begin();
        for(auto& g : mesh->geometries){
            if(g->do_topopt){
                this->add_geometry(mesh, g, rho, pc, psi);
            } else {
                this->add_geometry(mesh, g);
            }
        }
    }
    if(mesh->springs->size() > 0){
        this->add_springs(mesh);
    }
    for(size_t i = 0; i < mesh->load_vector.size(); ++i){
        this->add_to_matrix(i, i, this->K_MIN);
    }
}

void GlobalStiffnessMatrix::add_geometry(const Meshing* const mesh, const Geometry* const g){
    const double t = mesh->thickness;
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();

    std::vector<long> u_pos(dof*node_num);
    for(auto& e : g->mesh){
        const gp_Pnt c = e->get_centroid();
        const auto D = g->materials.get_D(e.get(), c);
        for(size_t i = 0; i < node_num; ++i){
            const auto& n = e->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                u_pos[i*dof + j] = n->u_pos[j];
            }
        }
        const std::vector<double> k = e->get_k(D, t);
        this->insert_element_matrix(k, u_pos);
    }
}

void GlobalStiffnessMatrix::add_springs(const Meshing * const mesh){
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
                    u_pos[i*dof + j] = n->u_pos[j];
                }
            }
            const std::vector<double> R = e->get_R(K, t, points);
            this->insert_element_matrix(R, u_pos);
        }
    }
}

void GlobalStiffnessMatrix::add_geometry(const Meshing* const mesh, const Geometry* const g, std::vector<double>::const_iterator& rho, const double pc, const double psi){
    const double t = mesh->thickness;
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();

    std::vector<long> u_pos(dof*node_num);
    const size_t num_den = g->number_of_densities_needed();
    const size_t k_size = mesh->elem_info->get_k_dimension();
    const size_t s_size = mesh->elem_info->get_D_dimension();
    auto D = std::vector<double>(s_size*s_size, 0);
    if(g->with_void){
        for(const auto& e:g->mesh){
            if(*rho > 0){
                for(size_t i = 0; i < node_num; ++i){
                    const auto& n = e->nodes[i];
                    for(size_t j = 0; j < dof; ++j){
                        u_pos[i*dof + j] = n->u_pos[j];
                    }
                }
                const gp_Pnt c = e->get_centroid();
                g->materials.get_D(rho, psi, e.get(), c, D);
                std::vector<double> k = e->get_k(D, t);
                cblas_dscal(k_size*k_size, std::pow(*rho, pc), k.data(), 1);
                this->insert_element_matrix(k, u_pos);
            }

            rho += num_den;
        }
    } else {
        for(const auto& e:g->mesh){
            for(size_t i = 0; i < node_num; ++i){
                const auto& n = e->nodes[i];
                for(size_t j = 0; j < dof; ++j){
                    u_pos[i*dof + j] = n->u_pos[j];
                }
            }
            const gp_Pnt c = e->get_centroid();
            g->materials.get_D(rho, psi, e.get(), c, D);
            const std::vector<double> k = e->get_k(D, t);
            this->insert_element_matrix(k, u_pos);

            rho += num_den;
        }
    }
}

void GlobalStiffnessMatrix::calculate_dimensions(const Meshing* const mesh, const std::vector<double>& load){
    const size_t k_dim    = mesh->elem_info->get_k_dimension();
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();

    W = load.size();
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
                    pos.push_back(n->u_pos[j]);
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
