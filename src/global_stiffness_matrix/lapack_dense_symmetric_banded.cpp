/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
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

#include "global_stiffness_matrix/lapack_dense_symmetric_banded.hpp"
#include "logger.hpp"

namespace global_stiffness_matrix{

void LAPACKDenseSymmetricBanded::calculate_dimensions(const Meshing* const mesh, const std::vector<double>& load){
    const size_t k_dim    = mesh->elem_info->get_k_dimension();
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();

    this->recalculated_dimensions = true;

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
}

void LAPACKDenseSymmetricBanded::generate(const Meshing* const mesh, const std::vector<double>& density, const double pc, const double psi){
        if(this->recalculated_dimensions){
        this->K.clear();
        this->K.resize(W*N,0);
        this->recalculated_dimensions = false;
    } else {
        std::fill(this->K.begin(), this->K.end(), 0);
    }
    logger::quick_log("Generating stiffness matrix...");
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
    logger::quick_log("Done.");
}

void LAPACKDenseSymmetricBanded::add_geometry(const Meshing* const mesh, const Geometry* const g){
    const auto D = g->materials.get_D();
    const double t = mesh->thickness;
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();

    for(auto& e : g->mesh){
        std::vector<long> u_pos;
        u_pos.reserve(dof*node_num);
        for(size_t i = 0; i < node_num; ++i){
            const auto& n = e->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                u_pos.push_back(n->u_pos[j]);
            }
        }
        std::vector<double> k = e->get_k(D, t);
        this->insert_element_matrix(k, u_pos, N);
    }
}

void LAPACKDenseSymmetricBanded::add_geometry(const Meshing* const mesh, const Geometry* const g, std::vector<double>::const_iterator& rho, const double pc, const double psi){
    const double t = mesh->thickness;
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();

    //const size_t num_den = g->number_of_densities_needed();
    const size_t num_mat = g->number_of_materials();
    if(num_mat == 1){
        const auto D = g->materials.get_D();
        auto rhoD = D;
        for(const auto& e:g->mesh){
            std::vector<long> u_pos;
            u_pos.reserve(dof*node_num);
            for(size_t i = 0; i < node_num; ++i){
                const auto& n = e->nodes[i];
                for(size_t j = 0; j < dof; ++j){
                    u_pos.push_back(n->u_pos[j]);
                }
            }
            const double rhop = K_MIN + (1-K_MIN)*std::pow(*rho, pc);
            for(size_t i = 0; i < D.size(); ++i){
                rhoD[i] = rhop*D[i];
            }
            const std::vector<double> k = e->get_k(rhoD, t);
            this->insert_element_matrix(k, u_pos, N);
            ++rho;
        }
    } else {
        auto D = g->materials.get_D();
        for(const auto& e:g->mesh){
            std::vector<long> u_pos;
            u_pos.reserve(dof*node_num);
            for(size_t i = 0; i < node_num; ++i){
                const auto& n = e->nodes[i];
                for(size_t j = 0; j < dof; ++j){
                    u_pos.push_back(n->u_pos[j]);
                }
            }
            g->materials.get_D(rho, g->with_void, pc, this->K_MIN, psi, D);
            const std::vector<double> k = e->get_k(D, t);
            this->insert_element_matrix(k, u_pos, N);
        }
    }
}

}
