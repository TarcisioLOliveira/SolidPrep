/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
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

#include <cblas.h>
#include "finite_element.hpp"
#include "logger.hpp"
#include "project_data.hpp"

std::vector<double> FiniteElement::calculate_forces(const Meshing* const mesh, const std::vector<double>& displacements) const{
    logger::quick_log("Calculating forces...");
    const size_t dof       = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    std::vector<double> results(mesh->node_list.size()*dof, 0);
  
   for(auto& g:mesh->geometries){ 
        const auto D = g->get_D(0); 
        for(auto& e:g->mesh){
            auto f = e->get_internal_loads(D, mesh->thickness, displacements);
            for(size_t n = 0; n < node_num; ++n){
                auto& node = e->nodes[n];
                for(size_t i = 0; i < dof; ++i){
                    results[node->id*dof + i] += f[n*dof+i];
                }
            }
        }
    }
    logger::quick_log("Done.");

    return results;
}

void FiniteElement::calculate_dimensions(const Meshing* const mesh, const std::vector<double>& load){
    const size_t k_dim    = mesh->elem_info->get_k_dimension();
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();

    this->recalculated_dimensions = true;

    W = load.size();
    N = k_dim;

    size_t ei = 0;
   for(auto& g:mesh->geometries){ 
        const auto D = g->get_D(0); 
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
            ++ei;
        }
    }
}

void FiniteElement::generate_K(const Meshing* const mesh, const std::vector<double>& density, const double pc){
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
            this->add_geometry_to_K(mesh, g);
        }
    } else {
        auto rho = density.begin();
        for(auto& g : mesh->geometries){
            if(g->do_topopt){
                this->add_geometry_to_K(mesh, g, rho, pc);
            } else {
                this->add_geometry_to_K(mesh, g);
            }
        }
    }
}

void FiniteElement::add_geometry_to_K(const Meshing* const mesh, const Geometry* const g){
    const auto D = g->get_D(0);
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
        this->insert_element_matrix(K, k, u_pos, N);
    }
}

void FiniteElement::add_geometry_to_K(const Meshing* const mesh, const Geometry* const g, std::vector<double>::const_iterator rho, const double pc){
    const auto D = g->get_D(0);
    const double t = mesh->thickness;
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();

    if(g->alternate_materials.empty()){
        for(auto& e : g->mesh){
            std::vector<long> u_pos;
            u_pos.reserve(dof*node_num);
            for(size_t i = 0; i < node_num; ++i){
                const auto& n = e->nodes[i];
                for(size_t j = 0; j < dof; ++j){
                    u_pos.push_back(n->u_pos[j]);
                }
            }
            auto rhoD = D;
            double rho_scal = this->K_MIN + (1-this->K_MIN)*std::pow(*rho, pc);
            for(auto& d:rhoD){
                d *= rho_scal;
            }
            std::vector<double> k = e->get_k(rhoD, t);
            this->insert_element_matrix(K, k, u_pos, N);
            ++rho;
        }
    } else if(g->alternate_materials.size() == 1){
        const auto D2 = g->get_D(1);
        for(auto& e : g->mesh){
            std::vector<long> u_pos;
            u_pos.reserve(dof*node_num);
            for(size_t i = 0; i < node_num; ++i){
                const auto& n = e->nodes[i];
                for(size_t j = 0; j < dof; ++j){
                    u_pos.push_back(n->u_pos[j]);
                }
            }
            std::vector<double> rhoD(D.size());
            double r = std::pow(*rho, pc);
            for(size_t i = 0; i < rhoD.size(); ++i){
                rhoD[i] = r*D[i] + (1-r)*D2[i];
            }
            std::vector<double> k = e->get_k(rhoD, t);
            this->insert_element_matrix(K, k, u_pos, N);
            ++rho;
        }
    } else {
        // TODO
        logger::log_assert(false, logger::ERROR, "use of more than two materials for topology optimization is not currently supported.");
    }
}

void FiniteElement::insert_element_matrix(std::vector<double>& K, const std::vector<double>& k, const std::vector<long>& pos, const size_t n) const{
    const size_t w = pos.size();
    for(size_t i = 0; i < w; ++i){
        for(size_t j = i; j < w; ++j){
            if(pos[i] > -1 && pos[j] > -1){
                K[utils::to_lower_band(pos[i], pos[j], n)] += k[w*i + j];
            }
        }
    }
}
