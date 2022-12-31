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

#include "density_filter/averaging.hpp"

namespace density_filter{

void Averaging::initialize(const Meshing* const mesh, const size_t x_size){
    (void)x_size;
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    this->mesh = mesh;
    const double t = this->mesh->thickness;

    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            this->id_mapping.reserve(this->id_mapping.size()+g->mesh.size());
            for(const auto& e:g->mesh){
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n = e->nodes[i];
                    this->id_mapping[n->id] = 0;
                }
            }
        }
    }
    size_t len = this->id_mapping.size();
    size_t id = 0;
    for(auto it = this->id_mapping.begin(); it != this->id_mapping.end(); ++it){
        it->second = id;
        ++id;
    }
    this->nodal_densities.resize(len,0);
    this->nodal_gradient.resize(len,0);
    this->D.resize(len,0);
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                const double V = e->get_volume(t)/num_nodes;
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n = e->nodes[i];
                    this->D[this->id_mapping[n->id]] += V;
                }
            }
        }
    }
}

void Averaging::filter_densities(const std::vector<double>& x, std::vector<double>& new_x){
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    const double t = this->mesh->thickness;
    if(new_x.size() < x.size()){
        new_x.resize(x.size());
    }
    std::fill(this->nodal_densities.begin(), this->nodal_densities.end(), 0);
    auto x_it = x.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                const double V = e->get_volume(t)/num_nodes;
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n = e->nodes[i];
                    this->nodal_densities[this->id_mapping[n->id]] += V*(*x_it);
                }
                ++x_it;
            }
        }
    }
    for(size_t i = 0; i < this->nodal_densities.size(); ++i){
        this->nodal_densities[i] /= this->D[i];
    }
    auto newx_it = new_x.begin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                *newx_it = 0;
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n = e->nodes[i];
                    *newx_it += this->nodal_densities[this->id_mapping[n->id]];
                }
                *newx_it /= num_nodes;
                ++newx_it;
            }
        }
    }
}


void Averaging::filter_gradient(const std::vector<double>& df, std::vector<double>& new_df){
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    const double t = this->mesh->thickness;
    if(new_df.size() < df.size()){
        new_df.resize(df.size());
    }
    std::fill(this->nodal_gradient.begin(), this->nodal_gradient.end(), 0);
    auto x_it = df.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                const double V = e->get_volume(t)/num_nodes;
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n = e->nodes[i];
                    this->nodal_gradient[this->id_mapping[n->id]] += V*(*x_it);
                }
                ++x_it;
            }
        }
    }
    for(size_t i = 0; i < this->nodal_gradient.size(); ++i){
        this->nodal_gradient[i] /= this->D[i];
    }
    auto newx_it = new_df.begin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                *newx_it = 0;
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n = e->nodes[i];
                    *newx_it += this->nodal_gradient[this->id_mapping[n->id]];
                }
                *newx_it /= num_nodes;
                ++newx_it;
            }
        }
    }
}


}
