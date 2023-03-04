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
#include "logger.hpp"
#include "utils.hpp"

namespace density_filter{

void Averaging::initialize(const Meshing* const mesh, const size_t x_size){
    (void)x_size;
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    this->mesh = mesh;
    const double t = this->mesh->thickness;

    size_t geom_id = 0;
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            // this->id_mapping.reserve(this->id_mapping.size()+g->mesh.size());
            for(const auto& e:g->mesh){
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n = e->nodes[i];
                    this->id_mapping[std::make_pair(geom_id, n->id)] = num_den;
                }
            }
        }
        ++geom_id;
    }
    size_t id = 0;
    for(auto it = this->id_mapping.begin(); it != this->id_mapping.end(); ++it){
        const size_t num_den = it->second;
        it->second = id;
        id += num_den;
    }
    size_t len = id;
    geom_id = 0;
    this->nodal_densities.resize(len,0);
    this->nodal_gradient.resize(len,0);
    this->D.resize(len,0);
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                const double V = e->get_volume(t)/num_nodes;
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n = e->nodes[i];
                    for(size_t j = 0; j < num_den; ++j){
                        this->D[this->id_mapping[std::make_pair(geom_id, n->id)]+j] += V;
                    }
                }
            }
        }
        ++geom_id;
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
    size_t geom_id = 0;
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                const double V = e->get_volume(t)/num_nodes;
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n = e->nodes[i];
                    for(size_t j = 0; j < num_den; ++j){
                        this->nodal_densities[this->id_mapping[std::make_pair(geom_id, n->id)]+j] += V*(*(x_it+j));
                    }
                }
                x_it += num_den;
            }
        }
        ++geom_id;
    }
    for(size_t i = 0; i < this->nodal_densities.size(); ++i){
        this->nodal_densities[i] /= this->D[i];
    }
    geom_id = 0;
    auto newx_it = new_x.begin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                *newx_it = 0;
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n = e->nodes[i];
                    for(size_t j = 0; j < num_den; ++j){
                        *(newx_it + j) += this->nodal_densities[this->id_mapping[std::make_pair(geom_id, n->id)]+j];
                    }
                }
                for(size_t j = 0; j < num_den; ++j){
                    *newx_it /= num_nodes;
                    ++newx_it;
                }
            }
        }
        ++geom_id;
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
    size_t geom_id = 0;
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                const double V = e->get_volume(t)/num_nodes;
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n = e->nodes[i];
                    for(size_t j = 0; j < num_den; ++j){
                        this->nodal_gradient[this->id_mapping[std::make_pair(geom_id, n->id)]+j] += V*(*(x_it+j));
                    }
                }
                x_it += num_den;
            }
        }
        ++geom_id;
    }
    for(size_t i = 0; i < this->nodal_gradient.size(); ++i){
        this->nodal_gradient[i] /= this->D[i];
    }
    geom_id = 0;
    auto newx_it = new_df.begin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                *newx_it = 0;
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n = e->nodes[i];
                    for(size_t j = 0; j < num_den; ++j){
                        *(newx_it+j) += this->nodal_gradient[this->id_mapping[std::make_pair(geom_id, n->id)]+j];
                    }
                }
                for(size_t j = 0; j < num_den; ++j){
                    *newx_it /= num_nodes;
                    ++newx_it;
                }
            }
        }
        ++geom_id;
    }
}

void Averaging::filter_gradient_nodal(const std::vector<double>& df, std::vector<double>& new_df){
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    const double t = this->mesh->thickness;
    if(new_df.size() < df.size()){
        new_df.resize(df.size());
    }
    std::copy(df.begin(), df.end(), this->nodal_gradient.begin());
    for(size_t i = 0; i < this->nodal_gradient.size(); ++i){
        this->nodal_gradient[i] /= this->D[i];
    }
    size_t geom_id = 0;
    auto newx_it = new_df.begin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                *newx_it = 0;
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n = e->nodes[i];
                    for(size_t j = 0; j < num_den; ++j){
                        *(newx_it+j) += this->nodal_gradient[this->id_mapping[std::make_pair(geom_id, n->id)]+j];
                    }
                }
                for(size_t j = 0; j < num_den; ++j){
                    *newx_it /= num_nodes;
                    ++newx_it;
                }
            }
        }
        ++geom_id;
    }
}

void Averaging::get_gradient(std::vector<double>& gradx) const{
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    const double t = this->mesh->thickness;
    size_t N = 0;
    if(this->mesh->elem_info->get_problem_type() == utils::PROBLEM_TYPE_2D){
        N = 2;
    } else if(this->mesh->elem_info->get_problem_type() == utils::PROBLEM_TYPE_3D){
        N = 3;
    }
    size_t geom_id = 0;
    std::fill(gradx.begin(), gradx.end(), 0);
    std::vector<double> nx(num_nodes,0);
    auto g_it = gradx.begin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                auto grad = e->get_nodal_density_gradient(e->get_centroid());
                //const double V = e->get_volume(t)/num_nodes;
                for(size_t j = 0; j < num_den; ++j){
                    for(size_t i = 0; i < num_nodes; ++i){
                        const auto& n = e->nodes[i];
                        //nx[i] = V*this->nodal_densities[this->id_mapping.at(std::make_pair(geom_id, n->id))+j];
                        nx[i] = this->nodal_densities[this->id_mapping.at(std::make_pair(geom_id, n->id))+j];
                    }
                    for(size_t i = 0; i < N; ++i){
                        for(size_t k = 0; k < num_nodes; ++k){
                            *g_it += grad[i*num_nodes + k]*nx[k];
                        }
                        ++g_it;
                    }
                }
            }
        }
        ++geom_id;
    }
}

}
