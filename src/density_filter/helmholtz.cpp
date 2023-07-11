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

#include <lapacke.h>
#include <set>
#include "density_filter/helmholtz.hpp"
#include "logger.hpp"

namespace density_filter{

Helmholtz::Helmholtz(const double radius)
    : radius(radius/(2*std::sqrt(3))){}

void Helmholtz::initialize(const Meshing* const mesh, const size_t x_size){
    (void)x_size;
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    this->mesh = mesh;
    const double t = this->mesh->thickness;

    this->NN_kd = 0;
    size_t geom_id = 0;
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            // this->id_mapping.reserve(this->id_mapping.size()+g->mesh.size());
            for(const auto& e:g->mesh){
                size_t top = 0;
                size_t bot = std::numeric_limits<size_t>::max();
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n = e->nodes[i];
                    if(n->id > top){
                        top = n->id;
                    }
                    if(n->id < bot){
                        bot = n->id;
                    }
                    this->id_mapping[std::make_pair(geom_id, n->id)] = num_den;
                }
                if((top - bot + 1) > NN_kd){
                    NN_kd = top - bot + 1;
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
    this->NN_n = len;
    this->nodal_densities.resize(this->NN_n,0);
    this->nodal_gradient.resize(this->NN_n,0);
    this->NN = std::vector<double>(NN_n*NN_kd,0);
    this->id_mapping_linear.resize(x_size*num_nodes);
    auto id_it = id_mapping_linear.begin();
    const std::vector<double> I{1, 0, 0,
                                0, 1, 0,
                                0, 0, 1};
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                const auto M = radius*radius*e->diffusion_1dof(t, I) +
                               e->absorption_1dof(t);
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n1 = e->nodes[i];
                    *id_it = this->id_mapping[std::make_pair(geom_id, n1->id)];
                    const size_t id1 = *id_it;
                    for(size_t k = 0; k <= i; ++k){
                        const auto& n2 = e->nodes[k];
                        const size_t id2 = this->id_mapping[std::make_pair(geom_id, n2->id)];
                        for(size_t j = 0; j < num_den; ++j){
                            this->NN[utils::to_band(id1, id2, NN_kd)] += M(i, k);
                        }
                    }
                    ++id_it;
                }
            }
        }
        ++geom_id;
    }
    int info = LAPACKE_dpbtrf_work(LAPACK_COL_MAJOR, 'L', NN_n, NN_kd-1, NN.data(), NN_kd);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while factoring Helmholtz matrix.", info);
    this->id_mapping.clear();
}

void Helmholtz::filter_densities(const std::vector<double>& x, std::vector<double>& new_x){
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    const double t = this->mesh->thickness;
    if(new_x.size() < x.size()){
        new_x.resize(x.size());
    }
    std::fill(this->nodal_densities.begin(), this->nodal_densities.end(), 0);
    auto x_it = x.cbegin();
    auto id_it = id_mapping_linear.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                auto N = e->source_1dof(t);
                for(size_t i = 0; i < num_nodes; ++i){
                    for(size_t j = 0; j < num_den; ++j){
                        this->nodal_densities[*id_it+j] += N[i]*(*(x_it+j));
                    }
                    ++id_it;
                }
                x_it += num_den;
            }
        }
    }
    int info = LAPACKE_dpbtrs_work(LAPACK_COL_MAJOR, 'L', NN_n, NN_kd-1, 1, NN.data(), NN_kd, this->nodal_densities.data(), NN_n);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating nodal densities.", info);
    auto newx_it = new_x.begin();
    id_it = id_mapping_linear.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(size_t it = 0; it < g->mesh.size(); ++it){
                *newx_it = 0;
                for(size_t i = 0; i < num_nodes; ++i){
                    for(size_t j = 0; j < num_den; ++j){
                        *(newx_it + j) += this->nodal_densities[*id_it+j];
                    }
                    ++id_it;
                }
                for(size_t j = 0; j < num_den; ++j){
                    *newx_it /= num_nodes;
                    ++newx_it;
                }
            }
        }
    }
}


void Helmholtz::filter_gradient(const std::vector<double>& df, std::vector<double>& new_df){
    if(new_df.size() < df.size()){
        new_df.resize(df.size());
    }
    const double t = this->mesh->thickness;
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    std::fill(this->nodal_gradient.begin(), this->nodal_gradient.end(), 0);
    auto x_it = df.cbegin();
    auto id_it = id_mapping_linear.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                auto N = e->source_1dof(t);
                for(size_t i = 0; i < num_nodes; ++i){
                    for(size_t j = 0; j < num_den; ++j){
                        this->nodal_gradient[*id_it+j] += N[i]*(*(x_it+j));
                    }
                    ++id_it;
                }
                x_it += num_den;
            }
        }
    }
    int info = LAPACKE_dpbtrs_work(LAPACK_COL_MAJOR, 'L', NN_n, NN_kd-1, 1, NN.data(), NN_kd, this->nodal_gradient.data(), NN_n);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating nodal densities' gradient.", info);
    auto newx_it = new_df.begin();
    id_it = id_mapping_linear.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(size_t it = 0; it < g->mesh.size(); ++it){
                *newx_it = 0;
                for(size_t i = 0; i < num_nodes; ++i){
                    for(size_t j = 0; j < num_den; ++j){
                        *(newx_it+j) += this->nodal_gradient[*id_it+j];
                    }
                    ++id_it;
                }
                for(size_t j = 0; j < num_den; ++j){
                    *newx_it /= num_nodes;
                    ++newx_it;
                }
            }
        }
    }
}


void Helmholtz::filter_gradient_nodal(const std::vector<double>& df, std::vector<double>& new_df){
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    if(new_df.size() < df.size()){
        new_df.resize(df.size());
    }
    std::copy(df.begin(), df.end(), this->nodal_gradient.begin());
    int info = LAPACKE_dpbtrs_work(LAPACK_COL_MAJOR, 'L', NN_n, NN_kd-1, 1, NN.data(), NN_kd, this->nodal_gradient.data(), NN_n);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating nodal densities' gradient.", info);
    auto newx_it = new_df.begin();
    auto id_it = id_mapping_linear.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(size_t it = 0; it < g->mesh.size(); ++it){
                *newx_it = 0;
                for(size_t i = 0; i < num_nodes; ++i){
                    for(size_t j = 0; j < num_den; ++j){
                        *(newx_it+j) += this->nodal_gradient[*id_it+j];
                    }
                    ++id_it;
                }
                for(size_t j = 0; j < num_den; ++j){
                    *newx_it /= num_nodes;
                    ++newx_it;
                }
            }
        }
    }
}

void Helmholtz::get_gradient(std::vector<double>& gradx) const{
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    size_t N = 0;
    if(this->mesh->elem_info->get_problem_type() == utils::PROBLEM_TYPE_2D){
        N = 2;
    } else if(this->mesh->elem_info->get_problem_type() == utils::PROBLEM_TYPE_3D){
        N = 3;
    }
    std::fill(gradx.begin(), gradx.end(), 0);
    std::vector<double> nx(num_nodes,0);
    auto g_it = gradx.begin();
    auto id_it = id_mapping_linear.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                auto grad = e->get_nodal_density_gradient(e->get_centroid());
                for(size_t j = 0; j < num_den; ++j){
                    for(size_t i = 0; i < num_nodes; ++i){
                        nx[i] = this->nodal_densities[*(id_it+i)+j];
                    }
                    for(size_t i = 0; i < N; ++i){
                        for(size_t k = 0; k < num_nodes; ++k){
                            *g_it += grad[i*num_nodes + k]*nx[k];
                        }
                        ++g_it;
                    }
                }
                id_it += num_nodes;
            }
        }
    }
}

}
