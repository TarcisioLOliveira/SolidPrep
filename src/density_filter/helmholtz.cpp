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

    std::set<size_t> involved_nodes;
    this->NN_kd = 0;
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                size_t top = 0;
                size_t bot = std::numeric_limits<size_t>::max();
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n = e->nodes[i];
                    involved_nodes.emplace(n->id);
                    if(n->id > top){
                        top = n->id;
                    }
                    if(n->id < bot){
                        bot = n->id;
                    }
                }
                if((top - bot + 1) > NN_kd){
                    NN_kd = top - bot + 1;
                }
            }
        }
    }
    this->NN_n = involved_nodes.size();
    involved_nodes.clear();
    this->nodal_densities.resize(this->NN_n,0);
    this->nodal_gradient.resize(this->NN_n,0);
    this->NN = std::vector<double>(NN_n*NN_kd,0);
    const std::vector<double> I{1, 0, 0,
                                0, 1, 0,
                                0, 0, 1};
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                const auto M = radius*radius*e->diffusion_1dof(t, I) +
                               e->absorption_1dof(t);
                for(size_t i = 0; i < num_nodes; ++i){
                    size_t n1 = e->nodes[i]->id;
                    for(size_t j = i; j < num_nodes; ++j){
                        size_t n2 = e->nodes[j]->id;
                        this->NN[utils::to_band(n1, n2, NN_kd)] += M(i, j);
                    }
                }
            }
        }
    }

    int info = LAPACKE_dpbtrf_work(LAPACK_COL_MAJOR, 'L', NN_n, NN_kd-1, NN.data(), NN_kd);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while factoring N^T*N matrix.", info);
}

void Helmholtz::filter_densities(const std::vector<double>& x, std::vector<double>& new_x){
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
                auto N = e->source_1dof(t);
                for(size_t j = 0; j < num_nodes; ++j){
                    size_t id = e->nodes[j]->id;
                    this->nodal_densities[id] += *x_it*N[j];
                }
                ++x_it;
            }
        }
    }

    int info = LAPACKE_dpbtrs_work(LAPACK_COL_MAJOR, 'L', NN_n, NN_kd-1, 1, NN.data(), NN_kd, this->nodal_densities.data(), NN_n);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating nodal densities.", info);
    auto newx_it = new_x.begin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                *newx_it = 0;
                for(size_t j = 0; j < num_nodes; ++j){
                    size_t id = e->nodes[j]->id;
                    *newx_it += this->nodal_densities[id];
                }
                *newx_it /= num_nodes;
                ++newx_it;
            }
        }
    }
}


void Helmholtz::filter_gradient(const std::vector<double>& df, std::vector<double>& new_df){
    if(new_df.size() < df.size()){
        new_df.resize(df.size());
    }
    const double t = this->mesh->thickness;
    std::fill(this->nodal_gradient.begin(), this->nodal_gradient.end(), 0);
    // This works only for a linear interpolation matrix, so try to use that.
    // The alternative is a bit too expensive computationally and in terms of
    // implementation
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    auto df_it = df.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                auto N = e->source_1dof(t);
                for(size_t j = 0; j < num_nodes; ++j){
                    size_t id = e->nodes[j]->id;
                    this->nodal_gradient[id] += *df_it*N[j];
                }
                ++df_it;
            }
        }
    }
    int info = LAPACKE_dpbtrs_work(LAPACK_COL_MAJOR, 'L', NN_n, NN_kd-1, 1, NN.data(), NN_kd, this->nodal_gradient.data(), NN_n);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating nodal densities' gradient.", info);
    auto newdf_it = new_df.begin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                *newdf_it = 0;
                for(size_t j = 0; j < num_nodes; ++j){
                    size_t id = e->nodes[j]->id;
                    *newdf_it += this->nodal_gradient[id];
                }
                *newdf_it /= num_nodes;
                ++newdf_it;
            }
        }
    }
}

}
