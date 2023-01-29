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

#include "density_filter/convolution.hpp"

namespace density_filter{

Convolution::Convolution(const double radius):
    radius(radius){}

void Convolution::initialize(const Meshing* const mesh, const size_t x_size){
    // Caching positions, because this calculation is expensive.
    this->p.resize(x_size*3,0);
    this->w.resize(x_size,0);
    auto p_it = this->p.begin();
    size_t neighbors_size = 0;
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            neighbors_size += g->mesh.size()*num_den;
            for(const auto& e:g->mesh){
                gp_Pnt c = e->get_centroid();
                *p_it = c.X();
                ++p_it;
                *p_it = c.Y();
                ++p_it;
                *p_it = c.Z();
                ++p_it;
            }
        }
    }
    this->neighbors.resize(neighbors_size);
    size_t offset = 0;
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t elem_num = g->mesh.size();
            const size_t num_den = g->number_of_densities_needed();
            const size_t end = elem_num*num_den + offset;
            for(size_t i = offset; i < end; i += num_den){
                for(size_t j = i; j < end; j += num_den){
                    double dist = this->get_distance(i, j);
                    if(dist <= this->radius){
                        double w = 1 - dist/this->radius;
                        for(size_t k = 0; k < num_den; ++k){
                            this->neighbors[i+k].push_back(j+k);
                            this->neighbors[j+k].push_back(i+k);
                            this->w[i+k] += w;
                            this->w[j+k] += w;
                        }
                    }
                }
            }
            offset += elem_num*num_den;
        }
    }
}
void Convolution::filter_densities(const std::vector<double>& x, std::vector<double>& new_x){
    if(new_x.size() < x.size()){
        new_x.resize(x.size());
    }
    #pragma omp parallel for
    for(size_t i = 0; i < x.size(); ++i){
        new_x[i] = 0;
        for(const auto& j:this->neighbors[i]){
             double dist = this->get_distance(i, j);
             double wj = 1 - dist/this->radius;
             new_x[i] += wj*x[j];
        }
        new_x[i] /= this->w[i];
    }
}

void Convolution::filter_gradient(const std::vector<double>& df, std::vector<double>& new_df){
    if(new_df.size() < df.size()){
        new_df.resize(df.size());
    }
    #pragma omp parallel for
    for(size_t i = 0; i < df.size(); ++i){
        new_df[i] = 0;
        for(const auto& j:this->neighbors[i]){
            double dist = this->get_distance(i, j);
            double wj = 1 - dist/this->radius;
            new_df[i] += wj*df[j]/this->w[j];
        }
    }
}

}
