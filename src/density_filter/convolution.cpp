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
    this->neighbors.resize(x_size);
    this->p.resize(x_size*3,0);
    this->w.resize(x_size,0);
    auto p_it = this->p.begin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
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
    for(size_t i = 0; i < x_size; ++i){
        for(size_t j = i; j < x_size; ++j){
            double dist = this->get_distance(i, j);
            if(dist <= this->radius){
                this->neighbors[i].push_back(j);
                this->neighbors[j].push_back(i);
                double w = 1 - dist/this->radius;
                this->w[i] += w;
                this->w[j] += w;
            }
        }
    }
}
std::vector<double> Convolution::filter_densities(const std::vector<double>& x) const{
    std::vector<double> new_x(x.size(),0);

    #pragma omp parallel for
    for(size_t i = 0; i < x.size(); ++i){
        for(const auto& j:this->neighbors[i]){
             double dist = this->get_distance(i, j);
             double wj = 1 - dist/this->radius;
             new_x[i] += wj*x[j];
        }
        new_x[i] /= this->w[i];
    }

    return new_x;
}

std::vector<double> Convolution::filter_gradient(const std::vector<double>& df) const{
    std::vector<double> grad(df.size(),0);

    #pragma omp parallel for
    for(size_t i = 0; i < df.size(); ++i){
         for(const auto& j:this->neighbors[i]){
             double dist = this->get_distance(i, j);
             double wj = 1 - dist/this->radius;
             grad[i] += wj*df[j]/this->w[j];
         }
    }
    return grad;
}

}
