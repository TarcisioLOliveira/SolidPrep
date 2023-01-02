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

#include "function/volume.hpp"

namespace function{

Volume::Volume(const Meshing* const mesh):
    mesh(mesh){}

void Volume::initialize(){
    size_t x_size = 0;
    for(const auto& g:this->mesh->geometries){
        if(g->do_topopt){
            x_size += g->mesh.size();
        }
    }

    this->grad_V.resize(x_size, 0);

    auto V_it = this->grad_V.begin();
    for(const auto& g:this->mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                *V_it = e->get_volume(mesh->thickness);
                this->max_V += *V_it;
                ++V_it;
            }
        }
    }
}

double Volume::calculate(const std::vector<double>& u, const std::vector<double>& x){
    (void)u;
    double V = 0;

    #pragma omp parallel for reduction(+:V)
    for(size_t i = 0; i < x.size(); ++i){
        V += x[i]*this->grad_V[i]/this->max_V;
    }

    return V;
}
double Volume::calculate_with_gradient(const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
    (void)u;
    double V = 0;

    #pragma omp parallel for reduction(+:V)
    for(size_t i = 0; i < x.size(); ++i){
        grad[i] = this->grad_V[i]/this->max_V;
        V += x[i]*this->grad_V[i]/this->max_V;
    }

    return V;
}

}

