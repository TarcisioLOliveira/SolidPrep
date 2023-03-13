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
#include <mpich-x86_64/mpi.h>
#include <numeric>

namespace function{

Volume::Volume(const Meshing* const mesh):
    mesh(mesh){}

void Volume::initialize(const Optimizer* const op){
    auto v = op->get_volumes();
    auto v_it = v.cbegin();
    for(auto& g:this->mesh->geometries){
        if(g->do_topopt){
            this->max_V += std::accumulate(v_it, v_it + g->mesh.size(), 0);
        }
        v_it += g->mesh.size();
    }
}

double Volume::calculate(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x){
    (void)u;
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return 0;
    }
    double V = 0;
    auto v = op->get_volumes();

    auto x_it = x.cbegin();
    auto v_it = v.cbegin();
    for(auto& g:this->mesh->geometries){
        if(g->do_topopt){
            for(auto xi = x_it; xi < x_it+g->mesh.size(); ++xi, ++v_it){
                V += *xi*(*v_it);
            }
            x_it += g->mesh.size();
        } else {
            v_it += g->mesh.size();
        }
    }
    V /= this->max_V;

    return V;
}
double Volume::calculate_with_gradient(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
    (void)u;
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return 0;
    }
    double V = 0;
    auto v = op->get_volumes();

    auto x_it = x.cbegin();
    auto v_it = v.cbegin();
    auto grad_it = grad.begin();
    for(auto& g:this->mesh->geometries){
        if(g->do_topopt){
            for(auto xi = x_it; xi < x_it+g->mesh.size(); ++xi, ++v_it, ++grad_it){
                V += *xi*(*v_it);
                *grad_it = *v_it/this->max_V;
            }
            x_it += g->mesh.size();
        } else {
            v_it += g->mesh.size();
        }
    }
    V /= this->max_V;

    return V;
}

}

