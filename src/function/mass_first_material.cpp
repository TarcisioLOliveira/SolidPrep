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

#include "function/mass_first_material.hpp"
#include "logger.hpp"
#include <mpich-x86_64/mpi.h>
#include <numeric>

namespace function{

MassFirstMaterial::MassFirstMaterial(const Meshing* const mesh):
    mesh(mesh){}

void MassFirstMaterial::initialize(const Optimizer* const op){
    (void)op;
    for(auto& g:this->mesh->geometries){
        if(g->do_topopt){
            if(g->with_void && g->number_of_materials() > 1){
                logger::log_assert(false, logger::ERROR, "MassFirstMaterial does not support void when using more than 1 material");
            }
        }
    }
}

double MassFirstMaterial::calculate(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x){
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
            const size_t num_den = g->number_of_densities_needed();
            const double d = g->materials.get_materials()[0]->density;
            for(auto xi = x_it; xi < x_it+g->mesh.size(); xi += num_den, ++v_it){
                V += (*xi)*d*(*v_it);
            }
            x_it += g->mesh.size()*num_den;
        } else {
            v_it += g->mesh.size();
        }
    }
    V /= 1e9;

    return V;
}
double MassFirstMaterial::calculate_with_gradient(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
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
            const size_t num_den = g->number_of_densities_needed();
            const double d = g->materials.get_materials()[0]->density;
            for(auto xi = x_it; xi < x_it+g->mesh.size(); xi += num_den, ++v_it, grad_it += num_den){
                V += (*xi)*d*(*v_it);
                *grad_it = d*(*v_it);
            }
            x_it += g->mesh.size()*num_den;
        } else {
            v_it += g->mesh.size();
        }
    }
    V /= 1e9;

    return V;
}

}

