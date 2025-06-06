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

#include <mpich-x86_64/mpi.h>
#include "function/density_based/mass.hpp"
#include "optimizer.hpp"
#include "project_data.hpp"

namespace function::density_based{

Mass::Mass(const projspec::DataMap& data):
    mesh(data.proj->topopt_mesher.get()){}

void Mass::initialize(const DensityBasedOptimizer* const op){
    (void) op;
}

double Mass::calculate(const DensityBasedOptimizer* const op, const std::vector<double>& u, const std::vector<double>& x){
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
        const size_t num_den = g->number_of_densities_needed();
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                const gp_Pnt p = e->get_centroid();
                const double d = g->materials.get_density(x_it, e.get(), p);
                V += d*(*v_it);
                ++v_it;
                x_it += num_den;
            }
        } else {
            v_it += g->mesh.size();
        }
    }
    V *= 1e3;
    V /= 1e9;

    return V;
}
double Mass::calculate_with_gradient(const DensityBasedOptimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
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
        const size_t num_den = g->number_of_densities_needed();
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                auto grad_it2 = grad_it;
                const gp_Pnt p = e->get_centroid();
                const auto d = g->materials.get_density_deriv(x_it, e.get(), p, grad_it2);
                V += d*(*v_it);
                while(grad_it < grad_it2){
                    *grad_it *= *v_it*1e3/1e9;
                    ++grad_it;
                }
                ++v_it;
                x_it += num_den;
            }
        } else {
            v_it += g->mesh.size();
        }
    }
    V *= 1e3;
    V /= 1e9;

    return V;
}

using namespace projspec;
const bool Mass::reg = Factory<DensityBasedFunction>::add(
    [](const DataMap& data){
        return std::make_unique<Mass>(data);
    },
    ObjectRequirements{
        "mass",
        {
        }
    }
);

}

