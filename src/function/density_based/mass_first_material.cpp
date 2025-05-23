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
#include "function/density_based/mass_first_material.hpp"
#include "logger.hpp"
#include "optimizer.hpp"
#include "project_data.hpp"

namespace function::density_based{

MassFirstMaterial::MassFirstMaterial(const projspec::DataMap& data):
    mesh(data.proj->topopt_mesher.get()){}

void MassFirstMaterial::initialize(const DensityBasedOptimizer* const op){
    (void)op;
    for(auto& g:this->mesh->geometries){
        if(g->do_topopt){
            if(g->with_void && g->number_of_materials() > 1){
                logger::log_assert(false, logger::ERROR, "MassFirstMaterial does not support void when using more than 1 material");
            }
        }
    }
}

double MassFirstMaterial::calculate(const DensityBasedOptimizer* const op, const std::vector<double>& u, const std::vector<double>& x){
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
            if(g->with_void){
                for(const auto& e:g->mesh){
                    const gp_Pnt p = e->get_centroid();
                    const double d = g->materials.get_materials()[0]->get_density(e.get(), p);
                    if(num_den == 1){
                        V += (*x_it)*d*(*v_it);
                    } else {
                        V += (*x_it)*(*(x_it+1))*d*(*v_it);
                    }
                    ++v_it;
                    x_it += num_den;
                }
            } else {
                for(const auto& e:g->mesh){
                    const gp_Pnt p = e->get_centroid();
                    const double d = g->materials.get_materials()[0]->get_density(e.get(), p);
                    V += (*x_it)*d*(*v_it);
                    ++v_it;
                    x_it += num_den;
                }
            }
        } else {
            v_it += g->mesh.size();
        }
    }
    V *= 1e3;
    V /= 1e9;

    return V;
}
double MassFirstMaterial::calculate_with_gradient(const DensityBasedOptimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
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
            if(g->with_void){
                for(const auto& e:g->mesh){
                    const gp_Pnt p = e->get_centroid();
                    const double d = g->materials.get_materials()[0]->get_density(e.get(), p);
                    if(num_den == 1){
                        V += (*x_it)*d*(*v_it);
                        *grad_it = d*(*v_it)*1e3/1e9;
                    } else {
                        V += (*x_it)*(*(x_it+1))*d*(*v_it);
                        *grad_it = (*(x_it+1))*d*(*v_it)*1e3/1e9;
                        *(grad_it+1) = (*x_it)*d*(*v_it)*1e3/1e9;
                    }
                    ++v_it;
                    x_it += num_den;
                    grad_it += num_den;
                }
            } else {
                for(const auto& e:g->mesh){
                    const gp_Pnt p = e->get_centroid();
                    const double d = g->materials.get_materials()[0]->get_density(e.get(), p);
                    V += (*x_it)*d*(*v_it);
                    *grad_it = d*(*v_it)*1e3/1e9;
                    ++v_it;
                    x_it += num_den;
                    grad_it += num_den;
                }
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
const bool MassFirstMaterial::reg = Factory<DensityBasedFunction>::add(
    [](const DataMap& data){
        return std::make_unique<MassFirstMaterial>(data);
    },
    ObjectRequirements{
        "mass_first_material",
        {
        }
    }
);

}

