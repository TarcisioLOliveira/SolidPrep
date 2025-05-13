/*
 *   Copyright (C) 2024 Tarcísio Ladeia de Oliveira.
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
#include <numeric>
#include "function/node_shape_based/volume.hpp"
#include "optimizer.hpp"
#include "project_specification/data_map.hpp"
#include "project_data.hpp"

namespace function::node_shape_based{

Volume::Volume(const projspec::DataMap& data):
    mesh(data.proj->topopt_mesher.get()){}

void Volume::initialize(const NodeShapeBasedOptimizer* const op){
    auto v = op->get_volumes();
    auto v_it = v.cbegin();
    for(auto& g:this->mesh->geometries){
        this->max_V += std::accumulate(v_it, v_it + g->mesh.size(), 0.0);
    }
}

double Volume::calculate(const NodeShapeBasedOptimizer* const op, const std::vector<double>& u){
    (void)u;
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return 0;
    }
    double V = 0;
    auto v = op->get_volumes();

    auto v_it = v.cbegin();
    for(auto& g:this->mesh->geometries){
        for(auto vi = v_it; vi < v_it+g->mesh.size(); ++vi){
            V += *vi;
        }
    }
    V /= this->max_V;

    return V;
}
double Volume::calculate_with_gradient(const NodeShapeBasedOptimizer* const op, const std::vector<double>& u, std::vector<double>& grad){
    (void)u;
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return 0;
    }

    const size_t dof = this->mesh->elem_info->get_dof_per_node();
    std::fill(grad.begin(), grad.end(), 0);

    auto nodes = op->shape_handler.get_nodes();

    for(size_t i = 0; i < nodes.size(); ++i){
        const auto& shn = nodes[i];
        for(const auto& e:shn.elements){
            for(size_t j = 0; j < dof; ++j){
                grad[i*dof + j] += e.e->get_dV_sh(this->mesh->thickness, e.node_num, j);
            }
        }
    }
    for(auto& gi:grad){
        gi /= this->max_V;
    }

    return this->calculate(op, u);
}

using namespace projspec;
const bool Volume::reg = Factory<NodeShapeBasedFunction>::add(
    [](const DataMap& data){
        return std::make_unique<Volume>(data);
    },
    ObjectRequirements{
        "volume",
        {
        }
    }
);

}

