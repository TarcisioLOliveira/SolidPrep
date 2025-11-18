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

#include <cblas.h>
#include <mpich-x86_64/mpi.h>
#include "function/node_shape_based/compliance.hpp"
#include "math/matrix.hpp"
#include "optimizer.hpp"
#include "project_data.hpp"

namespace function::node_shape_based{

Compliance::Compliance(const projspec::DataMap& data):
    mesh(data.proj->topopt_mesher.get()),
    fem(data.proj->topopt_fea.get())
    {}

double Compliance::calculate(const NodeShapeBasedOptimizer* const op, const std::vector<double>& u){
    (void)op;
    const double c = cblas_ddot(u.size(), this->mesh->global_load_vector.data(), 1, u.data(), 1);

    return c;
}
double Compliance::calculate_with_gradient(const NodeShapeBasedOptimizer* const op, const std::vector<double>& u, std::vector<double>& grad){
    (void)op;
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return 0;
    }

    const size_t dof = this->mesh->elem_info->get_dof_per_node();
    const size_t node_num = this->mesh->elem_info->get_nodes_per_element();
    const size_t kw = this->mesh->elem_info->get_k_dimension();
    math::Matrix D;
    std::fill(grad.begin(), grad.end(), 0);

    std::vector<double> ue(kw, 0);
    std::vector<double> dKu(kw, 0);

    auto nodes = op->shape_handler.get_nodes();

    for(size_t i = 0; i < nodes.size(); ++i){
        const auto& shn = nodes[i];
        for(const auto& e:shn.elements){
            const auto& g = this->mesh->elem_geom_mapping.at(e.e);
            const auto c = e.e->get_centroid();
            D = g->materials.get_D(e.e, c);
            for(size_t n = 0; n < node_num; ++n){
                for(size_t j = 0; j < dof; ++j){
                    ue[n*dof + j] = u[e.e->nodes[n]->u_pos[j]];
                }
            }
            for(size_t j = 0; j < dof; ++j){
                const auto dK = e.e->get_dk_sh(D, this->mesh->thickness, e.node_num, j);
                cblas_dgemv(CblasRowMajor, CblasNoTrans, kw, kw, 1, dK.data(), kw, ue.data(), 1, 0, dKu.data(), 1);
                const double udKu = cblas_ddot(kw, ue.data(), 1, dKu.data(), 1);

                grad[i*dof + j] -= udKu;
            }
        }
    }
    op->contact_gradient(mesh, fem, u, u, grad);
    const double c = cblas_ddot(u.size(), this->mesh->global_load_vector.data(), 1, u.data(), 1);

    return c;
}

using namespace projspec;
const bool Compliance::reg = Factory<NodeShapeBasedFunction>::add(
    [](const DataMap& data){
        return std::make_unique<Compliance>(data);
    },
    ObjectRequirements{
        "compliance",
        {
        }
    }
);

}
