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

#include <mpi.h>
#include "function/node_shape_based/global_stress_heaviside.hpp"
#include "logger.hpp"
#include "optimizer.hpp"

namespace function::node_shape_based{

GlobalStressHeaviside::GlobalStressHeaviside(const Meshing* const mesh, SolverManager* fem, double max_stress, double C):
    mesh(mesh), fem(fem), max_stress(max_stress), C(C){}

double GlobalStressHeaviside::calculate(const NodeShapeBasedOptimizer* const op, const std::vector<double>& u){
    (void) op;
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
    if(mpi_id != 0){
        return 0;
    }

    double result = 0;
    for(auto& g:this->mesh->geometries){
        for(const auto& e:g->mesh){
            const auto c = e->get_centroid();
            const auto D = g->materials.get_D(e.get(), c);
            const double Se = e->get_stress_at(D, c, u, this->vm_eps);
            const double H = this->heaviside(Se);
            result += H*Se;
        }
    }

    return result;
}

double GlobalStressHeaviside::calculate_with_gradient(const NodeShapeBasedOptimizer* const op, const std::vector<double>& u, std::vector<double>& grad){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    std::vector<std::vector<double>> fl(this->mesh->sub_problems->size());
    std::vector<double> l(mesh->max_dofs,0);
    for(size_t i = 0; i < fl.size(); ++i){
        fl[i].resize(mesh->max_dofs, 0);
    }

    auto nodes = op->shape_handler.get_nodes();
    const size_t node_num = this->mesh->elem_info->get_nodes_per_element();
    const size_t kw = this->mesh->elem_info->get_k_dimension();
    const size_t dof = this->mesh->elem_info->get_dof_per_node();
    std::fill(grad.begin(), grad.end(), 0);
    math::Matrix D;

    double result = 0;

    if(mpi_id == 0){
        for(auto& g:this->mesh->geometries){
            for(const auto& e:g->mesh){
                const auto c = e->get_centroid();
                const auto D = g->materials.get_D(e.get(), c);
                const double Se = e->get_stress_at(D, c, u, this->vm_eps);
                const double H = this->heaviside(Se);
                const double dH = this->heaviside_grad(Se);
                result += H*Se;
                e->get_virtual_load(D, dH + H/Se, e->get_centroid(), u, fl[0]);
            }
        }
        logger::quick_log("Calculating adjoint problem...{");
    }
    this->fem->calculate_displacements_adjoint(this->mesh, fl, l);
    if(mpi_id == 0){
        logger::quick_log("} Done.");

        logger::quick_log("Calculating stress gradient...");

        std::vector<double> ue(kw, 0);
        std::vector<double> le(kw, 0);
        std::vector<double> dKu(kw, 0);

        for(size_t i = 0; i < nodes.size(); ++i){
            const auto& shn = nodes[i];
            for(const auto& e:shn.elements){
                const auto& g = this->mesh->elem_geom_mapping.at(e.e);
                const auto c = e.e->get_centroid();
                D = g->materials.get_D(e.e, c);
                const double Se = e.e->get_stress_at(D, c, u, this->vm_eps);
                const double H = this->heaviside(Se);
                const double dH = this->heaviside_grad(Se);
                const double mult = dH + H/Se;
                for(size_t n = 0; n < node_num; ++n){
                    for(size_t j = 0; j < dof; ++j){
                        ue[n*dof + j] = u[e.e->nodes[n]->u_pos[j]];
                        le[n*dof + j] = l[e.e->nodes[n]->u_pos[j]];
                    }
                }
                for(size_t j = 0; j < dof; ++j){
                    const auto dK = e.e->get_dk_sh(D, this->mesh->thickness, e.node_num, j);
                    cblas_dgemv(CblasRowMajor, CblasNoTrans, kw, kw, 1, dK.data(), kw, ue.data(), 1, 0, dKu.data(), 1);
                    const double ldKu = cblas_ddot(kw, le.data(), 1, dKu.data(), 1);

                    const double dS = e.e->von_Mises_derivative_sh(D, mult, c, u, e.node_num, j);


                    grad[i*dof + j] += dS - ldKu;
                }
            }
        }
        logger::quick_log("Done.");
    }

    return result;
}

}
