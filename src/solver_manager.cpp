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

#include "solver_manager.hpp"
#include "logger.hpp"

std::vector<double> SolverManager::calculate_displacements(const Meshing* const mesh, const std::vector<double>& density, double pc, double psi){
    return this->calculate_displacements(mesh, mesh->load_vector, density, pc, psi);
}

std::vector<double> SolverManager::calculate_displacements(const Meshing* const mesh, const std::vector<std::vector<double>>& load, const std::vector<double>& density, double pc, double psi){
    std::vector<double> u(mesh->max_dofs, 0);
    this->split_u.resize(mesh->sub_problems->size());
    for(size_t i = 0; i < mesh->sub_problems->size(); ++i){
        auto& l = load[i];
        auto& n = mesh->node_positions[i];
        this->split_u[i].resize(mesh->max_dofs, 0);
        auto& ui = this->split_u[i];
        if(l.size() == mesh->max_dofs){
            std::vector<double> load_pos(mesh->dofs_per_subproblem[i]);
            #pragma omp parallel for
            for(size_t j = 0; j < n.size(); ++j){
                const long p_new = n[j];
                if(p_new >= 0){
                    load_pos[p_new] = l[j];
                }
            }
            auto u_tmp = this->solvers[i]->calculate_displacements(mesh, n, load_pos, density, pc, psi);
            #pragma omp parallel for
            for(size_t j = 0; j < n.size(); ++j){
                const long p_new = n[j];
                if(p_new >= 0){
                    u[j] += u_tmp[p_new];
                    ui[j] = u_tmp[p_new];
                }
            }
        } else if(l.size() == mesh->dofs_per_subproblem[i]){
            auto u_tmp = this->solvers[i]->calculate_displacements(mesh, n, l, density, pc, psi);
            #pragma omp parallel for
            for(size_t j = 0; j < n.size(); ++j){
                const long p_new = n[j];
                if(p_new >= 0){
                    u[j] += u_tmp[p_new];
                    ui[j] = u_tmp[p_new];
                }
            }
        } else {
            logger::log_assert(false, logger::ERROR, "incorrect load vector size, must be equal to mesh->max_dofs OR mesh->dofs_per_subproblem[i]");
        }
    }
    return u;
}

