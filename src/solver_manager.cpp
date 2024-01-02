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

void SolverManager::generate_matrix(const Meshing* const mesh, const std::vector<double>& density, double pc, double psi){
    for(size_t i = 0; i < mesh->sub_problems->size(); ++i){
        auto& n = mesh->node_positions[i];
        auto& l = mesh->load_vector[i];
        this->solvers[i]->generate_matrix(mesh, l.size(), n, density, pc, psi);
    }
}

void SolverManager::calculate_displacements_global(const Meshing* const mesh, std::vector<std::vector<double>>& load, std::vector<double>& u){
    u.resize(mesh->max_dofs, 0);
    std::fill(u.begin(), u.end(), 0);
    this->split_u.resize(mesh->sub_problems->size());
    for(size_t i = 0; i < mesh->sub_problems->size(); ++i){
        auto& l = load[i];
        auto& n = mesh->node_positions[i];
        this->split_u[i].resize(mesh->max_dofs, 0);
        auto& ui = this->split_u[i];
        this->solvers[i]->calculate_displacements(l);
        #pragma omp parallel for
        for(size_t j = 0; j < n.size(); ++j){
            const long p_new = n[j];
            if(p_new >= 0){
                u[j] += l[p_new];
                ui[j] = l[p_new];
            }
        }
    }
}

void SolverManager::calculate_displacements_adjoint(const Meshing* const mesh, std::vector<std::vector<double>>& load, std::vector<double>& u){
    for(size_t i = 0; i < mesh->sub_problems->size(); ++i){
        auto& l = load[i];
        auto& n = mesh->node_positions[i];
        this->split_u[i].resize(mesh->max_dofs, 0);
        if(l.size() == mesh->max_dofs){
            std::vector<double> load_pos(mesh->dofs_per_subproblem[i]);
            #pragma omp parallel for
            for(size_t j = 0; j < n.size(); ++j){
                const long p_new = n[j];
                if(p_new >= 0){
                    load_pos[p_new] = l[j];
                }
            }
            this->solvers[i]->calculate_displacements(load_pos);
            #pragma omp parallel for
            for(size_t j = 0; j < n.size(); ++j){
                const long p_new = n[j];
                if(p_new >= 0){
                    u[j] += load_pos[p_new];
                    l[p_new] = load_pos[p_new];
                }
            }
        } else if(l.size() == mesh->dofs_per_subproblem[i]){
            this->solvers[i]->calculate_displacements(l);
            #pragma omp parallel for
            for(size_t j = 0; j < n.size(); ++j){
                const long p_new = n[j];
                if(p_new >= 0){
                    u[j] += l[p_new];
                }
            }
        } else {
            logger::log_assert(false, logger::ERROR, "incorrect load vector size, must be equal to mesh->max_dofs OR mesh->dofs_per_subproblem[i]");
        }
    }
}

