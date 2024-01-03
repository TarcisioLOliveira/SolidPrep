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

#ifndef SOLVER_MANAGER_HPP
#define SOLVER_MANAGER_HPP

#include "finite_element.hpp"

class SolverManager{
    public:
    SolverManager(std::vector<std::unique_ptr<FiniteElement>> solvers):
        solvers(std::move(solvers)){}

    void generate_matrix(const Meshing* const mesh, const std::vector<double>& density = std::vector<double>(), double pc = 3, double psi = 0.5);

    void calculate_displacements_global(const Meshing* const mesh, std::vector<std::vector<double>>& load, std::vector<double>& u);
    void calculate_displacements_adjoint(const Meshing* const mesh, std::vector<std::vector<double>>& load, std::vector<double>& u);

    const std::vector<std::vector<double>>& sub_u = this->split_u;
    const std::vector<std::vector<double>>& D_vec = this->D_matrices;

    private:
    std::vector<std::vector<double>> split_u;
    std::vector<std::unique_ptr<FiniteElement>> solvers;
    std::vector<double> old_densities;
    std::vector<std::vector<double>> D_matrices;

    void update_D_matrices(const Meshing* const mesh, const std::vector<double>& density = std::vector<double>(), double pc = 3, double psi = 0.5);

    inline bool modified_densities(const std::vector<double>& density, size_t start, size_t range){
        // Using equality to compare densities because the tolerance would
        // have to be below the optimizer's xtol_abs, otherwise it could
        // prematurely stop optimization.
        bool modified = false;
        for(size_t i = 0; i < range; ++i){
            if(old_densities[start + i] != density[start + i]){
                modified = true;
                break;
            }
        }

        return modified;
    }
    inline void update_densities(const std::vector<double>& density, size_t start, size_t range){
        for(size_t i = 0; i < range; ++i){
            old_densities[start + i] = density[start + i];
        }
    }
};

#endif
