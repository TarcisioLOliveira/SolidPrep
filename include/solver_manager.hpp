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

    std::vector<double> calculate_displacements(const Meshing* const mesh, const std::vector<double>& density = std::vector<double>(), double pc = 3, double psi = 0.1);

    std::vector<double> calculate_displacements(const Meshing* const mesh, const std::vector<std::vector<double>>& load, const std::vector<double>& density = std::vector<double>(), double pc = 3, double psi = 0.1);

    inline void set_steps(size_t s){
        for(auto& solv:solvers){
            solv->set_steps(s);
        }
    }

    const std::vector<std::vector<double>>& sub_u = this->split_u;

    private:
    std::vector<std::vector<double>> split_u;
    std::vector<std::unique_ptr<FiniteElement>> solvers;
};

#endif
