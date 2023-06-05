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

#include "finite_element/cusolver.hpp"
#include "finite_element/cusolver_impl.hpp"
#include "logger.hpp"

namespace finite_element{

cuSolver::cuSolver():
    solver(std::make_unique<cuda::cuSolverImpl>()){
}

std::vector<double> cuSolver::calculate_displacements(const Meshing* const mesh, std::vector<double> load, const std::vector<double>& density, double pc, double psi){
    return this->solver->calculate_displacements(mesh, std::move(load), density, pc, psi);
}

}
