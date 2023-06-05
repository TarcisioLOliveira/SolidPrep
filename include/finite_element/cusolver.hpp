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

#ifndef SP_CUSOLVER_HPP
#define SP_CUSOLVER_HPP

#include <memory>
#include "finite_element.hpp"

namespace finite_element::cuda{
    class cuSolverImpl;
}

namespace finite_element{

class cuSolver : public FiniteElement{
    public:
    cuSolver();
    virtual ~cuSolver() = default;
    virtual std::vector<double> calculate_displacements(const Meshing* const mesh, std::vector<double> load, const std::vector<double>& density = std::vector<double>(), double pc = 3, double psi = 0.1) override;

    private:
    std::unique_ptr<cuda::cuSolverImpl> solver;
};

}

#endif
