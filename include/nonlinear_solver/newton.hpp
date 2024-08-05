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


#ifndef NONLINEAR_SOLVER_NEWTON_HPP
#define NONLINEAR_SOLVER_NEWTON_HPP

#include "nonlinear_solver.hpp"

namespace nonlinear_solver{

class Newton : public NonlinearSolver{
    public:
    Newton(double MULT, double rtol_abs);

    virtual ~Newton() = default;

    virtual void setup(size_t N);

    virtual bool update(double* x, double f, const double* dfdx);

    virtual SolverClass get_class() const{
        return SolverClass::HESSIAN_BASED;
    }

    virtual SolverType get_type() const{
        return SolverType::NEWTON;
    }

    private:
    const double MULT;
    size_t N;
    size_t it = 0;
};

}

#endif
