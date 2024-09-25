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


#ifndef NONLINEAR_SOLVER_MMA_HPP
#define NONLINEAR_SOLVER_MMA_HPP

#include "nonlinear_solver.hpp"

namespace nonlinear_solver{

class MMA : public NonlinearSolver{
    public:
    MMA(double INIT, double INC, double DEC, double delta_min, double delta_max, double xmin, double xmax, double xtol_abs);

    virtual ~MMA() = default;

    virtual void setup(size_t N);

    virtual bool update(double* x, double f, const double* dfdx);

    virtual SolverClass get_class() const{
        return SolverClass::GRADIENT_BASED;
    }

    virtual SolverType get_type() const{
        return SolverType::MMA;
    }

    private:
    const double INIT, INC, DEC;
    double fold1 = 0, fold2 = 0;
    size_t N;
    size_t it = 0;
    const double delta_min;
    const double delta_max;
    const double xmin;
    const double xmax;
    std::vector<double> xold1, xold2, delta_e;
};

}

#endif
