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


#ifndef NONLINEAR_SOLVER_HPP
#define NONLINEAR_SOLVER_HPP

#include <cstddef>
#include "logger.hpp"

class NonlinearSolver{
    public:
    enum class SolverClass{
        LINEAR,
        GRADIENT_BASED,
        HESSIAN_BASED
    };
    enum class SolverType{
        LINEAR,
        STEEPEST_DESCENT,
        MMA,
        NEWTON
    };

    NonlinearSolver(double xtol_abs, double rtol_abs):
        xtol_abs(xtol_abs), rtol_abs(rtol_abs){}

    const double xtol_abs;
    const double rtol_abs;

    virtual ~NonlinearSolver() = default;

    virtual void setup(size_t N) = 0;

    /**
     * LINEAR is a stub object, this should not be called.
     * GRADIENT_BASED takes all expected parameters and updates x based on
     * gradient.
     * HESSIAN_BASED takes the solution to the linear problem of a iteration
     * in dfdx and updates x based on it, f, and internal parameters.
     *
     * Return true when it reaches stop criterium.
     */
    virtual bool update(double* x, double f, const double* dfdx) = 0;

    virtual SolverClass get_class() const = 0; 

    virtual SolverType get_type() const = 0; 
};

namespace nonlinear_solver{

// Stub class
class Linear : public NonlinearSolver{
    public:
    ~Linear() = default;

    Linear():NonlinearSolver(0,0){}

    virtual void setup(size_t N){
        (void)N;
        logger::log_assert(false, logger::ERROR, "Linear solver does not support nonlinear solver methods");
    }

    virtual bool update(double* x, double f, const double* dfdx){
        (void)x;
        (void)f;
        (void)dfdx;
        logger::log_assert(false, logger::ERROR, "Linear solver does not support nonlinear solver methods");

        return true;
    }

    virtual SolverClass get_class() const{
        return SolverClass::LINEAR;
    }

    virtual SolverType get_type() const{
        return SolverType::LINEAR;
    }
};

}

#endif
