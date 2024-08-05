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

#include "nonlinear_solver/newton.hpp"

namespace nonlinear_solver{

Newton::Newton(double MULT, double rtol_abs):
    NonlinearSolver(0, rtol_abs), MULT(MULT){

}

void Newton::setup(size_t N){
    this->N = N;
    this->it = 0;
}

bool Newton::update(double* x, double f, const double* dfdx){
    f = std::min(f, MULT);
    for(size_t i = 0; i < N; ++i){
        x[i] += f*dfdx[i];
    }

    ++this->it;

    return false;
}

}
