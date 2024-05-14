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

#include "nonlinear_solver/steepest_descent.hpp"

namespace nonlinear_solver{

SteepestDescent::SteepestDescent(double INIT, double INC, double DEC, double xtol_abs):
    NonlinearSolver(xtol_abs, 0), INIT(INIT), INC(INC), DEC(DEC), STEP(INIT){

}

void SteepestDescent::setup(size_t N){
    this->N = N;
    this->xold.resize(N);
    this->it = 0;
    std::fill(this->xold.begin(), this->xold.end(), 0);
}

bool SteepestDescent::update(double* x, double f, const double* dfdx){
    if(this->it >= 2){
        if((fold2 - fold1)*(fold1 - f) <= 0){
            STEP *= DEC;
        } else {
            STEP *= INC;
            STEP = std::min(1.0, STEP);
        }
    }
    fold2 = fold1;
    fold1 = f;

    double maxl = 0;
    for(size_t i = 0; i < N; ++i){
        maxl = std::max(std::abs(dfdx[i]), maxl);
    }

    double ch = 0.0;
    for(size_t i = 0; i < N; ++i){
        xold[i] = x[i];
        x[i] -= STEP*dfdx[i]/maxl;
        ch = std::max(std::abs(x[i] - xold[i]), ch);
    }

    ++this->it;

    return ch < this->xtol_abs;
}

}
