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

#include "nonlinear_solver/mma.hpp"
#include "logger.hpp"

namespace nonlinear_solver{

MMA::MMA(double INIT, double INC, double DEC, double delta_min, double delta_max, double xmin, double xmax, double xtol_abs):
    NonlinearSolver(xtol_abs, 0), INIT(INIT), INC(INC), DEC(DEC),
        delta_min(delta_min), delta_max(delta_max), xmin(xmin), xmax(xmax){

}

void MMA::setup(size_t N){
    this->N = N;
    this->xold1.resize(N);
    this->xold2.resize(N);
    this->delta_e.resize(N);
    this->it = 0;
    std::fill(this->xold1.begin(), this->xold1.end(), 0);
    std::fill(this->xold2.begin(), this->xold2.end(), 0);
    std::fill(this->delta_e.begin(), this->delta_e.end(), INIT);
}

bool MMA::update(double* x, double f, const double* dfdx){
    (void)f;
    double ch = 0;
    for(size_t i = 0; i < this->N; ++i){
        const double S = dfdx[i];
        double d_e = (x[i] - xold1[i])*(xold1[i] - xold2[i]);
        if(d_e < 0){
            delta_e[i] = std::max(DEC*delta_e[i], delta_min);
            //delta_e[i] *= DEC;
        } else if(d_e > 0){
            delta_e[i] = std::min(INC*delta_e[i], delta_max);
            //delta_e[i] *= INC;
        }
        double x_inf = x[i] - delta_e[i];
        double x_sup = x[i] + delta_e[i];
        double L = x[i] - 2*delta_max;
        double U = x[i] + 2*delta_max;
        double a = std::max({x_inf, xmin, 0.9*L+0.1*x[i]});
        double b = std::min({x_sup, xmax, 0.9*U+0.1*x[i]});

        double p = std::pow(U-x[i], 2)*(std::max(S, 0.0) + 0.001*std::abs(S) + 0.5*1e-6*(U-L));
        double q = std::pow(x[i]-L, 2)*(std::max(-S, 0.0) + 0.001*std::abs(S) + 0.5*1e-6*(U-L));
        double sqrtp = std::sqrt(p);
        double sqrtq = std::sqrt(q);
        xold2[i] = xold1[i];
        xold1[i] = x[i];
        x[i] = std::max(a, std::min(b, (L*sqrtp+U*sqrtq)/(sqrtp+sqrtq)));
        ch = std::max(std::abs(x[i] - xold1[i]), ch);
    }
    logger::quick_log("dx:", ch);

    ++this->it;

    return ch < this->xtol_abs;
}

}
