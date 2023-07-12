/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
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

#include <algorithm>
#include <cblas.h>
#include <mpich-x86_64/mpi.h>
#include "finite_element/gradient_descent.hpp"
#include "logger.hpp"
#include "optimization/MMASolver.hpp"

namespace finite_element{

GradientDescent::GradientDescent(const double eps, Solver solver):
    eps(eps), displacement(1), solver(solver){}

std::vector<double> GradientDescent::calculate_displacements(const Meshing* const mesh, std::vector<double> load, const std::vector<double>& density, double pc, double psi){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return std::vector<double>();
    }

    const size_t& W = this->gsm.get_W();
    const size_t& N = this->gsm.get_N();
    std::vector<double>& K = this->gsm.get_K();

    if(this->current_step == 0){
        this->gsm.generate(mesh, density, pc, psi);
        if(this->first_time){
            for(auto& u:this->displacement){
                u.resize(W,0);
            }
            this->first_time = false;
        }
    }
    logger::quick_log("Done.");
    logger::quick_log("Calculating displacements...");
    logger::quick_log("W: ",W," N: ", N);

    auto& u = this->displacement[this->current_step];
    auto d = load;
    size_t it = 0;
    auto comp_abs = [](const double a, const double b){
        return std::abs(a) < std::abs(b);
    };

    if(this->solver == Solver::STANDARD){
        double alpha = 0;
        double top, bot;
        auto d2 = d;
        while(std::abs(*std::max_element(d.begin(), d.end(),  comp_abs)) > this->eps){
            cblas_dcopy(W, load.data(), 1, d.data(), 1);
            cblas_dsbmv(CblasColMajor, CblasLower, W, N-1, -1.0, K.data(), N, u.data(), 1, 1.0, d.data(), 1);
            
            top = cblas_ddot(W, d.data(), 1, d.data(), 1);
            cblas_dsbmv(CblasColMajor, CblasLower, W, N-1, 1.0, K.data(), N, d.data(), 1, 0.0, d2.data(), 1);
            bot = cblas_ddot(W, d.data(), 1, d2.data(), 1);
            alpha = top/bot;

            cblas_daxpy(W, alpha, d.data(), 1, u.data(), 1);
            ++it;
        }
    } else if(this->solver == Solver::MMA){
        optimization::MMASolver mma(W, 0, 0, 0, 0); //1e5
        mma.SetAsymptotes(1e7, 0.08, 1.05);
        std::vector<double> xmin(W, -1000);
        std::vector<double> xmax(W,  1000);

        while(std::abs(*std::max_element(d.begin(), d.end(),  comp_abs)) > this->eps){
            cblas_dcopy(W, load.data(), 1, d.data(), 1);
            cblas_dsbmv(CblasColMajor, CblasLower, W, N-1, 1.0, K.data(), N, u.data(), 1, -1.0, d.data(), 1);

            mma.Update(u.data(), d.data(), nullptr, nullptr, xmin.data(), xmax.data());
            ++it;
        }
    } else if(this->solver == Solver::LAGRANGE_MMA){

        const double xmin = -100;
        const double xmax =  100;

        const double delta_max = 1e20;
        const double delta_min = 1e-20;
        const double delta = 1e8;
        std::vector<double> delta_e(W, delta);

        auto uold1 = u;
        auto uold2 = u;

        const double INC = 2;
        const double DEC = 0.08;

        while(std::abs(*std::max_element(d.begin(), d.end(),  comp_abs)) > this->eps){
            cblas_dcopy(W, load.data(), 1, d.data(), 1);
            cblas_dsbmv(CblasColMajor, CblasLower, W, N-1, 1.0, K.data(), N, u.data(), 1, -1.0, d.data(), 1);

            #pragma omp parallel for
            for(size_t i = 0; i < W; ++i){
                double d_e = (u[i] - uold1[i])*(uold1[i] - uold2[i]);
                if(d_e < 0){
                    delta_e[i] = std::max(DEC*delta_e[i], delta_min);
                } else if(d_e > 0){
                    delta_e[i] = std::min(INC*delta_e[i], delta_max);
                }
                double x_inf = u[i] - delta_e[i];
                double x_sup = u[i] + delta_e[i];
                double L = u[i] - 2*delta_max;
                double U = u[i] + 2*delta_max;
                double a = std::max({x_inf, xmin, 0.9*L+0.1*u[i]});
                double b = std::min({x_sup, xmax, 0.9*U+0.1*u[i]});

                double p = std::pow(U-u[i], 2)*(std::max(d[i], 0.0) + 0.001*std::abs(d[i]) + 0.5*1e-6*(U-L));
                double q = std::pow(u[i]-L, 2)*(std::max(-d[i], 0.0) + 0.001*std::abs(d[i]) + 0.5*1e-6*(U-L));
                double sqrtp = std::sqrt(p);
                double sqrtq = std::sqrt(q);
                uold2[i] = uold1[i];
                uold1[i] = u[i];
                u[i] = std::max(a, std::min(b, (L*sqrtp+U*sqrtq)/(sqrtp+sqrtq)));
            }
            logger::quick_log(std::abs(*std::max_element(d.begin(), d.end(),  comp_abs)));
            
            ++it;
        }
    }

    logger::quick_log("Required iterations: ", it);
    logger::quick_log("");
    logger::quick_log("Done.");

    this->current_step = (this->current_step + 1) % this->steps;
   
    return u;
}

}
