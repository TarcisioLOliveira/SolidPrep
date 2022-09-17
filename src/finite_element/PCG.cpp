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
#include "finite_element/PCG.hpp"
#include "logger.hpp"

namespace finite_element{

PCG::PCG(const double eps, const Preconditioner precond):
    eps(eps), precond(precond), displacement(1), P(){}

std::vector<double> PCG::calculate_displacements(const Meshing* const mesh, std::vector<double> load, const std::vector<double>& density, double pc){

    if(this->W == 0 || this->N == 0){
        this->calculate_dimensions(mesh, load);
        for(auto& u:this->displacement){
            u.resize(W,0);
        }
        this->P.resize(W, 0);
    }

    if(this->current_step == 0){
        this->generate_K(mesh, density, pc);
        this->generate_P();
    }
    logger::quick_log("Done.");
    logger::quick_log("Calculating displacements...");
    logger::quick_log("W: ",W," N: ", N);

    auto& u = this->displacement[this->current_step];


    auto r = load;
    auto z = r;
    cblas_dsbmv(CblasColMajor, CblasLower, W, N-1, -1.0, this->K.data(), N, u.data(), 1, 1.0, r.data(), 1);
    if(this->precond == Preconditioner::JACOBI){
        cblas_dsbmv(CblasColMajor, CblasLower, W, 0, 1.0, this->P.data(), 1, r.data(), 1, 0.0, z.data(), 1);
    } else if(this->precond == Preconditioner::SSOR){

    }
    double rho1 = cblas_ddot(W, z.data(), 1, r.data(), 1);
    double rho2 = 0;
    auto p = z;
    auto q = p;
    cblas_dsbmv(CblasColMajor, CblasLower, W, N-1, 1.0, this->K.data(), N, p.data(), 1, 0.0, q.data(), 1);
    double alpha = rho1/cblas_ddot(W, p.data(), 1, q.data(), 1);
    cblas_daxpy(W, alpha, p.data(), 1, u.data(), 1);
    cblas_daxpy(W, -alpha, q.data(), 1, r.data(), 1);

    double beta = 0;
    size_t it = 1;
    auto comp_abs = [](const double a, const double b){
        return std::abs(a) < std::abs(b);
    };
    while(std::abs(*std::max_element(r.begin(), r.end(),  comp_abs)) > this->eps){
        if(this->precond == Preconditioner::JACOBI){
            cblas_dsbmv(CblasColMajor, CblasLower, W, 0, 1.0, this->P.data(), 1, r.data(), 1, 0.0, z.data(), 1);
        } else if(this->precond == Preconditioner::SSOR){

        }
        rho2 = rho1;
        rho1 = cblas_ddot(W, z.data(), 1, r.data(), 1);
        beta = rho1/rho2;
        cblas_daxpy(W, beta, p.data(), 1, z.data(), 1);
        std::swap(z, p);

        cblas_dsbmv(CblasColMajor, CblasLower, W, N-1, 1.0, this->K.data(), N, p.data(), 1, 0.0, q.data(), 1);
        alpha = rho1/cblas_ddot(W, p.data(), 1, q.data(), 1);
        cblas_daxpy(W, alpha, p.data(), 1, u.data(), 1);
        cblas_daxpy(W, -alpha, q.data(), 1, r.data(), 1);

        ++it;
    }

    logger::quick_log("Required iterations: ", it);
    logger::quick_log("");
    logger::quick_log("Done.");

    this->current_step = (this->current_step + 1) % this->steps;
   
    return load; 
}


void PCG::generate_P(){
    if(this->precond == Preconditioner::JACOBI){
        #pragma omp parallel for
        for(size_t i = 0; i < W; ++i){
            this->P[i] = 1.0/this->K[i*N];
        }
    } else if(this->precond == Preconditioner::SSOR){

    }
}

}
