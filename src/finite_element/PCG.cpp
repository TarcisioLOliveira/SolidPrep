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
#include <lapacke.h>
#include <mpich-x86_64/mpi.h>
#include "finite_element/PCG.hpp"
#include "logger.hpp"

namespace finite_element{

PCG::PCG(const double eps, const Preconditioner precond):
    eps(eps), precond(precond), P(){}

void PCG::generate_matrix(const Meshing* const mesh, const size_t L, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache){
    this->gsm.generate(mesh, node_positions, L, topopt, D_cache);
    this->setup = false;
}

void PCG::calculate_displacements(std::vector<double>& load){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return;
    }

    const size_t& W = this->gsm.get_W();
    const size_t& N = this->gsm.get_N();
    std::vector<double>& K = this->gsm.get_K();

    std::vector<double> u(load.size(), 0);
    auto r = load;
    auto z = r;
    auto z2 = z;
    if(!this->setup){
        if(this->first_time){
            this->P.resize(W, 0);
            this->first_time = false;
        }
        this->generate_P();
        if(this->precond == Preconditioner::JACOBI){
            cblas_dsbmv(CblasColMajor, CblasLower, W, 0, 1.0, this->P.data(), 1, r.data(), 1, 0.0, z.data(), 1);
        } else if(this->precond == Preconditioner::SSOR){
            LAPACKE_dtbtrs_work(LAPACK_COL_MAJOR, 'L', 'N', 'N', W, N-1, 1, K.data(), N, z.data(), W);
            cblas_dsbmv(CblasColMajor, CblasLower, W, 0, 1.0, this->P.data(), 1, z2.data(), 1, 0.0, z.data(), 1);
            LAPACKE_dtbtrs_work(LAPACK_COL_MAJOR, 'L', 'T', 'N', W, N-1, 1, K.data(), N, z.data(), W);
        }
    }
    logger::quick_log("Done.");
    logger::quick_log("Calculating displacements...");
    logger::quick_log("W: ",W," N: ", N);

    cblas_dsbmv(CblasColMajor, CblasLower, W, N-1, -1.0, K.data(), N, u.data(), 1, 1.0, r.data(), 1);
    double rho1 = cblas_ddot(W, z.data(), 1, r.data(), 1);
    double rho2 = 0;
    auto p = z;
    auto q = p;
    cblas_dsbmv(CblasColMajor, CblasLower, W, N-1, 1.0, K.data(), N, p.data(), 1, 0.0, q.data(), 1);
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
            cblas_dcopy(W, r.data(), 1, z.data(), 1);
            LAPACKE_dtbtrs_work(LAPACK_COL_MAJOR, 'L', 'N', 'N', W, N-1, 1, K.data(), N, z.data(), W);
            cblas_dcopy(W, z.data(), 1, z2.data(), 1);
            cblas_dsbmv(CblasColMajor, CblasLower, W, 0, 1.0, this->P.data(), 1, z2.data(), 1, 0.0, z.data(), 1);
            LAPACKE_dtbtrs_work(LAPACK_COL_MAJOR, 'L', 'T', 'N', W, N-1, 1, K.data(), N, z.data(), W);
        }
        rho2 = rho1;
        rho1 = cblas_ddot(W, z.data(), 1, r.data(), 1);
        beta = rho1/rho2;
        cblas_daxpy(W, beta, p.data(), 1, z.data(), 1);
        std::swap(z, p);

        cblas_dsbmv(CblasColMajor, CblasLower, W, N-1, 1.0, K.data(), N, p.data(), 1, 0.0, q.data(), 1);
        alpha = rho1/cblas_ddot(W, p.data(), 1, q.data(), 1);
        cblas_daxpy(W, alpha, p.data(), 1, u.data(), 1);
        cblas_daxpy(W, -alpha, q.data(), 1, r.data(), 1);

        ++it;
    }

    logger::quick_log("Required iterations: ", it);
    logger::quick_log("");
    logger::quick_log("Done.");
    
    cblas_dcopy(load.size(), u.data(), 1, load.data(), 1);
}


void PCG::generate_P(){
    const size_t& W = this->gsm.get_W();
    const size_t& N = this->gsm.get_N();
    std::vector<double>& K = this->gsm.get_K();
    if(this->precond == Preconditioner::JACOBI || this->precond == Preconditioner::SSOR){
        #pragma omp parallel for
        for(size_t i = 0; i < W; ++i){
            this->P[i] = 1.0/K[i*N];
        }
    }
}

}
