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

#include "general_solver/petsc_general_pcg.hpp"

namespace general_solver{

PETScGeneralPCG::~PETScGeneralPCG(){
    VecDestroy(&this->f);
    KSPDestroy(&this->ksp);
}


void PETScGeneralPCG::initialize_matrix(bool spd, size_t L){
    this->spd = spd;
    this->L = L;

    this->M.initialize(L);
    auto mat = M.get_matrix();

    long n = 0, m = 0;
    MatGetLocalSize(mat, &n, &m);
    MPI_Bcast(&L, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    VecCreateMPI(PETSC_COMM_WORLD, m, L, &this->f);
    VecSetType(this->f, VECCUDA);
    VecSetUp(this->f);

    VecCreateMPI(PETSC_COMM_WORLD, m, L, &this->u);
    VecSetType(this->u, VECCUDA);
    VecSetUp(this->u);

    KSPCreate(PETSC_COMM_WORLD, &this->ksp);
    if(spd){
        KSPSetType(this->ksp, KSPCG);
    } else {
        KSPSetType(this->ksp, KSPMINRES);
    }
    KSPSetInitialGuessNonzero(this->ksp, PETSC_FALSE);
    KSPSetNormType(this->ksp, KSP_NORM_UNPRECONDITIONED);
    
    KSPGetPC(this->ksp, &this->pc);
    PCSetType(this->pc, PCJACOBI);
    PCJacobiSetType(this->pc, PC_JACOBI_DIAGONAL);
    PCFactorSetUseInPlace(this->pc, PETSC_TRUE);
    PCFactorSetReuseOrdering(this->pc, PETSC_TRUE);

    KSPSetOperators(this->ksp, mat, mat);
}

void PETScGeneralPCG::solve(std::vector<double>& b){
    int mpi_id = 0;
    int mpi_size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    long begin = 0, end = 0;
    VecGetOwnershipRange(this->f, &begin, &end);
    double* f_data = nullptr;
    VecGetArray(this->f, &f_data);
    for(auto d = 0; d < end - begin; ++d){
        *(f_data+d) = b[begin+d];
    }
    VecRestoreArray(this->f, &f_data);

    VecAssemblyBegin(this->u);
    VecAssemblyEnd(this->u);

    VecAssemblyBegin(this->f);
    VecAssemblyEnd(this->f);

    KSPSetUp(this->ksp);
    PCSetUp(this->pc);

    KSPSolve(this->ksp, this->f, this->u);
    KSPConvergedReason r;
    KSPGetConvergedReason(this->ksp, &r);
    logger::quick_log("Converged?", r);
    {
        double tmp = 0;
        PetscInt its = 0;
        KSPGetResidualNorm(this->ksp, &tmp);
        KSPGetTotalIterations(this->ksp, &its);
        logger::quick_log("Residual norm", tmp);
        logger::quick_log("Total iterations", its);
    }

    auto mat = M.get_matrix();
    long n = 0, m = 0;
    MatGetLocalSize(mat, &n, &m);
    VecGetOwnershipRange(this->u, &begin, &end);

    const double* load_data;

    VecGetArrayRead(u, &load_data);

    if(mpi_size > 1){
        if(mpi_id == 0){
            std::fill(b.begin(), b.end(), 0);
            std::copy(load_data, load_data + m, b.begin());
            std::vector<double> load_data2(2*m,0);
            long l = 0;
            long step = m;
            for(int i = 1; i < mpi_size; ++i){
                MPI_Status mpi_status;
                MPI_Recv(&l, 1, MPI_DOUBLE, i, 111, MPI_COMM_WORLD, &mpi_status);
                MPI_Recv(load_data2.data(), l, MPI_DOUBLE, i, 111, MPI_COMM_WORLD, &mpi_status);
                for(long j = 0; j < m; ++j){
                    b[j+step] += load_data2[j];
                }
                step += l;
            }
        } else {
            long l = end - begin;
            MPI_Send(&l, 1, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD);
            MPI_Send(load_data, m, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD);
        }
    } else {
        std::copy(load_data, load_data + this->L, b.begin());
    }
    MPI_Barrier(MPI_COMM_WORLD);

    VecRestoreArrayRead(this->u, &load_data);
}

}
