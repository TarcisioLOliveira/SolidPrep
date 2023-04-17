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


void PETScGeneralPCG::initialize(general_global_matrix::PETScGlobalSparse* M, size_t L){
    this->M = M;
    this->L = L;

    auto mat = M->get_matrix();

    long n = 0, m = 0;
    MatGetLocalSize(mat, &n, &m);
    MPI_Bcast(&L, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    VecCreateMPI(PETSC_COMM_WORLD, m, L, &this->f);
    VecSetType(this->f, VECSTANDARD);
    VecSetUp(this->f);

    VecCreateMPI(PETSC_COMM_WORLD, m, L, &this->u);
    VecSetType(this->u, VECSTANDARD);
    VecSetUp(this->u);

    KSPCreate(PETSC_COMM_WORLD, &this->ksp);
    KSPSetType(this->ksp, KSPCGNE);
    KSPSetInitialGuessNonzero(this->ksp, PETSC_FALSE);

    KSPGetPC(this->ksp, &this->pc);
    PCSetType(this->pc, PCJACOBI);

    KSPSetOperators(this->ksp, mat, mat);

    KSPSetUp(this->ksp);
    PCSetUp(this->pc);
}

void PETScGeneralPCG::set_rhs(std::vector<double> b){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        b.resize(L);
    }
    MPI_Bcast(b.data(), b.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

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

}

void PETScGeneralPCG::solve(std::vector<double>& x){
    int mpi_id = 0;
    int mpi_size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    long begin = 0, end = 0;
    VecGetOwnershipRange(this->u, &begin, &end);
    double* u_data = nullptr;
    VecGetArray(this->u, &u_data);
    for(auto d = 0; d < end - begin; ++d){
        *(u_data+d) = x[begin+d];
    }
    VecRestoreArray(this->u, &u_data);

    KSPSolve(this->ksp, this->f, this->u);

    auto mat = M->get_matrix();
    long n = 0, m = 0;
    MatGetLocalSize(mat, &n, &m);
    VecGetOwnershipRange(this->f, &begin, &end);

    const double* load_data;

    VecGetArrayRead(u, &load_data);

    if(mpi_size > 1){
        if(mpi_id == 0){
            std::fill(x.begin(), x.end(), 0);
            std::copy(load_data, load_data + m, x.begin());
            std::vector<double> load_data2(2*m,0);
            long l = 0;
            long step = m;
            for(int i = 1; i < mpi_size; ++i){
                MPI_Status mpi_status;
                MPI_Recv(&l, 1, MPI_DOUBLE, i, 111, MPI_COMM_WORLD, &mpi_status);
                MPI_Recv(load_data2.data(), l, MPI_DOUBLE, i, 111, MPI_COMM_WORLD, &mpi_status);
                for(long j = 0; j < m; ++j){
                    x[j+step] += load_data2[j];
                }
                step += l;
            }
        } else {
            long l = end - begin;
            MPI_Send(&l, 1, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD);
            MPI_Send(load_data, m, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD);
        }
    } else {
        std::copy(load_data, load_data + this->L, x.begin());
    }
    MPI_Barrier(MPI_COMM_WORLD);

    VecRestoreArrayRead(this->u, &load_data);
}

}
