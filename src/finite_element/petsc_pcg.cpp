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


#include <mpich-x86_64/mpi.h>
#include "finite_element/petsc_pcg.hpp"
#include "logger.hpp"

namespace finite_element{

PETScPCG::PETScPCG(PETScBackend backend):
    gsm(nullptr), u(1){
    switch(backend){
        case PETScBackend::CPU:
            this->vec_type = VECSTANDARD;
            this->gsm = std::make_unique<global_stiffness_matrix::PETScSparseSymmetricCPU>();
            break;
        case PETScBackend::CUDA:
            this->vec_type = VECCUDA;
            this->gsm = std::make_unique<global_stiffness_matrix::PETScSparseSymmetricCUDA>();
            break;
    }
}

PETScPCG::~PETScPCG(){
    VecDestroy(&this->f);
    for(auto& v:this->u){
        VecDestroy(&v);
    }
    KSPDestroy(&this->ksp);
}

std::vector<double> PETScPCG::calculate_displacements(const Meshing* const mesh, std::vector<double> load, const std::vector<double>& density, double pc, double psi){
    int mpi_id = 0;
    int mpi_size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    if(this->current_step == 0){
        this->gsm->generate(mesh, density, pc, psi);
    }

    auto K = this->gsm->get_K();

    long M = load.size();
    long n = 0, m = 0;
    MatGetLocalSize(K, &n, &m);
    MPI_Bcast(&M, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    if(this->first_time){
        // DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, M, 1, 1, NULL, &this->dm);
        // DMSetVecType(this->dm, VECSTANDARD);
        // DMSetUp(this->dm);

        VecCreateMPI(PETSC_COMM_WORLD, m, M, &this->u[0]);
        VecSetType(this->u[0], this->vec_type.c_str());
        //VecSetSizes(this->u[0], PETSC_DECIDE, M);
        //DMCreateGlobalVector(this->dm, &this->u[0]);
        VecSetUp(this->u[0]);

        VecCreateMPI(PETSC_COMM_WORLD, m, M, &this->f);
        VecSetType(this->f, this->vec_type.c_str());
        //VecSetSizes(this->f, PETSC_DECIDE, M);
        //DMCreateGlobalVector(this->dm, &this->f);
        VecSetUp(this->f);

        KSPCreate(PETSC_COMM_WORLD, &this->ksp);
        KSPSetType(this->ksp, KSPCG);
        KSPCGSetType(this->ksp, KSP_CG_HERMITIAN);
        KSPSetInitialGuessNonzero(this->ksp, PETSC_TRUE);

        KSPGetPC(this->ksp, &this->pc);
        PCSetType(this->pc, PCGAMG);

        this->first_time = false;
    }

    if(this->u[current_step] == 0){
        VecCreateMPI(PETSC_COMM_WORLD, m, M, &this->u[current_step]);
        VecSetType(this->u[current_step], this->vec_type.c_str());
        //VecSetSizes(this->u[current_step], PETSC_DECIDE, load.size());
        //DMCreateGlobalVector(this->dm, &this->u[current_step]);
        VecSetUp(this->u[current_step]);
    }

    if(mpi_id != 0){
        load.resize(M);
    }
    MPI_Bcast(load.data(), load.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    long begin = 0, end = 0;
    VecGetOwnershipRange(this->f, &begin, &end);

    double* f_data = nullptr;
    VecGetArray(this->f, &f_data);
    for(auto d = 0; d < end - begin; ++d){
        *(f_data+d) = load[begin+d];
    }
    VecRestoreArray(this->f, &f_data);

    // if(mpi_id == 0){
    //     VecSetValues(this->f, load.size(), indices.data(), load.data(), INSERT_VALUES);
    // }

    VecAssemblyBegin(this->u[current_step]);
    VecAssemblyEnd(this->u[current_step]);

    VecAssemblyBegin(this->f);
    VecAssemblyEnd(this->f);

    if(current_step == 0){
        KSPSetOperators(this->ksp, K, K);

        KSPSetUp(this->ksp);
        PCSetUp(this->pc);
    }

    KSPSolve(this->ksp, this->f, this->u[current_step]);

    const double* load_data;

    VecGetArrayRead(u[current_step], &load_data);

    if(mpi_size > 1){
        if(mpi_id == 0){
            std::fill(load.begin(), load.end(), 0);
            std::copy(load_data, load_data + m, load.begin());
            std::vector<double> load_data2(2*m,0);
            long l = 0;
            long step = m;
            for(int i = 1; i < mpi_size; ++i){
                MPI_Status mpi_status;
                MPI_Recv(&l, 1, MPI_DOUBLE, i, 111, MPI_COMM_WORLD, &mpi_status);
                MPI_Recv(load_data2.data(), l, MPI_DOUBLE, i, 111, MPI_COMM_WORLD, &mpi_status);
                for(long j = 0; j < m; ++j){
                    load[j+step] += load_data2[j];
                }
                step += l;
            }
        } else {
            long l = end - begin;
            MPI_Send(&l, 1, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD);
            MPI_Send(load_data, m, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD);
        }
    } else {
        std::copy(load_data, load_data + M, load.begin());
    }
    MPI_Barrier(MPI_COMM_WORLD);

    VecRestoreArrayRead(this->u[current_step], &load_data);

    this->current_step = (this->current_step + 1) % this->steps;
   
    return load; 
}

}
