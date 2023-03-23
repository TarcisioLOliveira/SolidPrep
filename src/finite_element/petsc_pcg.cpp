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

PETScPCG::PETScPCG():gsm(), u(1){
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
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return std::vector<double>();
    }

    if(this->indices.size() == 0){
        indices.resize(load.size());
        std::iota(indices.begin(), indices.end(), 0);
    }

    if(this->current_step == 0){
        this->gsm.generate(mesh, density, pc, psi);
    }

    if(this->first_time){
        VecCreate(PETSC_COMM_WORLD, &this->u[0]);
        VecSetType(this->u[0], VECSTANDARD);
        VecSetSizes(this->u[0], PETSC_DECIDE, load.size());
        VecSetUp(this->u[0]);

        VecCreate(PETSC_COMM_WORLD, &this->f);
        VecSetType(this->f, VECSTANDARD);
        VecSetSizes(this->f, PETSC_DECIDE, load.size());
        VecSetUp(this->f);

        KSPCreate(PETSC_COMM_WORLD, &this->ksp);
        KSPSetType(this->ksp, KSPCG);
        KSPCGSetType(this->ksp, KSP_CG_HERMITIAN);
        KSPSetInitialGuessNonzero(this->ksp, PETSC_TRUE);

        KSPGetPC(this->ksp, &this->pc);
        PCSetType(this->pc, PCJACOBI);

        this->first_time = false;
    }

    if(this->u[current_step] == 0){
        VecCreate(PETSC_COMM_WORLD, &this->u[current_step]);
        VecSetType(this->u[current_step], VECSTANDARD);
        VecSetSizes(this->u[current_step], PETSC_DECIDE, load.size());
        VecSetUp(this->u[current_step]);
    }

    auto K = this->gsm.get_K();
    VecSetValues(this->f, load.size(), indices.data(), load.data(), INSERT_VALUES);

    VecAssemblyBegin(this->u[current_step]);
    VecAssemblyEnd(this->u[current_step]);

    VecAssemblyBegin(this->f);
    VecAssemblyEnd(this->f);

    KSPSetOperators(this->ksp, K, K);

    KSPSetUp(this->ksp);
    PCSetUp(this->pc);

    KSPSolve(this->ksp, this->f, this->u[current_step]);

    const double* load_data;

    VecGetArrayRead(u[current_step], &load_data);

    std::copy(load_data, load_data+load.size(), load.begin());

    VecRestoreArrayRead(u[current_step], &load_data);

    this->current_step = (this->current_step + 1) % this->steps;
   
    return load; 
}

}
