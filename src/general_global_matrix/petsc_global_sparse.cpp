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

#include "general_global_matrix/petsc_global_sparse.hpp"

namespace general_global_matrix{

PETScGlobalSparse::~PETScGlobalSparse(){
    MatDestroy(&this->M);
}

void PETScGlobalSparse::initialize(size_t L){
    this->L = L;
    MatCreate(PETSC_COMM_WORLD, &this->M);
    MatSetType(this->M, MATMPIAIJ);

    MPI_Bcast(&L, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    MatSetSizes(this->M, PETSC_DECIDE, PETSC_DECIDE, L, L);
    MatSetOption(this->M, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE);
}

void PETScGlobalSparse::begin_preallocation(){
    MatCreate(PETSC_COMM_WORLD, &this->tmp);
    MatSetType(tmp, MATPREALLOCATOR);
    MatSetSizes(tmp, PETSC_DECIDE, PETSC_DECIDE, L, L);
    MatSetOption(tmp, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE);
    MatSetUp(tmp);

    std::swap(tmp, M);
}

void PETScGlobalSparse::end_preallocation(){
    std::swap(tmp, M);

    MatAssemblyBegin(tmp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(tmp, MAT_FINAL_ASSEMBLY);

    MatPreallocatorPreallocate(tmp, PETSC_TRUE, this->M);

    MatSetUp(this->M);

    MatDestroy(&tmp);
}

}
