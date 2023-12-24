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
#include "global_stiffness_matrix/petsc_sparse_symmetric.hpp"
#include "logger.hpp"


namespace global_stiffness_matrix{

PETScSparseSymmetric::~PETScSparseSymmetric(){
    MatDestroy(&this->K);
}

void PETScSparseSymmetric::generate(const Meshing* const mesh, const std::vector<long>& node_positions, const size_t matrix_width, const std::vector<double>& density, const double pc, const double psi){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id == 0){
        logger::quick_log("Generating stiffness matrix...");
    }

    if(first_time){
        this->preallocate(mesh, node_positions, matrix_width, density, pc, psi, mpi_id);
    } else {
        this->zero();
    }

    this->assemble_matrix(mesh, node_positions, density, pc, psi, mpi_id);

    this->first_time = false;
    if(mpi_id == 0){
        logger::quick_log("Done.");
    }
}

void PETScSparseSymmetricCPU::preallocate(const Meshing * const mesh, const std::vector<long>& node_positions, const size_t matrix_width, const std::vector<double>& density, const double pc, const double psi, const size_t mpi_id){
    MatCreate(PETSC_COMM_WORLD, &this->K);
    MatSetType(this->K, MATAIJ);

    long M = matrix_width;
    MPI_Bcast(&M, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    MatSetSizes(this->K, PETSC_DECIDE, PETSC_DECIDE, M, M);
    MatSetOption(this->K, MAT_SPD, PETSC_TRUE);
    MatSetOption(this->K, MAT_SPD_ETERNAL, PETSC_TRUE);
    MatSetOption(this->K, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    //MatSetOption(this->K, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);

    Mat tmp;
    MatCreate(PETSC_COMM_WORLD, &tmp);
    MatSetType(tmp, MATPREALLOCATOR);
    MatSetSizes(tmp, PETSC_DECIDE, PETSC_DECIDE, M, M);
    MatSetOption(tmp, MAT_SPD, PETSC_TRUE);
    MatSetOption(tmp, MAT_SPD_ETERNAL, PETSC_TRUE);
    MatSetOption(tmp, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    //MatSetOption(tmp, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
    MatSetUp(tmp);

    if(mpi_id == 0){
        std::swap(tmp, K);
        this->generate_base(mesh, node_positions, density, pc, psi);
        std::swap(tmp, K);
    }

    MatAssemblyBegin(tmp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(tmp, MAT_FINAL_ASSEMBLY);

    MatPreallocatorPreallocate(tmp, PETSC_TRUE, this->K);

    MatSetUp(this->K);

    MatDestroy(&tmp);
    //MatSetSizes(this->K, mesh->load_vector.size(), mesh->load_vector.size(), mesh->load_vector.size(), mesh->load_vector.size());
    //MatSetUp(this->K);
}

void PETScSparseSymmetricCPU::assemble_matrix(const Meshing * const mesh, const std::vector<long>& node_positions, const std::vector<double>& density, const double pc, const double psi, const size_t mpi_id){
    if(mpi_id == 0){
        this->generate_base(mesh, node_positions, density, pc, psi);
    }

    MatAssemblyBegin(this->K, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(this->K, MAT_FINAL_ASSEMBLY);
}
void PETScSparseSymmetricCUDA::preallocate(const Meshing * const mesh, const std::vector<long>& node_positions, const size_t matrix_width, const std::vector<double>& density, const double pc, const double psi, const size_t mpi_id){
    MatCreate(PETSC_COMM_WORLD, &this->K);
    MatSetType(this->K, MATAIJCUSPARSE);

    long M = matrix_width;
    MPI_Bcast(&M, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    MatSetSizes(this->K, PETSC_DECIDE, PETSC_DECIDE, M, M);
    MatSetOption(this->K, MAT_SPD, PETSC_TRUE);
    MatSetOption(this->K, MAT_SPD_ETERNAL, PETSC_TRUE);
    MatSetOption(this->K, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);

    if(mpi_id == 0){
        this->generate_base(mesh, node_positions, density, pc, psi);
    }

    this->K_coo.generate_coo(M);

    MatSetPreallocationCOO(this->K, this->K_coo.nvals, this->K_coo.rows.data(), this->K_coo.cols.data());

    MatSetUp(this->K);
}
void PETScSparseSymmetricCUDA::assemble_matrix(const Meshing * const mesh, const std::vector<long>& node_positions, const std::vector<double>& density, const double pc, const double psi, const size_t mpi_id){

    if(mpi_id == 0){
        if(!this->first_time){
            this->generate_base(mesh, node_positions, density, pc, psi);
        }
        MatSetValuesCOO(this->K, this->K_coo.vals.data(), INSERT_VALUES);
    }
}

}
