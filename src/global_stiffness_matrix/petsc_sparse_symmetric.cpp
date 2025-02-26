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

void PETScSparseSymmetric::generate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
    this->u_size = u_size;
    this->l_num = l_num;

    if(mpi_id == 0){
        logger::quick_log("Generating stiffness matrix...");
    }

    if(first_time){
        this->preallocate(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, lambda, type, mpi_id);
    } else {
        this->zero();
    }

    this->assemble_matrix(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, lambda, type, mpi_id);

    this->first_time = false;
    if(mpi_id == 0){
        logger::quick_log("Done.");
    }
}

void PETScSparseSymmetricCPU::preallocate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type, const size_t mpi_id){
    MatCreate(PETSC_COMM_WORLD, &this->K);
    MatSetType(this->K, MATAIJ);

    long M = u_size;
    if(type != FiniteElement::ContactType::RIGID){
        M += l_num;
    }
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
        this->generate_base(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, lambda, type);
        std::swap(tmp, K);
    }

    MatAssemblyBegin(tmp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(tmp, MAT_FINAL_ASSEMBLY);

    MatPreallocatorPreallocate(tmp, PETSC_TRUE, this->K);

    MatSetUp(this->K);

    MatDestroy(&tmp);
}

void PETScSparseSymmetricCPU::assemble_matrix(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type, const size_t mpi_id){
    if(mpi_id == 0){
        this->generate_base(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, lambda, type);
    }

    MatAssemblyBegin(this->K, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(this->K, MAT_FINAL_ASSEMBLY);
}

void PETScSparseSymmetricCPU::dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const{
    int mpi_id = 0;
    int mpi_size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    long M = v.size();
    long n = 0, m = 0;
    MatGetLocalSize(K, &n, &m);
    MPI_Bcast(&M, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    Vec v1, v2;
    VecCreateMPI(PETSC_COMM_WORLD, m, M, &v1);
    VecCreateMPI(PETSC_COMM_WORLD, m, M, &v2);
    VecSetType(v1, VECSTANDARD);
    VecSetType(v2, VECSTANDARD);
    VecSetUp(v1);
    VecSetUp(v2);

    long begin = 0, end = 0;
    VecGetOwnershipRange(v1, &begin, &end);

    double* v_data = nullptr;
    VecGetArray(v1, &v_data);
    std::copy(v.begin() + begin, v.begin() + end, v_data);
    VecRestoreArray(v1, &v_data);

    MatMult(this->K, v1, v2);

    const double* v_out_data;
    VecGetArrayRead(v2, &v_out_data);

    if(mpi_size > 1){
        if(mpi_id == 0){
            std::copy(v_out_data, v_out_data + m, v_out.begin());
            std::vector<double> v_out_data2(2*m,0);
            long l = 0;
            long step = m;
            for(int i = 1; i < mpi_size; ++i){
                MPI_Status mpi_status;
                MPI_Recv(&l, 1, MPI_DOUBLE, i, 111, MPI_COMM_WORLD, &mpi_status);
                MPI_Recv(v_out_data2.data(), l, MPI_DOUBLE, i, 111, MPI_COMM_WORLD, &mpi_status);
                for(long j = 0; j < m; ++j){
                    v_out[j+step] += v_out_data2[j];
                }
                step += l;
            }
        } else {
            long l = end - begin;
            MPI_Send(&l, 1, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD);
            MPI_Send(v_out_data, m, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD);
        }
    } else {
        std::copy(v_out_data, v_out_data + M, v_out.begin());
    }
    MPI_Barrier(MPI_COMM_WORLD);
    VecRestoreArrayRead(v2, &v_out_data);
}

void PETScSparseSymmetricCUDA::preallocate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type, const size_t mpi_id){
    MatCreate(PETSC_COMM_WORLD, &this->K);
    MatSetType(this->K, MATAIJCUSPARSE);

    long M = u_size;
    if(type != FiniteElement::ContactType::RIGID){
        M += l_num;
    }
    MPI_Bcast(&M, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    MatSetSizes(this->K, PETSC_DECIDE, PETSC_DECIDE, M, M);
    MatSetOption(this->K, MAT_SPD, PETSC_TRUE);
    MatSetOption(this->K, MAT_SPD_ETERNAL, PETSC_TRUE);
    MatSetOption(this->K, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);

    if(mpi_id == 0){
        this->generate_base(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, lambda, type);
    }

    this->K_coo.generate_coo(M);

    MatSetPreallocationCOO(this->K, this->K_coo.nvals, this->K_coo.rows.data(), this->K_coo.cols.data());

    MatSetUp(this->K);
}
void PETScSparseSymmetricCUDA::assemble_matrix(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type, const size_t mpi_id){

    if(mpi_id == 0){
        if(!this->first_time){
            this->generate_base(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, lambda, type);
        }
        if(type >= FiniteElement::ContactType::FRICTIONLESS_DISPL_SIMPLE){
            this->K_coo.backup_matrix();
        }
        MatSetValuesCOO(this->K, this->K_coo.vals.data(), INSERT_VALUES);
    }
}

void PETScSparseSymmetricCUDA::reset_hessian(){
    this->K_coo.restore_matrix();
}

}
