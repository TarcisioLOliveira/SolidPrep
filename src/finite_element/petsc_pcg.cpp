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
#include <vector>
#include "finite_element/petsc_pcg.hpp"
#include "math/matrix.hpp"
#include "logger.hpp"
#include "project_data.hpp"

namespace finite_element{

PETScPCG::PETScPCG(const projspec::DataMap& data):
    FiniteElement(data.proj->contact_data),
    gsm(nullptr)
{
    const std::string backend(data.get_string("backend"));
    if(backend == "cpu"){
        this->vec_type = VECSTANDARD;
        this->gsm = std::make_unique<global_stiffness_matrix::PETScSparseSymmetricCPU>(data.proj->contact_data);
    } else if(backend == "cuda"){
        this->vec_type = VECCUDA;
        this->gsm = std::make_unique<global_stiffness_matrix::PETScSparseSymmetricCUDA>(data.proj->contact_data);
    } else {
        logger::log_assert(false, logger::ERROR, "undefined PETSc backend: {}", backend);
    }
    this->set_global_matrix(this->gsm.get());
}

PETScPCG::~PETScPCG(){
    VecDestroy(&this->f);
    VecDestroy(&this->u);
    KSPDestroy(&this->ksp);
    PCDestroy(&this->pc);
}

void PETScPCG::generate_matrix_base(const Meshing* const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const ContactType type){
    this->gsm->generate(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, lambda, type);
    this->nonlinear_boundary_conditions = mesh->nonlinear_boundary_conditions;
    this->spd_matrix = (mesh->internal_loads2 == nullptr || mesh->internal_loads2->size() == 0);
    this->setup = false;
}

void PETScPCG::solve_init(size_t M){
    auto K = this->gsm->get_K();

    long n = 0, m = 0;
    MatGetLocalSize(K, &n, &m);
    MPI_Bcast(&M, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    VecCreateMPI(PETSC_COMM_WORLD, m, M, &this->u);
    VecSetType(this->u, this->vec_type.c_str());
    VecSetUp(this->u);

    VecCreateMPI(PETSC_COMM_WORLD, m, M, &this->f);
    VecSetType(this->f, this->vec_type.c_str());
    VecSetUp(this->f);

    KSPCreate(PETSC_COMM_WORLD, &this->ksp);
    KSPSetNormType(this->ksp, KSP_NORM_UNPRECONDITIONED);

    KSPGetPC(this->ksp, &this->pc);
    PCFactorSetUseInPlace(this->pc, PETSC_TRUE);
    KSPSetInitialGuessNonzero(this->ksp, PETSC_FALSE);
    PCFactorSetReuseOrdering(this->pc, PETSC_TRUE);

    PCSetType(this->pc, PCJACOBI);
    PCJacobiSetType(this->pc, PC_JACOBI_DIAGONAL);
    PCJacobiSetUseAbs(this->pc, PETSC_TRUE);
    PCJacobiSetFixDiagonal(this->pc, PETSC_TRUE);

    if(this->spd_matrix){
        KSPSetType(this->ksp, KSPCG);
        KSPCGSetType(this->ksp, KSP_CG_SYMMETRIC);
        KSPSetTolerances(this->ksp, 1e-5, 1e-50, 1e10, 1e5);
    } else {
        KSPSetTolerances(this->ksp, 1e-7, 1e-50, 1e10, 1e6);
        PCJacobiSetType(this->pc, PC_JACOBI_ROWMAX);
        KSPSetType(this->ksp, KSPMINRES);
    }

    //PetscViewerAndFormat *vf;
    //PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, &vf);
    //KSPMonitorSet(this->ksp, (PetscErrorCode (*)(KSP, PetscInt, PetscReal, void *))KSPMonitorResidual, vf, (PetscCtxDestroyFn *)PetscViewerAndFormatDestroy);
}

bool PETScPCG::solve_exec(std::vector<double>& load){
    int mpi_id = 0;
    int mpi_size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    auto K = this->gsm->get_K();

    long M = load.size();

    long begin = 0, end = 0;
    VecGetOwnershipRange(this->f, &begin, &end);

    double* f_data = nullptr;
    VecGetArray(this->f, &f_data);
    for(auto d = 0; d < end - begin; ++d){
        *(f_data+d) = load[begin+d];
    }
    VecRestoreArray(this->f, &f_data);

    VecAssemblyBegin(this->u);
    VecAssemblyEnd(this->u);

    VecAssemblyBegin(this->f);
    VecAssemblyEnd(this->f);

    //PetscBool is_sym;
    //MatIsSymmetric(K, 1e-7, &is_sym);
    //logger::quick_log("is symmetric?", is_sym);

    if(!this->setup){
        KSPSetOperators(this->ksp, K, K);

        KSPSetUp(this->ksp);
        PCSetUp(this->pc);
        this->setup = true;
    }

    KSPSolve(this->ksp, this->f, this->u);
    KSPConvergedReason r;
    KSPGetConvergedReason(this->ksp, &r);
    logger::quick_log("Converged?", r);
    {
        double tmp = 0;
        PetscInt total_its = 0;
        PetscInt its = 0;
        KSPGetResidualNorm(this->ksp, &tmp);
        KSPGetTotalIterations(this->ksp, &total_its);
        KSPGetIterationNumber(this->ksp, &its);
        logger::quick_log("Residual norm", tmp);
        logger::quick_log("Iterations", its);
        logger::quick_log("Total iterations", total_its);
    }

    if(r == KSP_CONVERGED_STEP_LENGTH){
        return false;
    }

    const double* load_data;

    VecGetArrayRead(u, &load_data);

    long n = 0, m = 0;
    MatGetLocalSize(K, &n, &m);
    MPI_Bcast(&M, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    if(mpi_size > 1){
        if(mpi_id == 0){
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

    VecRestoreArrayRead(this->u, &load_data);

    return true;
}

void PETScPCG::solve(std::vector<double>& load, const std::vector<double>& u0){
    long M = load.size();
    if(this->first_time){
        this->solve_init(M);

        KSPSetInitialGuessNonzero(this->ksp, PETSC_TRUE);
        
        this->first_time = false;
    }

    MPI_Bcast(load.data(), load.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    long begin = 0, end = 0;
    VecGetOwnershipRange(this->u, &begin, &end);
    double* u_data = nullptr;
    VecGetArray(this->u, &u_data);
    for(auto d = 0; d < end - begin; ++d){
        *(u_data+d) = u0[d + begin];
    }
    VecRestoreArray(this->u, &u_data);
    
    if(!this->solve_exec(load)){
        long begin = 0, end = 0;
        VecGetOwnershipRange(this->u, &begin, &end);
        double* u_data = nullptr;
        VecGetArray(this->u, &u_data);
        for(auto d = 0; d < end - begin; ++d){
            *(u_data+d) = 0;
        }
        VecRestoreArray(this->u, &u_data);

        this->solve_exec(load);
    }
}
void PETScPCG::solve(std::vector<double>& load){
    long M = load.size();
    if(this->first_time){
        this->solve_init(M);

        KSPSetInitialGuessNonzero(this->ksp, PETSC_FALSE);

        if(this->contact_data.contact_type == FiniteElement::ContactType::FRICTIONLESS_DISPL_LOG || this->nonlinear_boundary_conditions){
            KSPSetTolerances(this->ksp, 1e-4, 1e-50, 1e10, 1e5);
            KSPSetInitialGuessNonzero(this->ksp, PETSC_FALSE);
        } else if(this->contact_data.contact_type >= FiniteElement::ContactType::FRICTIONLESS_DISPL_SIMPLE){
            KSPSetTolerances(this->ksp, 1e-3, 1e-50, 1e10, 1e5);
            KSPSetInitialGuessNonzero(this->ksp, PETSC_FALSE);

            PCJacobiSetUseAbs(this->pc, PETSC_TRUE);
            PCJacobiSetFixDiagonal(this->pc, PETSC_TRUE);
            PCJacobiSetType(this->pc, PC_JACOBI_ROWMAX);

            KSPSetType(this->ksp, KSPMINRES);
        }

        this->first_time = false;
    }

    MPI_Bcast(load.data(), load.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    this->solve_exec(load);
}
void PETScPCG::reset_hessian(){
    this->gsm->reset_hessian();
    this->setup = false;
}

using namespace projspec;
const bool PETScPCG::reg = Factory<FiniteElement>::add(
    [](const DataMap& data){
        return std::make_unique<PETScPCG>(data);
    },
    ObjectRequirements{
        "petsc_pcg",
        {
            DataEntry{.name = "backend", .type = TYPE_STRING, .required = true}
        }
    }
);

}
