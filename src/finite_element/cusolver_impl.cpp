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

#include "finite_element/cusolver_impl.hpp"
#include "logger.hpp"
#include "utils/cuda_array.hpp"

namespace finite_element::cuda{
void cuSolverImpl::cusolver_status(cusolverStatus_t status, std::string function) const{
    if(status == CUSOLVER_STATUS_NOT_INITIALIZED){
        logger::log_assert(false, logger::ERROR, "cuSolver not initialized.");
    } else if(status == CUSOLVER_STATUS_ALLOC_FAILED){
        logger::log_assert(false, logger::ERROR, "function {}() could not allocate necessary memory.", function);
    } else if(status == CUSOLVER_STATUS_ARCH_MISMATCH){
        logger::log_assert(false, logger::ERROR, "the device only supports compute capability 2.0 and above.");
    } else if(status == CUSOLVER_STATUS_INVALID_VALUE){
        logger::log_assert(false, logger::ERROR, "invalid parameters were passed to function {}().", function);
    } else if(status == CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED){
        logger::log_assert(false, logger::ERROR, "the matrix type is not supported.");
    } else if(status == CUSOLVER_STATUS_INTERNAL_ERROR){
        logger::log_assert(false, logger::ERROR, "an internal operation failed in function {}().", function);
    } else if(status != CUSOLVER_STATUS_SUCCESS){
        logger::log_assert(false, logger::ERROR, "cusparse returned error in function {}(), code {}.", function, status);
    }
}

cuSolverImpl::cuSolverImpl(){
    cudaSetDevice(0);
    this->cusolver_status(cusolverSpCreate(&this->solver_handle), "cusolverSpCreate");
    cusparseCreateMatDescr(&this->mat_descr);
    cusparseSetMatType(this->mat_descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(this->mat_descr, CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatDiagType(this->mat_descr, CUSPARSE_DIAG_TYPE_NON_UNIT);
    cusparseSetMatFillMode(this->mat_descr, CUSPARSE_FILL_MODE_LOWER);
    cusolverSpCreateCsrcholInfo(&this->info);
}

cuSolverImpl::~cuSolverImpl(){
    cusolverSpDestroyCsrcholInfo(this->info);
    cusparseDestroyMatDescr(this->mat_descr);
    cusolverSpDestroy(this->solver_handle);
}

std::vector<double> cuSolverImpl::calculate_displacements(const Meshing* const mesh, std::vector<double> load, const std::vector<double>& density, double pc, double psi){
    this->x.set_zero();
    if(this->b.get_size() > 0){
        this->b.set_data(load);
    }
    if(this->current_step == 0){
        this->gsm.generate(mesh, density, pc, psi);

        this->workspace.set_zero();
        if(this->first_time){
            this->first_time = false;
            size_t data = 0;
            size_t workspace_size = 0;
            this->cusolver_status(cusolverSpXcsrcholAnalysis(this->solver_handle, gsm.get_n(), gsm.get_nnz(), this->mat_descr, gsm.get_rows(), gsm.get_cols(), this->info), "cusolverSpXcsrcholAnalysis");
            this->cusolver_status(cusolverSpDcsrcholBufferInfo(this->solver_handle, gsm.get_n(), gsm.get_nnz(), this->mat_descr,
                gsm.get_vals(), gsm.get_rows(), gsm.get_cols(), this->info, &data, &workspace_size), "cusolverSpDcsrcholBufferInfo");
            this->workspace.allocate(workspace_size, 0);
            this->x.allocate(load.size(), 0);
            this->b.allocate(load);
        }

        this->cusolver_status(cusolverSpDcsrcholFactor(this->solver_handle, gsm.get_n(), gsm.get_nnz(), this->mat_descr,
                    gsm.get_vals(), gsm.get_rows(), gsm.get_cols(), this->info, static_cast<void*>(this->workspace.pointer())), "cusolverSpDcsrcholFactor");
    }

    this->cusolver_status(cusolverSpDcsrcholSolve(this->solver_handle, load.size(), this->b.pointer(), 
                    this->x.pointer(), this->info, static_cast<void*>(this->workspace.pointer())), "cusolverSpDcsrcholSolve");
    std::vector<double> x_host(load.size());
    this->x.get_data(x_host);

    this->current_step = (this->current_step + 1) % this->steps;

    return x_host;
}

}
