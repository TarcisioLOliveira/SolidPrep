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

#ifndef CUSOLVER_IMPL_HPP
#define CUSOLVER_IMPL_HPP

#include <cusparse.h>
#include <cusolverSp.h>
#include <cusolverSp_LOWLEVEL_PREVIEW.h>

#include "finite_element.hpp"
#include "global_stiffness_matrix/cusolver_sparse_symmetric.hpp"
#include "utils/cuda_array.hpp"

namespace finite_element::cuda{

// Based on: https://github.com/gishi523/cusparse-cholesky-solver
//
// Usses low level Cholesky functions which are undocumented in the latest
// version of the documentation. Documentation on them can be found in the
// 2017 version of the docs.
//
// Hopefully it will last for a while. Although it does not seem to be as
// fast as I'd like (at least on my GTX 1650 Mobile).
class cuSolverImpl : public FiniteElement{
    public:
    cuSolverImpl();
    virtual ~cuSolverImpl();
    virtual std::vector<double> calculate_displacements(const Meshing* const mesh, std::vector<double> load, const std::vector<double>& density = std::vector<double>(), double pc = 3, double psi = 0.1) override;

    private:
    global_stiffness_matrix::cuda::cuSolverSparseSymmetric gsm;
    cusolverSpHandle_t solver_handle;
    cusparseMatDescr_t mat_descr;
    csrcholInfo_t info;
    utils::cuda::CUDAArray<unsigned char> workspace;
    utils::cuda::CUDAArray<double> x;
    utils::cuda::CUDAArray<double> b;
    bool first_time = true;

    void cusolver_status(cusolverStatus_t status, std::string function) const;

};


}

#endif
