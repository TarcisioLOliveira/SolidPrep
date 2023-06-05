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

#ifndef CUSOLVER_SPARSE_SYMMETRIC_HPP
#define CUSOLVER_SPARSE_SYMMETRIC_HPP

#include "utils/csr.hpp"
#include "utils/cuda_array.hpp"
#include "meshing.hpp"

namespace global_stiffness_matrix::cuda{

class cuSolverSparseSymmetric{
    public:
    const double K_MIN = 1e-6;

    cuSolverSparseSymmetric();

    virtual ~cuSolverSparseSymmetric() = default;

    virtual void generate(const Meshing * const mesh, const std::vector<double>& density, const double pc, const double psi);

    inline int get_n() const{
        return this->n;
    }
    inline int get_nnz() const{
        return this->nnz;
    }
    inline int* get_rows() const{
        return this->csrRowPtr.pointer();
    }
    inline int* get_cols() const{
        return this->csrColInd.pointer();
    }
    inline double* get_vals() const{
        return this->csrVal.pointer();
    }


    protected:
    utils::CSR sK;
    int n;
    int nnz;
    utils::cuda::CUDAArray<int> csrRowPtr;
    utils::cuda::CUDAArray<int> csrColInd;
    utils::cuda::CUDAArray<double> csrVal;
    bool first_time = true;

    virtual void add_geometry(const Meshing * const mesh, const Geometry * const g);

    virtual void add_geometry(const Meshing * const mesh, const Geometry * const g, std::vector<double>::const_iterator& rho, const double pc, const double psi);
};

}

#endif
