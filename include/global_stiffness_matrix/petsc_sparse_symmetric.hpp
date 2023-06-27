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

#ifndef PETSC_SPARSE_SYMMETRIC_HPP
#define PETSC_SPARSE_SYMMETRIC_HPP

#include <vector>
#include <petsc.h>
#include "meshing.hpp"
#include "global_stiffness_matrix.hpp"

namespace global_stiffness_matrix{

class PETScSparseSymmetric : public GlobalStiffnessMatrix{
    public:
    const double K_MIN = 1e-6;
    enum class Backend{
        CPU,
        CUDA
    };

    PETScSparseSymmetric() = default;
    PETScSparseSymmetric(Backend backend);

    virtual ~PETScSparseSymmetric();

    virtual void generate(const Meshing * const mesh, const std::vector<double>& density, const double pc, const double psi) override;

    Mat get_K() const{
        return this->K;
    }

    protected:
    Mat K = 0;
    bool first_time = true;
    std::string mat_type = MATMPIAIJ;

    inline virtual void insert_element_matrix(const std::vector<double>& k, const std::vector<long>& pos) override{
        // Requires 64-bit indices
        MatSetValues(this->K, pos.size(), pos.data(), pos.size(), pos.data(), k.data(), ADD_VALUES);
    }

    inline virtual void add_to_matrix(size_t i, size_t j, double val) override{
        // Requires 64-bit indices
        MatSetValue(this->K, i, j, val, ADD_VALUES);
    }
};

}

#endif
