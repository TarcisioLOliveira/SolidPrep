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

#ifndef PETSC_PCG_HPP
#define PETSC_PCG_HPP

#include "finite_element.hpp"
#include "global_stiffness_matrix/petsc_sparse_symmetric.hpp"

namespace finite_element{

class PETScPCG : public FiniteElement{
    public:
    typedef global_stiffness_matrix::PETScSparseSymmetric::Backend PETScBackend;

    PETScPCG(PETScBackend backend);

    virtual ~PETScPCG();

    virtual std::vector<double> calculate_displacements(const Meshing* const mesh, std::vector<double> load, const std::vector<double>& density = std::vector<double>(), double pc = 3, double psi = 0.1) override;

    inline virtual void set_steps(size_t s) override{
        this->steps = s;
        this->u.resize(steps);
    }

    private:
    std::unique_ptr<global_stiffness_matrix::PETScSparseSymmetric> gsm;
    std::vector<Vec> u;
    Vec f = 0;
    PC pc = 0;
    KSP ksp = 0;
    bool first_time = true;
    std::string vec_type = VECSTANDARD;
};

}

#endif
