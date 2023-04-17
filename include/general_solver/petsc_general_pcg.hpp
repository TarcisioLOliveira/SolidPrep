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

#ifndef PETSC_GENERAL_PCG_HPP
#define PETSC_GENERAL_PCG_HPP

#include <vector>
#include "general_global_matrix/petsc_global_sparse.hpp"

namespace general_solver{

class PETScGeneralPCG{
    public:
    ~PETScGeneralPCG();

    void initialize(general_global_matrix::PETScGlobalSparse* M, size_t L);

    void set_rhs(std::vector<double> b);

    void solve(std::vector<double>& x);

    private:
    general_global_matrix::PETScGlobalSparse* M;
    Vec f = 0;
    Vec u = 0;
    PC pc = 0;
    KSP ksp = 0;
    size_t L = 0;
};

}

#endif
