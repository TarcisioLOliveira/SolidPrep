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
#include "math/matrix.hpp"
#include "project_specification/data_map.hpp"

namespace finite_element{

class PETScPCG : public FiniteElement{
    public:
    typedef global_stiffness_matrix::PETScSparseSymmetric::Backend PETScBackend;

    PETScPCG(const projspec::DataMap& data);

    virtual ~PETScPCG();

    private:
    const static bool reg;
    virtual void generate_matrix_base(const Meshing* const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const ContactType type) override;

    virtual void solve(std::vector<double>& load) override;
    virtual void reset_hessian() override;

    std::unique_ptr<global_stiffness_matrix::PETScSparseSymmetric> gsm;
    Vec u;
    Vec f = 0;
    PC pc = 0;
    KSP ksp = 0;
    bool first_time = true;
    bool setup = false;
    std::string vec_type = VECSTANDARD;
};

}

#endif
