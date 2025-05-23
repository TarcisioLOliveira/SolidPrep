/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
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

#ifndef MUMPS_SOLVER_HPP
#define MUMPS_SOLVER_HPP

#include <dmumps_c.h>
#include "finite_element.hpp"
#include "global_stiffness_matrix/mumps_sparse_symmetric.hpp"
#include "math/matrix.hpp"
#include "project_specification/data_map.hpp"

// Recommended by the documentation
#define ICNTL( i ) icntl[ (i) - 1 ]
#define INFO( i ) info[ (i) - 1 ]

namespace finite_element{

class MUMPSSolver : public FiniteElement{
    public:
    MUMPSSolver(const projspec::DataMap& data);
    virtual ~MUMPSSolver();

    private:
    static const bool reg;
    virtual void generate_matrix_base(const Meshing* const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const ContactType type) override;

    virtual void solve(std::vector<double>& load) override;
    virtual void reset_hessian() override;

    DMUMPS_STRUC_C config;
    global_stiffness_matrix::MUMPSSparseSymmetric gsm;
    std::vector<double> buffer;
    size_t u_size = 0;
    size_t l_num = 0;
    bool first_time = true;
    bool factorized = false;
};

}

#endif
