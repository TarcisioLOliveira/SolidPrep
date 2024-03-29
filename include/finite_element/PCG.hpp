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

#ifndef PCG_HPP
#define PCG_HPP

#include "finite_element.hpp"
#include "global_stiffness_matrix/lapack_dense_symmetric_banded.hpp"

namespace finite_element{

class PCG : public FiniteElement{
    public:
    enum class Preconditioner{
        JACOBI,
        SSOR
    };

    PCG(const double eps, const Preconditioner precond);

    virtual void generate_matrix(const Meshing* const mesh, const size_t L, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache) override;

    virtual void calculate_displacements(std::vector<double>& load) override;

    private:
    bool first_time = true;
    bool setup = false;
    const double eps;
    const Preconditioner precond;
    std::vector<double> P;
    global_stiffness_matrix::LAPACKDenseSymmetricBanded gsm;

    void generate_P();
};

}

#endif
