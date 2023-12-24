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

#ifndef GRADIENT_DESCENT_HPP
#define GRADIENT_DESCENT_HPP

#include "finite_element.hpp"
#include "global_stiffness_matrix/lapack_dense_symmetric_banded.hpp"

namespace finite_element{

class GradientDescent : public FiniteElement{
    public:
    enum class Solver{
        STANDARD,
        MMA,
        LAGRANGE_MMA
    };

    GradientDescent(const double eps, Solver solver);

    inline virtual void set_steps(size_t s) override{
        this->steps = s;
        this->displacement.resize(s);
    }

    virtual std::vector<double> calculate_displacements(const Meshing* const mesh, const std::vector<long>& node_positions, std::vector<double> load, const std::vector<double>& density = std::vector<double>(), double pc = 3, double psi = 0.1) override;

    private:
    bool first_time = true;
    const double eps;
    std::vector<std::vector<double>> displacement;
    Solver solver;
    global_stiffness_matrix::LAPACKDenseSymmetricBanded gsm;

};

}

#endif
