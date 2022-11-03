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
#include "utils/sparse_matrix.hpp"

// Recommended by the documentation
#define ICNTL( i ) icntl[ (i) - 1 ]

namespace finite_element{

class MUMPSSolver : public FiniteElement{
    public:
    MUMPSSolver();
    virtual ~MUMPSSolver() = default;
    virtual std::vector<double> calculate_displacements(const Meshing* const mesh, std::vector<double> load, const std::vector<double>& density = std::vector<double>(), double pc = 3) override;

    private:
    utils::SparseMatrix sK;
    std::vector<double> rows;
    std::vector<double> cols;
    std::vector<double> vals;
    DMUMPS_STRUC_C config;

    virtual void add_geometry_to_K(const Meshing * const mesh, const Geometry * const g) override;

    virtual void add_geometry_to_K(const Meshing * const mesh, const Geometry * const g, std::vector<double>::const_iterator& rho, const double pc) override;
};

}

#endif
