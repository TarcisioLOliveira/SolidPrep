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

#ifndef MUMPS_SPARSE_SYMMETRIC_HPP
#define MUMPS_SPARSE_SYMMETRIC_HPP

#include "meshing.hpp"
#include "utils/sparse_matrix.hpp"

namespace global_stiffness_matrix{

class MUMPSSparseSymmetric{
    public:
    const double K_MIN = 1e-6;

    virtual ~MUMPSSparseSymmetric() = default;

    virtual void generate(const Meshing * const mesh, const std::vector<double>& density, const double pc);

    inline std::vector<int>& get_rows(){
        return rows;
    }
    inline std::vector<int>& get_cols(){
        return cols;
    }
    inline std::vector<double>& get_vals(){
        return vals;
    }

    protected:
    utils::SparseMatrix sK;
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> vals;

    virtual void add_geometry(const Meshing * const mesh, const Geometry * const g);

    virtual void add_geometry(const Meshing * const mesh, const Geometry * const g, std::vector<double>::const_iterator& rho, const double pc);
};

}

#endif
