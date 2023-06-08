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

#ifndef GLOBAL_STIFFNESS_MATRIX_HPP
#define GLOBAL_STIFFNESS_MATRIX_HPP

#include "meshing.hpp"

class GlobalStiffnessMatrix{
    public:
    virtual ~GlobalStiffnessMatrix() = default;
    const double K_MIN = 1e-6;
    virtual void generate(const Meshing * const mesh, const std::vector<double>& density, const double pc, const double psi) = 0;

    protected:
    size_t W, N;

    virtual void generate_base(const Meshing * const mesh, const std::vector<double>& density, const double pc, const double psi);

    virtual void calculate_dimensions(const Meshing * const mesh, const std::vector<double>& load);

    virtual void add_geometry(const Meshing * const mesh, const Geometry * const g);

    virtual void add_geometry(const Meshing * const mesh, const Geometry * const g, std::vector<double>::const_iterator& rho, const double pc, const double psi);

    virtual void insert_element_matrix(const std::vector<double>& k, const std::vector<long>& pos) = 0;
};

#endif
