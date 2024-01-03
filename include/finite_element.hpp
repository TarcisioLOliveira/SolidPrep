/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
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

#ifndef FINITE_ELEMENT_HPP
#define FINITE_ELEMENT_HPP

#include <memory>
#include <vector>
#include "element.hpp"
#include "meshing.hpp"
#include "geometry.hpp"

class ProjectData;

class FiniteElement{
    public:
    const double K_MIN = 1e-6;

    virtual ~FiniteElement() = default;

    virtual void generate_matrix(const Meshing* const mesh, const size_t L, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache) = 0;

    virtual void calculate_displacements(std::vector<double>& load) = 0;

    virtual std::vector<double> calculate_forces(const Meshing* const mesh, const std::vector<double>& displacements) const;
};

#endif
