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

#include <vector>
#include "element.hpp"
#include "meshing.hpp"

class ProjectData;

class FiniteElement{
    public:

    virtual std::vector<double> calculate_displacements(ProjectData* data, Meshing* mesh, const std::vector<double>& density = std::vector<double>(), double pc = 3, bool use_stored_matrix = false, const std::vector<double>& virtual_load = std::vector<double>()) = 0;

    virtual void calculate_stresses(Meshing* mesh, const std::vector<double>& displacements, const std::vector<double>& density = std::vector<double>()) const;
    virtual void calculate_forces(Meshing* mesh, const std::vector<double>& displacements) const;
};

#endif
