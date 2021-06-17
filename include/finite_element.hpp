/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   Foobar is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Foobar is distributed in the hope that it will be useful,
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

class FiniteElement{
    public:

    virtual std::vector<float> calculate_displacements(const std::vector<MeshElement*>& mesh, const std::vector<float>& loads, const std::vector<float>& density = std::vector<float>()) const = 0;

    virtual void calculate_stresses(const std::vector<MeshElement*>& mesh, const std::vector<float>& displacements) const;
    virtual void calculate_forces(const std::vector<MeshElement*>& mesh, const std::vector<float>& displacements) const;
};

#endif
