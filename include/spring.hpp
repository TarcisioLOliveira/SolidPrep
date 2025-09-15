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

#ifndef SPRING_HPP
#define SPRING_HPP

#include <array>
#include "element_factory.hpp"
#include "material.hpp"
#include "math/matrix.hpp"
#include "utils.hpp"
#include "cross_section.hpp"
#include "utils/delayed_pointer.hpp"

class BoundaryElement;

class Spring{
    public:

    // ASSUMES CROSS SECTION IS PLANE
    // ASSUMES this->normal == this->S.normal

    Spring(CrossSection cross_section, double thickness, gp_Dir normal, gp_Dir v, gp_Dir w, utils::DelayedPointerView<Material> mat, std::array<double, 3> L, MeshElementFactory* elem, utils::ProblemType type);
    Spring(Spring&&) = default;
    void generate_mesh(std::vector<BoundaryElement>& boundary_elements);

    math::Matrix get_K(const MeshElement* const e, const gp_Pnt& p) const;

    const CrossSection S;
    const double A;
    const double thickness;
    const math::Matrix rot2D;
    const math::Matrix rot3D;
    utils::DelayedPointerView<Material> mat;
    const gp_Dir normal;
    const MeshElementFactory* elem_info;
    std::vector<const BoundaryElement*> submesh;

    private:
    const gp_Dir v;
    const gp_Dir w;
    const std::array<double, 3> L;
    const utils::ProblemType type;
};

#endif
