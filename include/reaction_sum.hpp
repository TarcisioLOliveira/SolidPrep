/*
 *   Copyright (C) 2026 Tarcísio Ladeia de Oliveira.
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

#ifndef REACTION_SUM_HPP
#define REACTION_SUM_HPP

#include <array>
#include <unordered_map>
#include "math/matrix.hpp"
#include "utils.hpp"
#include "cross_section.hpp"
#include "element_factory.hpp"
#include "utils/delayed_pointer.hpp"

class Meshing;

class ReactionSum{
    public:

    // ASSUMES CROSS SECTION IS PLANE
    // ASSUMES this->normal == this->S.normal

    ReactionSum(CrossSection cross_section, double thickness, std::array<double, 3> F, utils::ProblemType type);
    ReactionSum(ReactionSum&&) = default;

    void initialize(const Meshing* const mesh, size_t& max_dofs);
    void apply_constraints(std::vector<double>& load) const;

    math::Matrix get_constr(const Meshing* const mesh, const BoundaryElement* b) const;

    const CrossSection S;
    const double thickness;
    const std::array<double, 3> F;
    const size_t F_dim;
    std::vector<const BoundaryElement*> submesh;

    size_t offset;

    private:
    const utils::ProblemType type;
    const MeshElementFactory* elem_info;
    size_t BOUNDARY_NODE_NUM;
    std::unordered_map<size_t, bool> affected_nodes;
};

#endif
