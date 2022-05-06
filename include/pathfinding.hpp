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
#ifndef PATHFINDING_HPP
#define PATHFINDING_HPP

#include <gp_Pnt.hxx>
#include <vector>
#include <TopoDS_Shape.hxx>
#include <cross_section.hpp>

class Pathfinding{
    public:
    virtual ~Pathfinding() = default;

    virtual std::vector<gp_Pnt> find_path(const CrossSection& begin, const CrossSection& end) = 0;
};

#endif
