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

#ifndef GROUND_STRUCTURE_HPP
#define GROUND_STRUCTURE_HPP

#include <string>
#include <gp_Pnt.hxx>
#include <TopoDS_Shape.hxx>
#include "utils.hpp"

class GroundStructure{
    public:
    GroundStructure(const std::string& path, double scale, utils::ProblemType type);

    bool is_inside(const gp_Pnt& p) const;

    const TopoDS_Shape shape;
    const double scale;
    private:
    utils::ProblemType type;

    TopoDS_Shape load_shape(const std::string& path, double scale) const;
    bool is_inside_2D(const gp_Pnt& p) const;
    bool is_inside_3D(const gp_Pnt& p) const;
};

#endif
