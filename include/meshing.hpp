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

#ifndef MESHING_HPP
#define MESHING_HPP

#include <TopoDS_Shape.hxx>
#include "element.hpp"
#include <vector>

struct ElementShape{
    std::vector<MeshNode*> nodes;
};

class Meshing{
    public:

    virtual std::vector<ElementShape> mesh(TopoDS_Shape s) = 0;

    protected:
    std::vector<std::unique_ptr<MeshNode>> node_list;
};

#endif
