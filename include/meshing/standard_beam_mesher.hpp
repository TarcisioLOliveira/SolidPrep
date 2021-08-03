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

#ifndef STANDARD_BEAM_MESHER_HPP
#define STANDARD_BEAM_MESHER_HPP

#include "meshing.hpp"
#include "utils.hpp"

namespace meshing{

class StandardBeamMesher : public Meshing{
    public:
    StandardBeamMesher(double size, int order, utils::ProblemType type, int algorithm = 5);

    virtual std::vector<ElementShape> mesh(TopoDS_Shape s) override;

    private:
    double size;
    int order;
    int dim;
    int algorithm;

    MeshNode* find_node(size_t id) const;
};

}

#endif
