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

#ifndef GMSH_HPP
#define GMSH_HPP

#include "meshing.hpp"
#include "utils.hpp"

class ProjectData;

namespace meshing{

class Gmsh : public Meshing{
    public:
    Gmsh(const std::vector<std::unique_ptr<Geometry>>& geometries,
         const MeshElementFactory* const elem_type,
         double size, int order, utils::ProblemType type,
         int algorithm = 6);

    virtual void mesh(const std::vector<Force>& forces, 
                      const std::vector<Support>& supports,
                      const double thickness) override;

    private:
    double size;
    int order;
    int dim;
    int algorithm;

    void gmsh_meshing(bool has_condition_inside, TopoDS_Shape sh, std::vector<size_t>& elem_tags, std::vector<size_t>& elem_node_tags, const MeshElementFactory* const elem_type);
};

}

#endif
