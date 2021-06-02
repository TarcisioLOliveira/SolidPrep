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

#include "meshing/gmsh.hpp"
#include "utils.hpp"
#include <gmsh.h>

namespace meshing{
Gmsh::Gmsh(double size, int order, utils::ProblemType type, int algorithm):
    size(size), order(order), dim(0), algorithm(algorithm){
    if(type == utils::PROBLEM_TYPE_2D){
        dim = 2;
    } else if(type == utils::PROBLEM_TYPE_3D){
        dim = 3;
    }
}

std::vector<ElementShape> Gmsh::mesh(TopoDS_Shape s){
    gmsh::initialize();

    gmsh::model::add("base");

    gmsh::vectorpair vec;
    gmsh::model::occ::importShapesNativePointer(static_cast<const void*>(&s), vec);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeFactor", this->size);

    gmsh::option::setNumber("Mesh.Algorithm", this->algorithm);

    gmsh::option::setNumber("Mesh.ElementOrder", this->order);
    gmsh::option::setNumber("Mesh.HighOrderOptimize", this->order);

    gmsh::model::mesh::generate(this->dim);

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, nodeParams;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, -1, -1, true);

    std::vector<int> elemTypes;
    std::vector<std::vector<std::size_t> > elemTags, elemNodeTags;
    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, -1, -1);

    for(auto n:nodeTags){
        gp_Pnt p(nodeCoords[n*3], nodeCoords[n*3+1], nodeCoords[n*3+2]);
        this->node_list.emplace_back(MeshNodeFactory::make_node(p, n, MeshNodeFactory::MESH_NODE_2D)); 
    }

    int node_per_elem = 0;
    if(this->dim == 2){
        if(this->order == 1){
            node_per_elem = 3;
        } else if(this->order == 2){
            node_per_elem = 6;
        } else {
            node_per_elem = 12;
        }
    } else if(this->dim == 3){
        if(this->order == 1){
            node_per_elem = 4;
        } else {
            node_per_elem = 10;
        }
    }

    // Gets only tris/quads
    auto& e = elemNodeTags[this->dim-1];
    std::vector<ElementShape> list;
    int i = node_per_elem;
    for(auto n:e){
        if(i == node_per_elem){
            list.emplace_back();
            i = 0;
        }
        list.back().nodes.push_back(this->node_list[n].get());
        ++i;
    }

    gmsh::clear();
    gmsh::finalize();

    return list;
}

}
