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

#include "meshing/standard_beam_mesher.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include <gmsh.h>
#include <algorithm>
#include <limits>

namespace meshing{

StandardBeamMesher::StandardBeamMesher(double size, int order, utils::ProblemType type, int algorithm):
    BeamMeshing(size), order(order), dim(0), algorithm(algorithm){
    if(type == utils::PROBLEM_TYPE_2D){
        dim = 2;
    } else if(type == utils::PROBLEM_TYPE_3D){
        dim = 3;
    }
}

std::vector<ElementShape> StandardBeamMesher::mesh(TopoDS_Shape s){
    this->shape = s;
    this->node_list.clear();
    gmsh::initialize();

    gmsh::model::add("base");

    gmsh::vectorpair vec;
    gmsh::model::occ::importShapesNativePointer(static_cast<const void*>(&s), vec);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMin", this->size);
    gmsh::option::setNumber("Mesh.MeshSizeMax", this->size);

    gmsh::option::setNumber("Mesh.Algorithm", this->algorithm);

    gmsh::option::setNumber("Mesh.ElementOrder", this->order);
    gmsh::option::setNumber("Mesh.HighOrderOptimize", 2);

    gmsh::model::mesh::generate(this->dim);

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, nodeParams;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, this->dim, -1, true);

    std::vector<std::size_t> boundNodeTags;
    std::vector<double> boundNodeCoords, boundNodeParams;
    gmsh::model::mesh::getNodes(boundNodeTags, boundNodeCoords, boundNodeParams, this->dim-1, -1, true);

    std::vector<int> elemTypes;
    std::vector<std::vector<size_t> > elemTags, elemNodeTags;
    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, this->dim, -1);

    std::vector<int> boundElemTypes;
    std::vector<std::vector<size_t> > boundElemTags, boundElemNodeTags;
    gmsh::model::mesh::getElements(boundElemTypes, boundElemTags, boundElemNodeTags, this->dim-1, -1);

    this->node_list.reserve(nodeTags.size());
    if(this->dim == 2){
        for(size_t i = 0; i < nodeTags.size(); ++i){
            gp_Pnt p(nodeCoords[i*3], nodeCoords[i*3+1], nodeCoords[i*3+2]);
            this->node_list.emplace_back(MeshNodeFactory::make_node(p, nodeTags[i], MeshNodeFactory::MESH_NODE_2D)); 
        }
    } else {
        // for(size_t i = 0; i < nodeTags.size(); ++i){
        //     gp_Pnt p(nodeCoords[i*3], nodeCoords[i*3+1], nodeCoords[i*3+2]);
        //     this->node_list.emplace_back(MeshNodeFactory::make_node(p, nodeTags[i], MeshNodeFactory::MESH_NODE_3D)); 
        // }
    }
    for(size_t i = 0; i < boundNodeTags.size(); ++i){
        auto get_id = [&boundNodeTags, i](const std::unique_ptr<MeshNode>& m)->bool{ return boundNodeTags[i] == m->id; };
        MeshNode* node = std::find_if(this->node_list.begin(), this->node_list.end(), get_id)->get();
        auto get_id2 = [&boundNodeTags, i](const BoundaryNode& m)->bool{ return boundNodeTags[i] == m.node->id; };
        if(std::find_if(this->boundary_nodes.begin(), this->boundary_nodes.end(), get_id2) == this->boundary_nodes.end()){
            boundary_nodes.push_back(node);
        }
    }
    std::vector<std::vector<size_t>> neighbors(boundary_nodes.size());
    size_t nodes_per_bound_elem = 0;
    if(this->dim == 2){
        nodes_per_bound_elem = 2;
    } else if(this->dim == 3){
        nodes_per_bound_elem = 3;
    }

    size_t idx = 0;
    std::vector<size_t> elem;
    for(auto& n:boundElemNodeTags[0]){
        auto get_id = [&n](const BoundaryNode& m)->bool{ return n == m.node->id; };
        size_t node_id = std::find_if(this->boundary_nodes.begin(), this->boundary_nodes.end(), get_id) - this->boundary_nodes.begin();
        elem.push_back(node_id);
        ++idx;

        if(idx == nodes_per_bound_elem){
            idx = 0;
            for(size_t i = 0; i < elem.size(); ++i){
                for(size_t j = i+1; j < elem.size(); ++j){
                    if(elem[i] != elem[j]){
                        neighbors[elem[i]].push_back(elem[j]);
                        neighbors[elem[j]].push_back(elem[i]);
                    }
                }
            }
            elem.clear();
        }
    }
    for(size_t i = 0; i < neighbors.size(); ++i){
        gp_Vec vec(0,0,0);
        for(auto& n:neighbors[i]){
            vec += gp_Vec(boundary_nodes[n].node->point, boundary_nodes[i].node->point);
        }
        // If colinear/coplanar
        if(vec.IsEqual(gp_Vec(0, 0, 0), Precision::Confusion(), Precision::Angular())){
            if(this->dim == 2){
                gp_Pnt p1 = boundary_nodes[neighbors[i][0]].node->point;
                gp_Pnt p2 = boundary_nodes[neighbors[i][1]].node->point;
                gp_Pnt p3 = boundary_nodes[i].node->point;
                vec = gp_Vec(p1, p2);
                gp_Ax1 ax(p3, gp_Dir(0, 0, 1));
                vec.Rotate(ax, M_PI/2);
            } else if(this->dim == 3){
                gp_Pnt p1 = boundary_nodes[neighbors[i][0]].node->point;
                gp_Pnt p2 = boundary_nodes[neighbors[i][1]].node->point;
                gp_Pnt p3 = boundary_nodes[i].node->point;
                gp_Vec v1(p3, p1);
                gp_Vec v2(p3, p2);
                vec = v1.Crossed(v2);
            }
        }
        gp_Dir dir(vec);
        gp_Pnt p = boundary_nodes[i].node->point.Translated(dir);
        bool outside = true;
        if(this->dim == 2){
            outside = !this->is_inside_2D(p, s);
        } else if(this->dim == 3){
            outside = !this->is_inside_3D(p, s);
        }
        if(outside){
            boundary_nodes[i].normal = std::move(dir);
        } else {
            dir.Reverse();
            boundary_nodes[i].normal = std::move(dir);
        }
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

    std::vector<ElementShape> list;
    list.reserve(elemTags.size());
    list.emplace_back();
    int i = 0;
    for(auto& n:elemNodeTags[0]){
        auto get_id = [&n](const std::unique_ptr<MeshNode>& m)->bool{ return n == m->id; };
        MeshNode* node = std::find_if(this->node_list.begin(), this->node_list.end(), get_id)->get();
        list.back().nodes.push_back(node);
        ++i;

        if(i == node_per_elem){
            auto& nodes = list.back().nodes;
            double Delta = 0;

            if(node_per_elem % 3 == 0){
                gp_Pnt p[3] = {nodes[0]->point, nodes[1]->point, nodes[2]->point};
                gp_Mat deltaM(1, p[0].X(), p[0].Y(), 1, p[1].X(), p[1].Y(), 1, p[2].X(), p[2].Y());
                Delta = std::abs(deltaM.Determinant());
            } else {
                // TODO: Checking for 3D
            }
            if(Delta < 1e-3){
                list.pop_back();
            }

            list.emplace_back();
            i = 0;
        }
    }
    list.pop_back();

    for(size_t n = 0; n < this->node_list.size(); ++n){
        this->node_list[n]->id = std::numeric_limits<size_t>::max();
    }

    size_t new_id = 0;
    for(auto& n : list[0].nodes){
        n->id = new_id;
        ++new_id;
    }
    bool changed = true;
    while(changed){
        changed = false;
        for(auto& e : list){
            bool any = false;
            bool all = true;
            for(auto& n : e.nodes){
                if(n->id < std::numeric_limits<size_t>::max()){
                    any = true;
                } else if(n->id == std::numeric_limits<size_t>::max()){
                    all = false;
                }
            }
            if(any && !all){
                for(auto& n : e.nodes){
                    if(n->id == std::numeric_limits<size_t>::max()){
                        n->id = new_id;
                        ++new_id;
                        changed = true;
                    }
                }
            }
        }
    }

    gmsh::clear();
    gmsh::finalize();

    return list;
}

bool StandardBeamMesher::is_inside_2D(gp_Pnt p, const TopoDS_Shape& t){
    BRepClass3d_SolidClassifier insider(t, p, 0.01);
    return insider.State() == TopAbs_ON;
}

bool StandardBeamMesher::is_inside_3D(gp_Pnt p, const TopoDS_Shape& t){
    BRepClass3d_SolidClassifier insider(t, p, 0.01);
    return insider.State() == TopAbs_IN;
}

}