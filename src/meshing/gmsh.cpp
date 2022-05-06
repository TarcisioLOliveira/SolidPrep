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

#include "meshing/gmsh.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include <gmsh.h>
#include <algorithm>
#include <limits>
#include <BOPAlgo_Splitter.hxx>
#include <BOPAlgo_Builder.hxx>
#include "project_data.hpp"
#include <BRepBuilderAPI_Copy.hxx>

namespace meshing{

Gmsh::Gmsh(double size, int order, utils::ProblemType type, ProjectData* data, int algorithm):
    Meshing(size), order(order), dim(0), algorithm(algorithm), data(data){
    if(type == utils::PROBLEM_TYPE_2D){
        dim = 2;
    } else if(type == utils::PROBLEM_TYPE_3D){
        dim = 3;
    }
}

std::vector<ElementShape> Gmsh::mesh(TopoDS_Shape s){
    this->shape = s;
    this->node_list.clear();

    bool has_condition_inside = false;

    TopoDS_Shape shape = BRepBuilderAPI_Copy(s);
    for(auto& f:this->data->forces){
        if(this->is_strictly_inside2D(f.S.get_centroid(), this->shape)){
            has_condition_inside = true;

            BOPAlgo_Splitter splitter;
            splitter.SetNonDestructive(true);
            splitter.AddArgument(shape);
            splitter.AddTool(f.S.get_shape());
            splitter.Perform();
            shape = splitter.Shape();
        }
    }
    for(auto& s:this->data->supports){
        if(this->is_strictly_inside2D(s.S.get_centroid(), this->shape)){
            has_condition_inside = true;

            BOPAlgo_Splitter splitter;
            splitter.SetNonDestructive(true);
            splitter.AddArgument(shape);
            splitter.AddTool(s.S.get_shape());
            splitter.Perform();
            shape = splitter.Shape();
        }
    }

    gmsh::initialize();

    gmsh::model::add("base");

    gmsh::vectorpair vec;
    gmsh::model::occ::importShapesNativePointer(static_cast<const void*>(&shape), vec);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMin", this->size);
    gmsh::option::setNumber("Mesh.MeshSizeMax", this->size);

    gmsh::option::setNumber("Mesh.Algorithm", this->algorithm);

    gmsh::option::setNumber("Mesh.ElementOrder", this->order);
    gmsh::option::setNumber("Mesh.HighOrderOptimize", 2);
    gmsh::option::setNumber("Mesh.Optimize", 1);
    gmsh::option::setNumber("Mesh.OptimizeNetgen", 1);

    gmsh::model::mesh::generate(this->dim);

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, nodeParams;
    if(has_condition_inside){
        gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, -1, -1, true);
    } else {
        gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, this->dim, -1, true);
    }

    std::vector<int> elemTypes;
    std::vector<std::vector<std::size_t> > elemTags, elemNodeTags;
    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, this->dim, -1);

    this->node_list.reserve(nodeTags.size());
    if(this->dim == 2){
        for(size_t i = 0; i < nodeTags.size(); ++i){
            gp_Pnt p(nodeCoords[i*3], nodeCoords[i*3+1], nodeCoords[i*3+2]);
            this->node_list.emplace_back(MeshNodeFactory::make_node(p, nodeTags[i], MeshNodeFactory::MESH_NODE_2D));
        }
        // gp_Pnt p(nodeCoords[0], nodeCoords[1], nodeCoords[2]);
        // this->node_list.emplace_back(MeshNodeFactory::make_node(p, nodeTags[0], MeshNodeFactory::MESH_NODE_2D)); 
        // for(size_t i = 1; i < nodeTags.size(); ++i){
        //     gp_Pnt p = gp_Pnt(nodeCoords[i*3], nodeCoords[i*3+1], nodeCoords[i*3+2]);
        //     auto get_id = [&p](const std::unique_ptr<MeshNode>& m)->bool{ return p.IsEqual(m->point, Precision::Confusion()); };
        //     auto it = std::find_if(this->node_list.begin(), this->node_list.end(), get_id);
        //     if(it == this->node_list.end()){
        //         this->node_list.emplace_back(MeshNodeFactory::make_node(p, nodeTags[i], MeshNodeFactory::MESH_NODE_2D)); 
        //     }
        // }
    } else {
        // for(size_t i = 0; i < nodeTags.size(); ++i){
        //     gp_Pnt p(nodeCoords[i*3], nodeCoords[i*3+1], nodeCoords[i*3+2]);
        //     this->node_list.emplace_back(MeshNodeFactory::make_node(p, nodeTags[i], MeshNodeFactory::MESH_NODE_3D)); 
        // }
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
    for(auto n:elemNodeTags[0]){
        // gp_Pnt p(nodeCoords[n*3], nodeCoords[n*3+1], nodeCoords[n*3+2]);
        // auto get_id = [&p](const std::unique_ptr<MeshNode>& m)->bool{ return p.IsEqual(m->point, Precision::Confusion()); };
        auto get_id = [n](const std::unique_ptr<MeshNode>& m)->bool{ return n == m->id; };
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
    // Prune unused nodes
    size_t new_id = 0;
    for(auto& n : this->node_list){
        n->id = new_id;
        ++new_id;
    }

    std::vector<bool> used(this->node_list.size(), false);
    for(size_t i = 0; i < list.size(); ++i){
        for(auto& n:list[i].nodes){
            used[n->id] = true;
        }   
    }

    bool renumber = false;
    auto it = this->node_list.begin();
    while(it < this->node_list.end()){
        if(!used[(*it)->id]){
            it = this->node_list.erase(it);
            renumber = true;
        } else {
            ++it;
        }
    }
    used.clear();

    if(renumber){
        new_id = 0;
        for(auto& n : this->node_list){
            n->id = new_id;
            ++new_id;
        }
    }
    this->reverse_cuthill_mckee(list);

    gmsh::clear();
    gmsh::finalize();

    return list;
}

}
