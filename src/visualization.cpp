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

#include "visualization.hpp"
#include "utils.hpp"
#include <gmsh.h>
#include <algorithm>
#include <vector>
#include "logger.hpp"

void Visualization::load_mesh(Meshing* mesh, utils::ProblemType type){
    gmsh::clear();
    gmsh::model::add(this->MODEL_NAME);

    this->mesh = mesh;
    this->type = type;

    std::vector<double> coords;
    coords.reserve(mesh->node_list.size()*3);
    std::vector<size_t> node_tags;
    node_tags.reserve(mesh->node_list.size());
    for(auto& n:mesh->node_list){
        coords.push_back(n->point.X());
        coords.push_back(n->point.Y());
        coords.push_back(n->point.Z());

        // It looks like node tags have to begin at 1 for some reason.
        node_tags.push_back(n->id+1);
    }

    if(type == utils::PROBLEM_TYPE_2D){
        this->mesh_tag = gmsh::model::addDiscreteEntity(2);
        gmsh::model::mesh::addNodes(2, this->mesh_tag, node_tags, coords);
    } else if(type == utils::PROBLEM_TYPE_3D){
        this->mesh_tag = gmsh::model::addDiscreteEntity(3);
        gmsh::model::mesh::addNodes(3, this->mesh_tag, node_tags, coords);
    }

    size_t s_size = 0;
    for(auto& g:mesh->geometries){
        s_size += g->mesh.size();
    }
    std::vector<size_t> elem_tags;
    elem_tags.reserve(s_size);
    std::vector<size_t> elem_nodes;
    elem_nodes.reserve(s_size*mesh->elem_info->get_nodes_per_element());
    size_t etag = 1;
    const size_t nodes_num = mesh->elem_info->get_nodes_per_element();
    for(auto& g:mesh->geometries){
        for(auto& e:g->mesh){
            elem_tags.push_back(etag++);
            for(size_t i = 0; i < nodes_num; ++i){
                const auto& n = e->nodes[i];
                elem_nodes.push_back(n->id+1);
            }
        }
    }
    gmsh::model::mesh::addElementsByType(this->mesh_tag, mesh->elem_info->get_gmsh_element_type(), elem_tags, elem_nodes);

    gmsh::option::setNumber("Mesh.SurfaceEdges", 0);
    gmsh::option::setNumber("Mesh.VolumeEdges", 0);
    gmsh::option::setNumber("Mesh.Triangles", 0);
    gmsh::option::setNumber("General.ColorScheme", 3);
    gmsh::option::setNumber("General.FltkColorScheme", 1);
    gmsh::option::setString("General.GraphicsFontEngine", "StringTexture");

}

ViewHandler* Visualization::add_view(const std::string& view_name, ViewHandler::ViewType view_type, ViewHandler::DataType data_type){
    ++this->last_view_tag;
    this->handler_list.emplace_back(std::make_unique<ViewHandler>(this->mesh, this->MODEL_NAME, view_name, view_type, data_type, this->type, this->last_view_tag));
    return this->handler_list.back().get();
}

void Visualization::show(){
    this->shown = true;

    gmsh::fltk::initialize();
}
void Visualization::wait(){
    auto checkForEvent = [=]() -> bool {
        std::vector<std::string> action;
        gmsh::onelab::getString("ONELAB/Action", action);
        if(action.size() && action[0] == "check") {
            gmsh::onelab::setString("ONELAB/Action", {""});
            gmsh::graphics::draw();
        }
        return true;
    };

    while(gmsh::fltk::isAvailable() && checkForEvent() && this->shown)
        gmsh::fltk::wait();
}
