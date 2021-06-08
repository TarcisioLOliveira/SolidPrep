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

#include "visualization.hpp"
#include <gmsh.h>

void Visualization::load_mesh(Meshing* mesh){
    gmsh::clear();
    gmsh::model::add(this->MODEL_NAME);

    this->mesh = mesh;

    std::vector<double> coords;
    coords.reserve(mesh->node_list.size()*3);
    std::vector<size_t> node_tags;
    node_tags.reserve(mesh->node_list.size());
    for(auto& n:mesh->node_list){
        coords.push_back(n->point.X());
        coords.push_back(n->point.Y());
        coords.push_back(n->point.Z());
        node_tags.push_back(n->id);
    }
    gmsh::model::mesh::addNodes(-1, -1, node_tags, coords);

    std::vector<size_t> types;
    types.reserve(mesh->element_list.size());
    std::vector<size_t> elem_tags;
    elem_tags.reserve(mesh->element_list.size());
    std::vector<size_t> elem_nodes;
    elem_nodes.reserve(mesh->element_list.size()*mesh->element_list[0]->nodes.size());
    size_t tag = 0;
    for(auto& e:mesh->element_list){
        types.push_back(e->get_gmsh_element_type());
        elem_tags.push_back(tag++);
        for(auto& n:e->nodes){
            elem_nodes.push_back(n->id);
        }
    }
}

void Visualization::update_view(){
    gmsh::view::remove(1);
    gmsh::view::add(this->STRESS_VIEW, 1);

    std::vector<std::vector<double>> stress;
    stress.reserve(mesh->node_list.size());
    std::vector<size_t> node_tags;
    node_tags.reserve(mesh->node_list.size());
    for(auto& n:mesh->node_list){
        node_tags.push_back(n->id);
        double* S = n->results;
        double s = std::sqrt(0.5*(S[0]*S[0] - S[0]*S[1] + S[1]*S[1] + 3*S[2]*S[2]));
        std::vector<double> tmp{s};
        stress.push_back(tmp);
    }
    gmsh::view::addModelData(1, 0, this->MODEL_NAME, "NodeData", node_tags, stress, 0, 1);

    if(this->shown){
        gmsh::graphics::draw();
    }
}

void Visualization::show(){
    this->shown = true;

    auto checkForEvent = [=]() -> bool {
        std::vector<std::string> action;
        gmsh::onelab::getString("ONELAB/Action", action);
        if(action.size() && action[0] == "check") {
            gmsh::onelab::setString("ONELAB/Action", {""});
            gmsh::graphics::draw();
        }
        return true;
    };

    gmsh::fltk::initialize();
    while(gmsh::fltk::isAvailable() && checkForEvent() && this->shown)
        gmsh::fltk::wait();
}
