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
#include "utils.hpp"
#include <gmsh.h>
#include <algorithm>
#include <vector>
#include "logger.hpp"

void Visualization::load_mesh(Meshing* mesh, utils::ProblemType type){
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

        // It looks like node tags have to begin at 1 for some reason.
        node_tags.push_back(n->id+1);
    }

    if(type == utils::PROBLEM_TYPE_2D){
        this->tag = gmsh::model::addDiscreteEntity(2);
        gmsh::model::mesh::addNodes(2, this->tag, node_tags, coords);
    } else if(type == utils::PROBLEM_TYPE_3D){
        this->tag = gmsh::model::addDiscreteEntity(3);
        gmsh::model::mesh::addNodes(3, this->tag, node_tags, coords);
    }

    std::vector<size_t> elem_tags;
    elem_tags.reserve(mesh->element_list.size());
    std::vector<size_t> elem_nodes;
    elem_nodes.reserve(mesh->element_list.size()*mesh->element_list[0]->nodes.size());
    size_t etag = 1;
    for(auto& e:mesh->element_list){
        elem_tags.push_back(etag++);
        for(auto& n:e->nodes){
            elem_nodes.push_back(n->id+1);
        }
    }
    gmsh::model::mesh::addElementsByType(this->tag, mesh->element_list[0]->get_gmsh_element_type(), elem_tags, elem_nodes);
}

void Visualization::update_stress_view(const std::vector<double>& s){
    logger::quick_log("Updating view...");
    gmsh::view::add(this->STRESS_VIEW, 1);

    std::vector<std::vector<double>> stress;
    stress.reserve(mesh->element_list.size());
    std::vector<size_t> elem_tags;
    elem_tags.reserve(s.size());
    for(size_t i = 0; i < s.size(); ++i){
        elem_tags.push_back(i+1);
        std::vector<double> tmp{s[i]};
        stress.push_back(tmp);
    }
    gmsh::view::addModelData(1, 0, this->MODEL_NAME, "ElementData", elem_tags, stress, 0, 1);

    if(this->shown){
        gmsh::graphics::draw();
    }
    logger::quick_log("Done.");
}

void Visualization::update_nodal_stress_view(const std::vector<double>& s){
    logger::quick_log("Updating view...");
    gmsh::view::add(this->STRESS_VIEW, 1);

    std::vector<std::vector<double>> stress;
    stress.reserve(mesh->node_list.size());
    std::vector<size_t> node_tags;
    node_tags.reserve(s.size());
    for(size_t i = 0; i < s.size(); ++i){
        node_tags.push_back(i+1);
        std::vector<double> tmp{s[i]};
        stress.push_back(tmp);
    }
    gmsh::view::addModelData(1, 0, this->MODEL_NAME, "NodeData", node_tags, stress, 0, 1);

    if(this->shown){
        gmsh::graphics::draw();
    }
    logger::quick_log("Done.");
}

void Visualization::update_density_view(const std::vector<double>& d){
    logger::quick_log("Updating view...");
    gmsh::view::add(this->DENSITY_VIEW, 2);

    std::vector<std::vector<double>> density;
    density.reserve(mesh->node_list.size());
    std::vector<size_t> elem_tags;
    elem_tags.reserve(mesh->element_list.size());
    for(size_t i = 0; i < mesh->element_list.size(); ++i){
        elem_tags.push_back(i+1);
        std::vector<double> tmp{d[i]};
        density.push_back(tmp);
    }
    gmsh::view::addModelData(2, 0, this->MODEL_NAME, "ElementData", elem_tags, density, 0, 1);

    if(this->shown){
        gmsh::graphics::draw();
    }
    logger::quick_log("Done.");
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
