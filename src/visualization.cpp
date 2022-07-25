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
        this->tag = gmsh::model::addDiscreteEntity(2);
        gmsh::model::mesh::addNodes(2, this->tag, node_tags, coords);
    } else if(type == utils::PROBLEM_TYPE_3D){
        this->tag = gmsh::model::addDiscreteEntity(3);
        gmsh::model::mesh::addNodes(3, this->tag, node_tags, coords);
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
    gmsh::model::mesh::addElementsByType(this->tag, mesh->elem_info->get_gmsh_element_type(), elem_tags, elem_nodes);

    // gmsh::option::setNumber("View.DrawLines", 0);
    // gmsh::option::setNumber("View.DrawPoints", 0);
    // gmsh::option::setNumber("View.DrawTrihedra", 0);
    // gmsh::option::setNumber("View.DrawTriangles", 0);
    gmsh::option::setNumber("Mesh.SurfaceEdges", 0);
    gmsh::option::setNumber("Mesh.VolumeEdges", 0);
    gmsh::option::setNumber("Mesh.Triangles", 0);
    gmsh::option::setNumber("General.ColorScheme", 3);
    gmsh::option::setNumber("General.FltkColorScheme", 1);
    gmsh::option::setString("General.GraphicsFontEngine", "StringTexture");
    

    gmsh::view::add(this->STRESS_VIEW, 1);
    gmsh::view::add("Normal Stress (X axis)", 2);
    gmsh::view::add("Normal Stress (Y axis)", 3);
    gmsh::view::add("Shear Stress", 4);
    gmsh::view::add(this->DENSITY_VIEW, 5);
    gmsh::view::option::setNumber(5, "ColormapNumber", 9); //grayscale
    gmsh::view::option::setNumber(5, "ColormapInvert", 1.0); //inverted
    gmsh::view::add("Displacement View", 6);

}

void Visualization::update_stress_view(const std::vector<double>& s, size_t id){
    logger::quick_log("Updating view...");

    std::vector<std::vector<double>> stress;
    size_t s_size = 0;
    for(auto& g:mesh->geometries){
        s_size += g->mesh.size();
    }
    stress.reserve(s_size);
    std::vector<size_t> elem_tags;
    elem_tags.reserve(s.size());
    for(size_t i = 0; i < s.size(); ++i){
        elem_tags.push_back(i+1);
        std::vector<double> tmp{s[i]};
        stress.push_back(tmp);
    }
    gmsh::view::addModelData(id, 0, this->MODEL_NAME, "ElementData", elem_tags, stress, 0, 1);

    if(this->shown){
        gmsh::graphics::draw();
    }
    logger::quick_log("Done.");
}

void Visualization::update_nodal_stress_view(const std::vector<double>& s){
    logger::quick_log("Updating view...");

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

    std::vector<std::vector<double>> density;
    density.reserve(mesh->node_list.size());
    std::vector<size_t> elem_tags;

    size_t s_size = 0;
    for(auto& g:mesh->geometries){
        s_size += g->mesh.size();
    }
    elem_tags.reserve(s_size);
    size_t i = 0;
    size_t di = 0;
    for(auto& g:mesh->geometries){
        if(g->do_topopt){
            for(auto& e:g->mesh){
                (void)e;
                elem_tags.push_back(i+1);
                std::vector<double> tmp{d[di]};
                density.push_back(tmp);
                ++i;
                ++di;
            }
        } else {
            for(auto& e:g->mesh){
                (void)e;
                elem_tags.push_back(i+1);
                std::vector<double> tmp{1.0};
                density.push_back(tmp);
                ++i;
            }
        }
    }
    gmsh::view::addModelData(5, 0, this->MODEL_NAME, "ElementData", elem_tags, density, 0, 1);

    if(this->shown){
        gmsh::graphics::draw();
    }
    logger::quick_log("Done.");
}

void Visualization::update_vector_view(const std::vector<std::unique_ptr<MeshNode>>& nodes, const std::vector<double>& values){
    logger::quick_log("Updating view...");

    std::vector<std::vector<double>> vecs;
    vecs.reserve(nodes.size());
    std::vector<size_t> node_tags;
    node_tags.reserve(nodes.size());
    if(this->type == utils::PROBLEM_TYPE_2D){
        for(size_t i = 0; i < nodes.size(); ++i){
            node_tags.push_back(i+1);
            std::vector<double> tmp;
            for(size_t j = 0; j < 2; ++j){
                if(nodes[i]->u_pos[j] > -1){
                    tmp.push_back(values[nodes[i]->u_pos[j]]);
                } else {
                    tmp.push_back(0);
                }
            }
            tmp.push_back(0);
            vecs.push_back(tmp);
        }
    }
    gmsh::view::addModelData(6, 0, this->MODEL_NAME, "NodeData", node_tags, vecs, 0, 3);

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
