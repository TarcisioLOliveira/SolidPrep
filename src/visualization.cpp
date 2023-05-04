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
#include <algorithm>
#include <vector>
#include "logger.hpp"

Visualization::Visualization():
    server(this->server_name()){

}

void Visualization::load_mesh(Meshing* mesh, utils::ProblemType type){
    this->mesh = mesh;
    this->type = type;

    std::vector<double> coords;
    coords.reserve(mesh->node_list.size()*3);
    for(auto& n:mesh->node_list){
        coords.push_back(n->point.X());
        coords.push_back(n->point.Y());
        coords.push_back(n->point.Z());
    }

    size_t s_size = 0;
    for(auto& g:mesh->geometries){
        s_size += g->mesh.size();
    }
    std::vector<size_t> elem_nodes;
    const size_t nodes_num = mesh->elem_info->get_nodes_per_element();
    elem_nodes.reserve(s_size*nodes_num);
    for(auto& g:mesh->geometries){
        for(auto& e:g->mesh){
            for(size_t i = 0; i < nodes_num; ++i){
                const auto& n = e->nodes[i];
                elem_nodes.push_back(n->id);
            }
        }
    }

    spview::defs::InitData init;
    if(type == utils::PROBLEM_TYPE_2D){
        init.model_type = spview::defs::MODEL_2D;
    } else if(type == utils::PROBLEM_TYPE_3D){
        init.model_type = spview::defs::MODEL_3D;
    }
    init.backend = spview::defs::GMSH;
    init.element = mesh->elem_info->get_spview_code();
    init.elem_num = s_size;
    init.node_num = mesh->node_list.size();
    init.mat_num = this->get_number_of_material_colors();

    this->server.init_client(init, coords, elem_nodes);

}

ViewHandler* Visualization::add_view(const std::string& view_name, spview::defs::ViewType view_type, spview::defs::DataType data_type){
    this->handler_list.emplace_back(std::make_unique<ViewHandler>(this->mesh, &this->server, view_name, view_type, data_type, this->type, this->last_view_tag));
    ++this->last_view_tag;
    return this->handler_list.back().get();
}

size_t Visualization::get_number_of_material_colors() const{
    size_t num = 0;
    for(const auto& g:this->mesh->geometries){
        if(g->do_topopt){
            if(g->with_void || g->number_of_materials() == 1){
                ++num;
            }
            num += g->number_of_materials();
        } else {
            ++num;
        }
    }
    return num;
}
