/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
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

#include "view_handler.hpp"

ViewHandler::ViewHandler(const Meshing* const mesh, const std::string& model_name, const std::string& view_name, const ViewType view_type, const DataType data_type, utils::ProblemType problem_type, const size_t view_id)
    : model_name(model_name), view_type(view_type), data_type(data_type), view_id(view_id), mesh(mesh),
      elem_num(this->get_number_of_elements()),
      node_num(this->get_number_of_nodes()),
      mat_color_num(this->get_number_of_material_colors()),
      problem_type(problem_type){

    gmsh::view::add(view_name, view_id);
    if(this->mat_color_num == 2 && data_type == MATERIAL){
        gmsh::view::option::setNumber(view_id, "ColormapNumber", 9); //grayscale
        gmsh::view::option::setNumber(view_id, "ColormapInvert", 1.0); //inverted
    }
}


void ViewHandler::update_view(const std::vector<double>& data, const std::vector<size_t>& geometries) const{
    if(this->removed){
        return;
    }
    std::vector<size_t> tags;
    // All geometries
    if(geometries.size() == 0 || this->view_type != ELEMENTAL){ // Default
        tags.resize(this->elem_num);
        std::iota(tags.begin(), tags.end(), 1);
    } else if(this->view_type == ELEMENTAL){ // Elemental with select geometries
        size_t size_offset = 0;
        size_t tag_offset = 1;
        for(size_t i = 0; i < this->mesh->geometries.size(); ++i){
            const auto pos = std::find(geometries.begin(), geometries.end(), i);
            const auto size = this->mesh->geometries.size();
            if(pos < geometries.end()){
                tags.resize(tags.size() + size);
                std::iota(tags.begin()+size_offset, tags.begin()+size_offset+size, tag_offset);
                size_offset += size;
            }
            tag_offset += size;
        }
    } else { // Not elemental
        tags.resize(this->node_num);
        std::iota(tags.begin(), tags.end(), 1);
    }

    // Elemental view shouldn't really be used for displacement
    if(this->view_type == ELEMENTAL && (this->data_type == STRESS || this->data_type == OTHER)){
        this->update_elemental(data, tags);
    } else if(this->view_type == NODAL && (this->data_type == STRESS || this->data_type == OTHER)){
        this->update_nodal(data, tags);
    } else if(this->view_type == NODAL && (this->data_type == STRESS || this->data_type == OTHER)){
        this->update_tensor(data, tags);
    } else if(this->view_type == ELEMENTAL && this->data_type == MATERIAL){
        // TODO
    } else if(this->view_type == NODAL && this->data_type == MATERIAL){
        // TODO
    } else if(this->data_type == DISPLACEMENT){ // Can only be VECTOR
        std::vector<double> vecs;
        vecs.reserve(tags.size()*3);
        if(this->problem_type == utils::PROBLEM_TYPE_2D){
            for(const auto i:tags){
                const auto& node = this->mesh->node_list[i];
                for(size_t j = 0; j < 2; ++j){
                    if(node->u_pos[j] > -1){
                        vecs.push_back(data[node->u_pos[j]]);
                    } else {
                        vecs.push_back(0);
                    }
                }
                vecs.push_back(0);
            }
        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
            for(const auto i:tags){
                const auto& node = this->mesh->node_list[i];
                for(size_t j = 0; j < 3; ++j){
                    if(node->u_pos[j] > -1){
                        vecs.push_back(data[node->u_pos[j]]);
                    } else {
                        vecs.push_back(0);
                    }
                }
            }
        }
        this->update_vector(vecs, tags);
    } else if(this->view_type == VECTOR && this->data_type == OTHER){
        this->update_vector(data, tags);
    }
}

size_t ViewHandler::get_number_of_elements() const{
    size_t elem_num = 0;
    if(this->view_type == ELEMENTAL){
        for(const auto& g:this->mesh->geometries){
            elem_num += g->mesh.size();
        }
    }
    return elem_num;
}
size_t ViewHandler::get_number_of_nodes() const{
    size_t node_num = 0;
    if(this->view_type != ELEMENTAL){
        if(this->mesh->node_list.size() > 0){
            node_num = this->mesh->node_list.size();
        } else {
            for(const auto& g:this->mesh->geometries){
                node_num += g->node_list.size();
            }
        }
    }
    return node_num;
}
size_t ViewHandler::get_number_of_material_colors() const{
    size_t num = 0;
    for(const auto& g:this->mesh->geometries){
        if(g->do_topopt){
            if(g->with_void){
                ++num;
            }
            num += g->number_of_materials();
        } else {
            ++num;
        }
    }
    return num;
}

