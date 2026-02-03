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

#include <defs.hpp>
#include <limits>
#include <set>
#include "view_handler.hpp"
#include "logger.hpp"

ViewHandler::ViewHandler(const Meshing* const mesh, spview::Server* server, const std::string& view_name, const spview::defs::ViewType view_type, const spview::defs::DataType data_type, utils::ProblemType problem_type, const size_t view_id)
    : view_type((data_type != spview::defs::MESH) ? view_type : spview::defs::ELEMENTAL),
      data_type(data_type), view_id(view_id), mesh(mesh),
      elem_num(this->get_number_of_elements()),
      node_num(this->get_number_of_nodes()),
      mat_color_num(this->get_number_of_material_colors()),
      problem_type(problem_type),
      server(server){

    this->server->add_view(view_type, data_type, view_name);
}


void ViewHandler::update_view(const std::vector<double>& data, const std::vector<size_t>& geometries) const{
    if(this->removed){
        return;
    }
    std::vector<size_t> tags;
    // All geometries
    if(geometries.size() == 0 && this->view_type == spview::defs::ELEMENTAL){ // Default
        tags.resize(this->elem_num);
        std::iota(tags.begin(), tags.end(), 0);
    } else if(this->view_type == spview::defs::ELEMENTAL){ // Elemental with select geometries
        size_t tag_size = 0;
        for(auto& gi:geometries){
            auto& g = this->mesh->geometries[gi];
            tag_size += g->mesh.size();
        }
        tags.reserve(tag_size);
        for(auto& gi:geometries){
            auto& g = this->mesh->geometries[gi];
            for(auto& e:g->mesh){
                tags.push_back(e->id);
            }
        }
    } else if(geometries.size() == 0){ // Not elemental
        tags.resize(this->node_num);
        std::iota(tags.begin(), tags.end(), 0);
    } else {
        std::set<size_t> node_ids;
        for(size_t i = 0; i < this->mesh->geometries.size(); ++i){
            const auto pos = std::find(geometries.begin(), geometries.end(), i);
            if(pos < geometries.end()){
                const auto& g = this->mesh->geometries[i];
                for(const auto& n:g->node_list){
                    node_ids.insert(n->id);
                }
            }
        }

        size_t num_nodes = node_ids.size();

        tags.resize(num_nodes);
        auto tag_it = tags.begin();
        for(auto& n:node_ids){
            *tag_it = n;
            ++tag_it;
        }
    }

    // Elemental view shouldn't really be used for displacement
    if(this->view_type == spview::defs::ELEMENTAL && (this->data_type == spview::defs::STRESS || this->data_type == spview::defs::OTHER)){
        this->server->update_data(this->view_id, tags, data);
    } else if(this->view_type == spview::defs::NODAL && (this->data_type == spview::defs::STRESS || this->data_type == spview::defs::OTHER)){
        this->server->update_data(this->view_id, tags, data);
    } else if(this->view_type == spview::defs::TENSOR && (this->data_type == spview::defs::STRESS || this->data_type == spview::defs::OTHER)){
        this->server->update_data(this->view_id, tags, data);
    } else if(this->data_type == spview::defs::DENSITY){
        if(this->mesh->geometries.size() == 1){
            if(view_type == spview::defs::ELEMENTAL){
                this->server->update_data(this->view_id, tags, data);
            } else if(view_type == spview::defs::NODAL){
                this->server->update_data(this->view_id, tags, data);
            }
        } else {
            if(view_type == spview::defs::ELEMENTAL){
                std::vector<double> mat(tags.size(), std::numeric_limits<double>::quiet_NaN());
                if(geometries.size() == 0){
                    auto m = mat.begin();
                    auto d = data.cbegin();
                    for(const auto& g:this->mesh->geometries){
                        if(g->do_topopt){
                            for(size_t i = 0; i < g->mesh.size(); ++i){
                                *m = *d;
                                ++m;
                                ++d;
                            }
                        }
                    }
                    this->server->update_data(this->view_id, tags, mat);
                } else {
                    this->server->update_data(this->view_id, tags, data);
                }
            } else if(view_type == spview::defs::NODAL){
                this->server->update_data(this->view_id, tags, data);
            }
        }
    } else if(this->view_type == spview::defs::ELEMENTAL && this->data_type == spview::defs::MATERIAL){
        if(this->mat_color_num == 2){
            this->server->update_data(this->view_id, tags, data);
        } else {
            std::vector<double> mat(tags.size());
            double cur_mat = 0.0;
            size_t mat_pos = 0;
            if(geometries.size() == 0){
                for(const auto& g:this->mesh->geometries){
                    if(g->number_of_densities_needed() == 1){
                        for(size_t i = mat_pos; i < mat_pos+g->mesh.size(); ++i){
                            mat[i] = cur_mat + data[i];
                        }
                        mat_pos += g->mesh.size();
                        cur_mat += 2.0;
                    } else {
                        // TODO
                        // Probably needs an overhaul of visualizaiton method,
                        // e.g. a custom viewer.
                    }
                }
            } else {
                // TODO
            }
            this->server->update_data(this->view_id, tags, mat);
        }
    } else if(this->data_type == spview::defs::MESH){
        std::vector<double> mat(tags.size());
        size_t mat_pos = 0;
        size_t geom_id = 0;
        for(const auto& g:this->mesh->geometries){
            std::fill(mat.begin()+mat_pos, mat.begin()+mat_pos+g->mesh.size(), geom_id);
            mat_pos += g->mesh.size();
            ++geom_id;
        }
        this->server->update_data(this->view_id, tags, mat);
    } else if(this->view_type == spview::defs::NODAL && this->data_type == spview::defs::MATERIAL){
        if(this->mat_color_num == 2){
            this->server->update_data(this->view_id, tags, data);
        } else {
            std::vector<double> mat(tags.size());
            double cur_mat = 0.0;
            size_t mat_pos = 0;
            if(geometries.size() == 0){
                for(const auto& g:this->mesh->geometries){
                    if(g->number_of_densities_needed() == 1){
                        for(size_t i = mat_pos; i < mat_pos+g->mesh.size(); ++i){
                            mat[i] = cur_mat + data[i];
                        }
                        mat_pos += g->mesh.size();
                        cur_mat += 2.0;
                    } else {
                        // TODO
                        // Probably needs an overhaul of visualizaiton method,
                        // e.g. a custom viewer.
                    }
                }
            } else {
                // TODO
            }
            this->server->update_data(this->view_id, tags, mat);
        }
    } else if(this->data_type == spview::defs::DISPLACEMENT){ // Can only be VECTOR
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
        this->server->update_data(this->view_id, tags, vecs);
    } else if(this->view_type == spview::defs::VECTOR && this->data_type == spview::defs::OTHER){
        this->server->update_data(this->view_id, tags, data);
    }
}

size_t ViewHandler::get_number_of_elements() const{
    size_t elem_num = 0;
    if(this->view_type == spview::defs::ELEMENTAL){
        for(const auto& g:this->mesh->geometries){
            elem_num += g->mesh.size();
        }
    }
    return elem_num;
}
size_t ViewHandler::get_number_of_nodes() const{
    size_t node_num = 0;
    if(this->view_type != spview::defs::ELEMENTAL){
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

