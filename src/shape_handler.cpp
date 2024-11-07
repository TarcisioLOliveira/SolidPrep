/*
 *   Copyright (C) 2024 Tarcísio Ladeia de Oliveira.
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

#include "shape_handler.hpp"
#include "project_data.hpp"

ShapeHandler::ShapeHandler(Meshing* mesh, std::vector<Geometry*> geometries):
    mesh(mesh), geometries(std::move(geometries)){

};

void ShapeHandler::obtain_affected_nodes(){
    const size_t node_num = this->mesh->elem_info->get_nodes_per_element();

    std::vector<bool> affected(this->mesh->boundary_node_list.size(), true);
    size_t affected_num = 0;
    // This is quite slow, so it's better to parallelize it
    #pragma omp parallel for reduction(+:affected_num)
    for(size_t i = 0; i < this->mesh->boundary_node_list.size(); ++i){
        const auto& b = this->mesh->boundary_node_list[i];
        for(const auto& f:this->mesh->proj_data->forces){
            if(f.S.is_inside(b->point)){
                affected[i] = false;
                break;
            }
        }
        if(affected[i]){
            for(const auto& f:this->mesh->proj_data->supports){
                if(f.S.is_inside(b->point)){
                    affected[i] = false;
                    break;
                }
            }
        }
        if(affected[i]){
            for(const auto& f:this->mesh->proj_data->springs){
                if(f.S.is_inside(b->point)){
                    affected[i] = false;
                    break;
                }
            }
        }
        if(affected[i]){
            for(const auto& f:this->mesh->proj_data->internal_loads){
                if(f.S.is_inside(b->point)){
                    affected[i] = false;
                    break;
                }
            }
        }
        if(affected[i]){
            ++affected_num;
        }
    }

    // TODO deduplicate superimposed nodes
    this->optimized_nodes.reserve(affected_num);
    for(size_t i = 0; i < affected.size(); ++i){
        if(affected[i]){
            std::vector<size_t> node_ids{i};
            std::vector<AffectedElement> elems;
            const auto& range = this->mesh->inverse_mesh.equal_range(i);
            size_t e_size = 0;
            for(auto it = range.first; it != range.second; ++it){
                ++e_size;
            }
            elems.reserve(e_size);
            for(auto it = range.first; it != range.second; ++it){
                auto e = it->second;
                size_t n;
                for(n = 0; n < node_num; ++n){
                    if(e->nodes[n]->id == i){
                        break;
                    }
                }
                elems.push_back(AffectedElement{e, n});
            }
            this->optimized_nodes.push_back(AffectedNode{std::move(node_ids), std::move(elems)});
        }
    }
}
    
void ShapeHandler::update_nodes(const std::vector<double>& dx){
    const size_t dof = this->mesh->elem_info->get_dof_per_node();
    const size_t bnum = this->mesh->elem_info->get_boundary_nodes_per_element();

    for(size_t i = 0; i < this->optimized_nodes.size(); ++i){
        auto& nids = this->optimized_nodes[i].node_ids;
        for(auto& nid:nids){
            auto& n = this->mesh->node_list[nid];
            for(size_t d = 0; d < dof; ++d){
                n->point.SetCoord(1+d, n->point.Coord(1+d) + dx[i*dof + d]);
            }
        }
    }
    for(auto& e:this->affected_elements){
        e->calculate_coefficients();
    }
    for(auto& b:this->mesh->boundary_elements){
        b.update_normal(bnum, this->mesh->proj_data->type);
    }
}
