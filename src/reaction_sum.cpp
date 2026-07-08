/*
 *   Copyright (C) 2026 Tarcísio Ladeia de Oliveira.
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

#include "reaction_sum.hpp"
#include "logger.hpp"
#include "meshing.hpp"
#include "utils.hpp"

ReactionSum::ReactionSum(CrossSection cross_section, double thickness, std::array<double, 3> F, utils::ProblemType type):
    S(std::move(cross_section)), thickness(thickness), F(std::move(F)),
    F_dim((type == utils::ProblemType::PROBLEM_TYPE_2D) ? 2 : 3),
    type(type)
{
}

void ReactionSum::initialize(const Meshing* const mesh, size_t& max_dofs){
    this->elem_info = mesh->elem_info;

    this->offset = max_dofs;
    max_dofs += this->F_dim;

    size_t bound_nodes_per_elem = mesh->elem_info->get_boundary_nodes_per_element();
    this->BOUNDARY_NODE_NUM = elem_info->get_boundary_nodes_per_element();
    const size_t node_num = elem_info->get_nodes_per_element();

    const auto& boundary_elements = mesh->boundary_elements;
    if(this->submesh.size() == 0){
        std::vector<bool> apply_spring(boundary_elements.size());
        size_t num_elems = 0;
        // Fill apply_spring concurrently
        // (CrossSection::is_inside() is pretty slow, unfortunately
        #pragma omp parallel for reduction(+:num_elems)
        for(size_t i = 0; i < boundary_elements.size(); ++i){
            const auto& e = boundary_elements[i];
            apply_spring[i] = this->S.is_inside(e.get_centroid(bound_nodes_per_elem));
            if(apply_spring[i]){
                ++num_elems;
            }
        }

        this->submesh.resize(num_elems);
        size_t cur_elem = 0;
        for(size_t i = 0; i < boundary_elements.size(); ++i){
            if(apply_spring[i]){
                this->submesh[cur_elem] = &boundary_elements[i];
                for(size_t n = 0; n < node_num; ++n){
                    this->affected_nodes[boundary_elements[i].parent->nodes[n]->id] = false;
                }
                for(size_t n = 0; n < BOUNDARY_NODE_NUM; ++n){
                    this->affected_nodes[boundary_elements[i].nodes[n]->id] = true;
                }

                ++cur_elem;
            }
        }
    }
}

void ReactionSum::apply_constraints(std::vector<double>& load) const{
    for(size_t i = 0; i < F_dim; ++i){
        load[this->offset + i] = this->F[i];
    }
}

math::Matrix ReactionSum::get_constr(const Meshing* const mesh, const BoundaryElement* b) const{

    const auto Ke = b->parent->get_k(mesh->geometries[b->geom_id]->materials, this->thickness);

    const auto KW = Ke.get_W();

    math::Matrix result(this->F_dim, Ke.get_W(), 0);
    
    const size_t node_num = this->elem_info->get_nodes_per_element();
    const size_t dof = this->elem_info->get_dof_per_node();
    for(size_t i = 0; i < node_num; ++i){
        if(this->affected_nodes.at(b->parent->nodes[i]->id)){
            for(size_t f = 0; f < this->F_dim; ++f){
                for(size_t j = 0; j < KW; ++j){
                    result(f, j) += Ke(i*dof + f, j);
                }
            }
        }
    }

    return result;
}
