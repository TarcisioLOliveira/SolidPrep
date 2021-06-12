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

#include "meshing.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include "project_data.hpp"
#include <algorithm>


void Meshing::prepare_for_FEM(const std::vector<ElementShape>& base_mesh,
                              MeshElementFactory::MeshElementType element_type,
                              ProjectData* data,
                              std::vector<MeshElement*>& final_mesh,
                              std::vector<double>& load_vector){
    logger::log_assert(this->node_list.size() > 0, logger::ERROR, "the object's node list is empty. Ensure you are using the same Meshing instance as the one used to obtain the list of ElementShape instances");
    auto test = base_mesh[0].nodes[0];
    bool correct = false;
    for(auto& n : this->node_list){
        if(n.get() == test){
            correct = true;
            break;
        }
    }
    logger::log_assert(correct, logger::ERROR, "object mismatch. Please ensure that the Meshing instance used to generate the list of ElementShape instances is the same as the one being used to prepare the mesh for FEM.");

    this->element_list.clear();
    this->element_list.reserve(base_mesh.size());
    final_mesh.reserve(base_mesh.size());

    auto comp = [](const std::unique_ptr<MeshNode>& a, const std::unique_ptr<MeshNode>& b){ return a->id < b->id;};
    std::sort(this->node_list.begin(), this->node_list.end(), comp);

    size_t dof = MeshElementFactory::get_dof_per_node(element_type);
    size_t k_size = MeshElementFactory::get_k_dimension(element_type);

    std::vector<long> world_pos;
    world_pos.resize(this->node_list.size()*dof);

    size_t current = 0;
    for(size_t i = 0; i < this->node_list.size(); ++i){
        auto& n = this->node_list[i];
        bool supported = false;
        size_t max_offset = 0;
        for(auto& s : data->supports){
            if(s.S.is_inside(n->point)){
                size_t offset = 0;
                std::vector<long> sup_pos = this->get_support_dof(offset, 0, s, element_type);
                for(size_t j = 0; j < dof; ++j){
                    if(sup_pos[j] > 0){
                        if(world_pos[i*dof + j] >= 0){
                            world_pos[i*dof + j] = sup_pos[j] + current;
                        }
                    } else {
                        world_pos[i*dof + j] = sup_pos[j];
                    }
                }
                max_offset = std::max(offset, max_offset);
                supported = true;
            }
        }
        if(!supported){
            for(size_t j = 0; j < dof; ++j){
                world_pos[i*dof + j] = j + current;
            }
            current += dof;
        } else {
            current += max_offset;
        }
    }

    load_vector.resize(current+1);

    current = 0;
    for(size_t i = 0; i < this->node_list.size(); ++i){
        auto& n = this->node_list[i];
        for(auto& f : data->forces){
            if(f.S.is_inside(n->point)){
                std::vector<double> f_vec = this->get_force_dof(f, element_type);
                for(size_t j = 0; j < dof; ++j){
                    if(world_pos[current] >= 0){
                        load_vector[world_pos[current]] += f_vec[j];
                        ++current;
                    }
                }
            }
        }
    }


    for(auto& e : base_mesh){
        std::vector<long> u_pos(k_size);
        for(size_t i = 0; i < e.nodes.size(); ++i){
            auto& n = e.nodes[i];
            for(size_t j = 0; j < dof; ++j){
                u_pos[i*dof + j] = world_pos[n->id*dof + j];
            }
        }

        this->element_list.emplace_back(MeshElementFactory::make_element(element_type, u_pos, e, data));
        final_mesh.push_back(this->element_list.back().get());
    }
}

std::vector<long> Meshing::get_support_dof(size_t& offset, size_t id, const Support& support, MeshElementFactory::MeshElementType type) const{
    size_t size = MeshElementFactory::get_dof_per_node(type);
    utils::ProblemType prob_type = MeshElementFactory::get_problem_type(type);
    std::vector<long> pos(size);
    id *= size;
    switch(size){
        case 6:
            pos[5] = support.MZ ? -1 : (id + offset++);
            pos[4] = support.MY ? -1 : (id + offset++);
            pos[3] = support.MX ? -1 : (id + offset++);
        case 3:
            if(prob_type == utils::PROBLEM_TYPE_2D){
                pos[2] = support.MZ ? -1 : (id + offset++);
            } else {
                pos[2] = support.Z ? -1 : (id + offset++);
            }
        case 2:
            pos[1] = support.Y ? -1 : (id + offset++);
            pos[0] = support.X ? -1 : (id + offset++);
    }

    return pos;
}

std::vector<double> Meshing::get_force_dof(const Force& force, MeshElementFactory::MeshElementType type) const{
    size_t size = MeshElementFactory::get_dof_per_node(type);
    utils::ProblemType prob_type = MeshElementFactory::get_problem_type(type);
    std::vector<double> f(size);
    switch(size){
        case 6:
        case 3:
            if(prob_type == utils::PROBLEM_TYPE_3D){
                f[2] = -force.vec.Z();
            }
        case 2:
            f[1] = -force.vec.Y();
            f[0] = -force.vec.X();
    }



    return f;
}
