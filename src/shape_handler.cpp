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
#include "general_solver/mumps_general.hpp"
#include "logger.hpp"
#include "math/matrix.hpp"
#include "project_data.hpp"
#include "utils.hpp"
#include <algorithm>
#include <set>

ShapeHandler::ShapeHandler(Meshing* mesh, std::vector<Geometry*> geometries, std::unique_ptr<shape_op::ShapeOp> root_op):
    mesh(mesh), geometries(std::move(geometries)),
    solver(std::make_unique<general_solver::MUMPSGeneral>()),
    root_op(std::move(root_op)){
};

void ShapeHandler::obtain_affected_nodes(){
    const size_t node_num = this->mesh->elem_info->get_nodes_per_element();
    const size_t bnode_num = this->mesh->elem_info->get_boundary_nodes_per_element();

    this->original_points.resize(this->mesh->node_list.size()*3, 0);
    this->shape_displacement.resize(this->original_points.size(), 0);

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

    auto short_list = this->apply_op(root_op.get());

    // TODO deduplicate superimposed nodes
    this->optimized_nodes.reserve(short_list.size());
    std::set<BoundaryElement*> elem_set;
    for(auto& ni:short_list){
        size_t af_id = 0;
        while(af_id < affected.size()){
            if(this->mesh->boundary_node_list[af_id]->id == ni->id){
                break;
            }
            ++af_id;
        }
        if(affected[af_id]){
            const size_t curr_id = ni->id;
            std::vector<size_t> node_ids{curr_id};
            std::vector<AffectedElement> elems;
            const auto& range = this->mesh->inverse_mesh.equal_range(curr_id);
            size_t e_size = 0;
            for(auto it = range.first; it != range.second; ++it){
                ++e_size;
            }
            elems.reserve(e_size);
            for(auto it = range.first; it != range.second; ++it){
                auto e = it->second;
                size_t n;
                for(n = 0; n < node_num; ++n){
                    if(e->nodes[n]->id == curr_id){
                        break;
                    }
                }
                elems.push_back(AffectedElement{e, n});
                auto& vec = this->elem_to_affected_node_mapping[e];
                vec.reserve(bnode_num);
                vec.push_back(n);
            }
            const auto& brange = this->mesh->boundary_inverse_mesh.equal_range(curr_id);
            for(auto it = brange.first; it != brange.second; ++it){
                auto e = it->second;
                elem_set.insert(e);
            }
            this->optimized_nodes.push_back(AffectedNode{std::move(node_ids), std::move(elems)});
        }
    }
    for(size_t i = 0; i < this->optimized_nodes.size(); ++i){
        for(const auto id:this->optimized_nodes[i].node_ids){
            this->optimized_nodes_mapping[id] = i;
        }
    }
    this->boundary_elements.reserve(elem_set.size());
    this->boundary_elements.insert(this->boundary_elements.begin(), elem_set.begin(), elem_set.end());

    const auto bcomp = [](const BoundaryElement* e1, const BoundaryElement* e2) -> bool{
        return e1->id < e2->id;
    };
    std::sort(this->boundary_elements.begin(), this->boundary_elements.end(), bcomp);

    // Generate shape elements
    const auto elem_maker = this->mesh->elem_info->get_shape_element_info();
    this->shape_elements.reserve(this->boundary_elements.size());
    std::vector<MeshNode*> nodes(bnode_num);
    for(size_t i = 0; i < bnode_num; ++i){
        nodes[i] = this->mesh->node_list[this->boundary_elements[0]->nodes[i]->id].get();
    }
    this->bound_to_shape_mapping[0] = 0;
    ElementShape es{std::move(nodes), this->boundary_elements[0]->normal};
    this->shape_elements.emplace_back(elem_maker->make_element(std::move(es)));

    size_t shape_id = 1;
    for(size_t j = 1; j < this->boundary_elements.size(); ++j){
        std::vector<MeshNode*> nodes(bnode_num);
        if(this->boundary_elements[j]->get_centroid(bnode_num).IsEqual(
                    this->boundary_elements[j-1]->get_centroid(bnode_num),
                    Precision::Confusion())){

            this->bound_to_shape_mapping[j] = shape_id - 1;
            continue;
        }
        for(size_t i = 0; i < bnode_num; ++i){
            nodes[i] = this->mesh->node_list[this->boundary_elements[j]->nodes[i]->id].get();
        }
        this->bound_to_shape_mapping[j] = shape_id;
        ElementShape es{std::move(nodes), this->boundary_elements[j]->normal};
        this->shape_elements.emplace_back(elem_maker->make_element(std::move(es)));
        ++shape_id;
    }


    // Save original coordinates (for visualization)
    for(size_t i = 0; i < this->mesh->node_list.size(); ++i){
        const auto& n = this->mesh->node_list[i];
        for(size_t j = 0; j < 3; ++j){
            this->original_points[3*i + j] = n->point.Coord(1+j);
        }
    }

    const auto& g = this->geometries[0];
    // Initialize id_mapping
    if(!this->full_boundary_optimization){
        for(const auto& n:g->node_list){
            this->id_mapping[n->id] = -1;
        }
    }
    // Generate linear problem
    // TODO: expand for multiple geometries
    const auto ecomp = [](const Node* e1, const Node* e2) -> bool{
        return e1->id < e2->id;
    };
    std::set<const Node*, std::function<bool(const Node*, const Node*)>> domain_nodes(ecomp);
    std::set<const Node*, std::function<bool(const Node*, const Node*)>> bound_nodes(ecomp);
    for(const auto& n:this->optimized_nodes){
        bool found = false;
        for(const auto& e:g->boundary_node_list){
            for(size_t nn:n.node_ids){
                if(nn == e->id){
                    bound_nodes.insert(e);
                    found = true;
                    break;
                }
            }
            if(found) break;
        }
    }
    domain_nodes.insert(g->node_list.begin(), g->node_list.end());
    for(const auto& b:g->boundary_node_list){
        auto it = domain_nodes.find(static_cast<const Node*>(b));
        auto keep = bound_nodes.find(static_cast<const Node*>(b));
        if(keep == bound_nodes.end() && it != domain_nodes.end()){
            domain_nodes.erase(it);
        }
    }

    // Create node to element inverse mapping to calculate nodal gradients 
    // (that is, displacement vectors)
    for(const auto& b:g->boundary_node_list){
        auto it = domain_nodes.find(static_cast<const Node*>(b));
        if(it != domain_nodes.end()){
            domain_nodes.erase(it);
        }
    }
    this->domain_nodes.reserve(domain_nodes.size());
    for(auto n:domain_nodes){
        this->domain_nodes.push_back(this->mesh->node_list[n->id].get());
    }

    if(this->geometries.size() > 1){
        return;
    }

    size_t idm = 0;
    // Optimized nodes have non-homogeneous Dirichlet.
    // Non-optimized nodes have homogeneous Dirichlet.
    // As a result, the displacement spreading only affects inner nodes, which
    // means that it doesn't modify the shape aside from the optimized region.
    //
    // I'm not sure how good of an idea this is, but it's a start.
    for(const auto& n:g->node_list){
        this->id_mapping[n->id] = -1;
    }
    for(const auto& n:domain_nodes){
        this->id_mapping[n->id] = idm;
        ++idm;
    }
    this->matrix_width = idm;
    this->solver->initialize_matrix(true, idm);
    this->b.resize(idm);

    const math::Matrix A({1, 0, 0,
                          0, 1, 0,
                          0, 0, 1}, 3, 3);

    std::vector<long> pos(node_num);
    if(this->full_boundary_optimization){
        for(const auto& e:g->mesh){
            for(size_t i = 0; i < node_num; ++i){
                pos[i] = this->id_mapping[e->nodes[i]->id];
            }
            const auto k = e->diffusion_1dof(this->mesh->thickness, A);// + 1e-3*e->absorption_1dof(this->mesh->thickness);
            this->solver->add_element(k, pos);
        }
    } else {
        for(const auto& e:g->mesh){
            for(size_t i = 0; i < node_num; ++i){
                pos[i] = this->id_mapping[e->nodes[i]->id];
            }
            const auto k = e->diffusion_1dof(this->mesh->thickness, A);
            this->solver->add_element(k, pos);
        }
    }
    this->solver->compute();
}
    
void ShapeHandler::update_nodes(const std::vector<double>& dx){
    const size_t dof = this->mesh->elem_info->get_dof_per_node();
    const size_t bnum = this->mesh->elem_info->get_boundary_nodes_per_element();
    const size_t num = this->mesh->elem_info->get_nodes_per_element();
    const size_t dim =
        (this->mesh->elem_info->get_problem_type() == utils::ProblemType::PROBLEM_TYPE_2D)
        ? 2 : 3;

    // Update boundary coordinates
    for(size_t i = 0; i < this->optimized_nodes.size(); ++i){
        auto& nids = this->optimized_nodes[i].node_ids;
        for(auto& nid:nids){
            auto& n = this->mesh->node_list[nid];
            for(size_t d = 0; d < dof; ++d){
                n->point.SetCoord(1+d, n->point.Coord(1+d) + dx[i*dof + d]);
            }
        }
    }
    for(auto& b:this->mesh->boundary_elements){
        b.update_normal(bnum, this->mesh->proj_data->type);
    }

    //const double mult = (this->full_boundary_optimization) ? 1e-3 : 0;
    const double mult = 0;

    const math::Matrix A({1, 0, 0,
                          0, 1, 0,
                          0, 0, 1}, 3, 3);

    if(this->geometries.size() == 1){

    // Update other nodes
    math::Vector bn_vals(num);
    for(size_t dim_i = 0; dim_i < dim; ++dim_i){
        std::fill(this->b.begin(), this->b.end(), 0);
        for(auto& it: this->elem_to_affected_node_mapping){
            const auto e = it.first;
            const auto& nodes = it.second;
            for(const auto ni: nodes){
                const auto opt_id = this->optimized_nodes_mapping.find(e->nodes[ni]->id);
                if(opt_id != this->optimized_nodes_mapping.end()){
                    bn_vals[ni] = dx[dim*(opt_id->second) + dim_i];
                }
            }
            const math::Vector Fe = (e->diffusion_1dof(this->mesh->thickness, A) + mult*e->absorption_1dof(this->mesh->thickness))*bn_vals;
            for(size_t j = 0; j < num; ++j){
                const long global_id = this->id_mapping[e->nodes[j]->id];
                if(global_id > -1){
                    this->b[global_id] -= Fe[j];
                }
            }
            bn_vals.fill(0);
        }
        this->solver->solve(b);
        for(auto ni:this->domain_nodes){
            const long global_id = this->id_mapping[ni->id];
            if(global_id > -1){
                ni->point.SetCoord(1+dim_i, ni->point.Coord(1+dim_i) + this->b[global_id]);
            }
        }
    }

    }
    
    // Update view
    for(size_t i = 0; i < this->mesh->node_list.size(); ++i){
        const auto& n = this->mesh->node_list[i];
        for(size_t j = 0; j < 3; ++j){
            this->shape_displacement[3*i + j] = n->point.Coord(1+j) - this->original_points[3*i + j];
        }
    }

    // Update all coefficients
    for(auto& g:this->mesh->geometries){
        #pragma omp parallel for
        for(auto& e:g->mesh){
            e->calculate_coefficients();
        }
    }
}

std::set<MeshNode*> ShapeHandler::apply_op(shape_op::ShapeOp* op) const{
    switch(op->get_type()){
        case shape_op::Code::GEOMETRY:{
            logger::log_assert(op->get_id() < this->geometries.size(), logger::ERROR,
                    "unknown geometry id: {}", op->get_id());
            const auto g = this->geometries[op->get_id()];
            std::set<MeshNode*> nodes(g->boundary_node_list.begin(), g->boundary_node_list.end());

            return nodes;
        }
        case shape_op::Code::UNION:{
            auto s1 = this->apply_op(op->first());
            auto s2 = this->apply_op(op->second());
            s1.insert(s2.begin(), s2.end());
            s2.clear();
            return s1;
        }
        case shape_op::Code::INTERSECTION:{
            auto s1 = this->apply_op(op->first());
            auto s2 = this->apply_op(op->second());
            std::vector<MeshNode*> result_vec(std::min(s1.size(), s2.size()));
            std::set<MeshNode*> result;
            auto input_end = std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), result_vec.begin());
            s1.clear();
            s2.clear();
            result.insert(result_vec.begin(), input_end);
            result_vec.clear();
            return result;
        }
        case shape_op::Code::DIFFERENCE:{
            auto s1 = this->apply_op(op->first());
            auto s2 = this->apply_op(op->second());
            std::vector<MeshNode*> result_vec(s1.size());
            std::set<MeshNode*> result;
            auto input_end = std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(), result_vec.begin());
            s1.clear();
            s2.clear();
            result.insert(result_vec.begin(), input_end);
            result_vec.clear();
            return result;
        }
        case shape_op::Code::SHELL:{
            logger::log_assert(false, logger::ERROR,
                    "SHELL shape operation is currently not implemented");
            return std::set<MeshNode*>();
        }
    }

    return std::set<MeshNode*>();
}
