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
#include "general_solver/petsc_general_pcg.hpp"
#include "logger.hpp"
#include "math/matrix.hpp"
#include "project_data.hpp"
#include "utils.hpp"
#include <algorithm>
#include <set>

ShapeHandler::ShapeHandler(Meshing* mesh, std::vector<Geometry*> geometries, std::unique_ptr<shape_op::ShapeOp> root_op, double smoothing_radius):
    mesh(mesh), geometries(std::move(geometries)),
    root_op(std::move(root_op)),
    r(smoothing_radius)
    {
};

void ShapeHandler::obtain_affected_nodes(){
    const bool rigid = (this->mesh->proj_data->contact_data.contact_type == FiniteElement::ContactType::RIGID);
    const size_t node_num = this->mesh->elem_info->get_nodes_per_element();
    const size_t bnode_num = this->mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = this->mesh->elem_info->get_dof_per_node();
    const size_t s_size = this->mesh->elem_info->get_D_dimension();
    const size_t dim =
        (this->mesh->elem_info->get_problem_type() == utils::ProblemType::PROBLEM_TYPE_2D)
        ? 2 : 3;

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

    // Find out which nodes should be considered the same one for the purposes
    // of this function if interface is non-rigid
    this->merged_nodes.reserve(this->mesh->boundary_node_list.size());
    if(rigid){
        size_t bid = 0;
        for(auto& bn:this->mesh->boundary_node_list){
            SuperimposedNodes s{bid, {bn}};
            this->merged_nodes.push_back(std::move(s));
            this->merged_nodes_mapping[bn->id] = &this->merged_nodes.back();
            ++bid;
        }
    } else {
        logger::log_assert(this->mesh->to_rigid_map.size() > 0,
                logger::ERROR,
                "to_rigid_map is being cleared somewhere");
        size_t bid = 0;
        for(auto& bn:this->mesh->boundary_node_list){
            MeshNode* rn = this->mesh->to_rigid_map[bn->id];
            if(bn == rn){
                auto existing = this->merged_nodes_mapping.find(bn->id);
                if(existing == this->merged_nodes_mapping.end()){
                    SuperimposedNodes s{bid, {bn}};
                    this->merged_nodes.push_back(std::move(s));
                    this->merged_nodes_mapping[bn->id] = &this->merged_nodes.back();
                    ++bid;
                } else {
                    existing->second->nodes.push_back(bn);
                    this->merged_nodes_mapping[bn->id] = existing->second;
                }
            } else {
                auto existing = this->merged_nodes_mapping.find(rn->id);
                if(existing == this->merged_nodes_mapping.end()){
                    SuperimposedNodes s{bid, {bn, rn}};
                    this->merged_nodes.push_back(std::move(s));
                    this->merged_nodes_mapping[bn->id] = &this->merged_nodes.back();
                    this->merged_nodes_mapping[rn->id] = &this->merged_nodes.back();
                    ++bid;
                } else {
                    existing->second->nodes.push_back(bn);
                    this->merged_nodes_mapping[bn->id] = existing->second;
                }
            }
        }
    }

    auto short_list = this->apply_op(root_op.get());

    std::set<SuperimposedNodes*> final_list;

    // Find out which nodes will have their position optimized
    this->optimized_nodes.reserve(short_list.size());
    std::set<BoundaryElement*> elem_set;
    for(auto& ni:short_list){
        size_t af_id = 0;
        while(af_id < affected.size()){
            if(this->mesh->boundary_node_list[af_id]->id == ni->nodes[0]->id){
                break;
            }
            ++af_id;
        }
        if(affected[af_id]){
            final_list.insert(ni);
            std::vector<size_t> node_ids;
            node_ids.reserve(ni->nodes.size());
            std::vector<AffectedElement> elems;
            size_t e_size = 0;
            for(auto curr_node:ni->nodes){
                const size_t curr_id = curr_node->id;
                node_ids.push_back(curr_id);
                const auto& range = this->mesh->inverse_mesh.equal_range(curr_id);
                for(auto it = range.first; it != range.second; ++it){
                    ++e_size;
                }
            }
            elems.reserve(e_size);
            for(auto curr_node:ni->nodes){
                const size_t curr_id = curr_node->id;
                const auto& range = this->mesh->inverse_mesh.equal_range(curr_id);
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
    const auto shape_elem = this->mesh->elem_info->get_shape_element_info();
    std::set<ShapeMeshElement*> shape_mesh_tmp;

    {
        std::vector<MeshNode*> nodes(bnode_num);
        for(size_t i = 0; i < bnode_num; ++i){
            nodes[i] = this->mesh->node_list[this->boundary_elements[0]->nodes[i]->id].get();
        }
        this->bound_to_shape_mapping[0] = 0;
        ElementShape es{std::move(nodes), this->boundary_elements[0]->normal};
        shape_mesh_tmp.insert(shape_elem->make_element(std::move(es)));
    }

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
        shape_mesh_tmp.insert(shape_elem->make_element(std::move(es)));
        ++shape_id;
    }
    this->shape_mesh.insert(this->shape_mesh.begin(), shape_mesh_tmp.begin(), shape_mesh_tmp.end());
    shape_mesh_tmp.clear();


    // Save original coordinates (for visualization)
    for(size_t i = 0; i < this->mesh->node_list.size(); ++i){
        const auto& n = this->mesh->node_list[i];
        for(size_t j = 0; j < 3; ++j){
            this->original_points[3*i + j] = n->point.Coord(1+j);
        }
    }

    // Generate geometry clusters for node displacement spreading within each
    // geometry
    const size_t N = this->geometries.size();
    const size_t mat_size = (N-1)*N/2;
    const auto to_int_pos = [N](size_t i, size_t j) -> size_t{
        const size_t n = N - 1;
        const size_t ni = n - i;
        --j;
        return (n*(n+1)/2 - ni*(ni+1)/2) + j;
    };
    std::vector<std::set<SuperimposedNodes*>> node_list(N);
    std::vector<std::set<SuperimposedNodes*>> intersection_matrix(mat_size);
    for(size_t i = 0; i < this->geometries.size(); ++i){
        const auto g = this->geometries[i];
        auto& nodes = node_list[i];
        for(auto bn:g->boundary_node_list){
            nodes.insert(this->merged_nodes_mapping.at(bn->id));
        }
    }
    for(size_t i = 0; i < this->geometries.size(); ++i){
        const auto& s1 = node_list[i];
        for(size_t j = i + 1; j < this->geometries.size(); ++j){
            const auto& s2 = node_list[j];
            std::vector<SuperimposedNodes*> result_vec(std::min(s1.size(), s2.size()));
            auto& result = intersection_matrix[to_int_pos(i, j)];
            auto input_end = std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), result_vec.begin());
            result.insert(result_vec.begin(), input_end);
            result_vec.clear();
        }
    }
    std::map<size_t, GeometryCluster*> cluster_mapping;
    for(size_t i = 0; i < this->geometries.size(); ++i){
        auto gcpos = cluster_mapping.find(i);
        GeometryCluster* current_cluster = nullptr;
        if(gcpos == cluster_mapping.end()){
            this->clusters.emplace_back();
            current_cluster = &this->clusters.back();
            current_cluster->geometries.insert(this->geometries[i]);
            cluster_mapping[i] = current_cluster;
        } else {
            current_cluster = gcpos->second;
        }
        for(size_t j = i + 1; j < this->geometries.size(); ++j){
            auto find_geom = current_cluster->geometries.find(this->geometries[j]);
            if(find_geom == current_cluster->geometries.end()){
                auto& curr_set = intersection_matrix[to_int_pos(i, j)];
                if(!std::includes(final_list.begin(), final_list.end(), curr_set.begin(), curr_set.end())){
                    auto gcpos2 = cluster_mapping.find(j);
                    if(gcpos2 == cluster_mapping.end()){
                        current_cluster->geometries.insert(this->geometries[j]);
                        cluster_mapping[j] = current_cluster;
                    } else if(gcpos2->second != current_cluster) {
                        // Merge clusters if there is a "path" from one to
                        // another due to lack of intersections
                        GeometryCluster* correct_cluster = gcpos2->second;
                        correct_cluster->geometries.insert(current_cluster->geometries.begin(),
                                                           current_cluster->geometries.end());

                        for(size_t k = 0; k < this->clusters.size(); ++k){
                            if(&this->clusters[k] == current_cluster){
                                this->clusters.erase(this->clusters.begin() + k);
                                break;
                            }
                        }
                        current_cluster = correct_cluster;
                    }
                }
            } 
        }
    }
    // Check if there is no further need to merge clusters. Algorithm for
    // merging them is currently not implemented, so this just results in
    // an error.
    for(size_t i = 0; i < this->clusters.size(); ++i){
        auto& s1 = this->clusters[i].geometries;
        for(size_t j = i + 1; j < this->clusters.size(); ++j){
            auto& s2 = this->clusters[j].geometries;
            std::vector<Geometry*> intersec(std::min(s1.size(), s2.size()));
            auto pos = std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                                    intersec.begin());
            std::set<Geometry*> result(intersec.begin(), pos);
            logger::log_assert(result.size() == 0, logger::ERROR,
                    "there is an intersection between clusters");
        }
    }

    // Find optimized nodes in boundary
    std::set<SuperimposedNodes*> inter_geom;
    for(auto& e:mesh->inter_geometry_boundary){
        for(size_t i = 0; i < bnode_num; ++i){
            inter_geom.insert(this->merged_nodes_mapping[e->nodes[i]->id]);
        }
    }

    const auto get_rel_node_id = 
        [](const Node *const *n, const size_t num, const gp_Pnt& p)->size_t{
            for(size_t i = 0; i < num; ++i){
                if(n[i]->point.IsEqual(p, Precision::Confusion())){
                    return i;
                }
            }
            return 0;
        };
    const auto get_rel_bnode_id = 
        [](const MeshNode *const *n, const size_t num, const gp_Pnt& p)->size_t{
            for(size_t i = 0; i < num; ++i){
                if(n[i]->point.IsEqual(p, Precision::Confusion())){
                    return i;
                }
            }
            return 0;
        };

    // If contact is not rigid, handle node duplication
    std::list<AffectedContactNode> contact_tmp;
    for(size_t i = 0; i < this->optimized_nodes.size(); ++i){
        const auto& a = this->optimized_nodes[i];
        const auto ni = a.node_ids[0];
        if(inter_geom.contains(this->merged_nodes_mapping[ni])){
            const gp_Pnt& p = mesh->node_list[ni]->point;
            std::vector<AffectedPairedBoundary> elems_tmp;
            elems_tmp.reserve(a.elements.size());
            for(auto& e:a.elements){
                for(auto& be:this->mesh->paired_boundary){
                    if(e.e == be.b1->parent || e.e == be.b2->parent){
                        AffectedPairedBoundary ab;
                        ab.e = &be;
                        ab.en1 = get_rel_node_id(be.b1->parent->nodes, node_num, p);
                        ab.en2 = get_rel_node_id(be.b2->parent->nodes, node_num, p);
                        ab.bn1 = get_rel_bnode_id(be.b1->nodes, bnode_num, p);
                        ab.bn2 = get_rel_bnode_id(be.b2->nodes, bnode_num, p);
                        elems_tmp.push_back(std::move(ab));
                        break;
                    }
                }
            }
            elems_tmp.shrink_to_fit();
            contact_tmp.emplace_back(i, std::move(elems_tmp));
        }
    }
    this->optimized_contact_nodes.reserve(contact_tmp.size());
    for(auto& c:contact_tmp){
        this->optimized_contact_nodes.emplace_back(c.id, c.elements);
    }
    contact_tmp.clear();

    // Generate linear problems
    const math::Matrix D(math::Matrix::identity(s_size));
    const math::Matrix A(math::Matrix::identity(dim));
    for(auto& cluster:this->clusters){
        cluster.solver = std::make_unique<general_solver::PETScGeneralPCG>();

        // Go through boundary (in parallel)
        //      Set as id -1 the ones on boundary conditions (but not on inter_geom)
        // Go through domain
        //      If id != -1:
        //          Set id

        size_t num_bound_nodes = 0;
        for(auto& g:cluster.geometries){ 
            num_bound_nodes += g->boundary_node_list.size();
        }

        std::vector<bool> fixed_node(num_bound_nodes, false);
        size_t fixed_offset = 0;
        for(auto& g:cluster.geometries){ 
            #pragma omp parallel for
            for(size_t i = 0; i < g->boundary_node_list.size(); ++i){
                const auto& b = g->boundary_node_list[i];
                const size_t bid = i + fixed_offset;
                for(const auto& f:this->mesh->proj_data->forces){
                    if(f.S.is_inside(b->point)){
                        fixed_node[bid] = true;
                        break;
                    }
                }
                if(!fixed_node[bid]){
                    for(const auto& f:this->mesh->proj_data->supports){
                        if(f.S.is_inside(b->point)){
                            fixed_node[bid] = true;
                            break;
                        }
                    }
                }
                if(!fixed_node[bid]){
                    for(const auto& f:this->mesh->proj_data->springs){
                        if(f.S.is_inside(b->point)){
                            fixed_node[bid] = true;
                            break;
                        }
                    }
                }
                if(!fixed_node[bid]){
                    for(const auto& f:this->mesh->proj_data->internal_loads){
                        if(f.S.is_inside(b->point)){
                            fixed_node[bid] = true;
                            break;
                        }
                    }
                }
                if(!fixed_node[bid]){
                    bool found = false;
                    for(const auto& f:this->optimized_nodes){
                        for(auto& n:f.node_ids){
                            if(b->id == n){
                                found = true;
                                fixed_node[bid] = true;
                                break;
                            }
                        }
                        if(found){
                            break;
                        }
                    }
                }
            }
            fixed_offset += g->boundary_node_list.size();
        }
        fixed_offset = 0;
        for(auto& g:cluster.geometries){ 
            for(size_t i = 0; i < g->boundary_node_list.size(); ++i){
                const auto& b = g->boundary_node_list[i];
                const size_t bid = i + fixed_offset;
                if(fixed_node[bid]){
                    cluster.id_mapping[b->id] = -1;
                }
            }
            fixed_offset += g->boundary_node_list.size();
        }
        fixed_node.clear();
        std::set<SuperimposedNodes*> geom_union;
        for(auto& g:cluster.geometries){
            geom_union.insert(node_list[g->id].begin(), node_list[g->id].end());
        }
        size_t idm = 0;
        for(const auto& n:geom_union){
            bool fixed = false;
            for(auto& mn:n->nodes){
                if(cluster.id_mapping[mn->id] < 0){
                    fixed = true;
                    break;
                }
            }
            if(!fixed){
                for(auto& mn:n->nodes){
                    cluster.id_mapping[mn->id] = idm;
                }
                ++idm;
            }
        }
        for(auto& g:cluster.geometries){
            for(auto& n:g->node_list){
                if(!cluster.id_mapping.contains(n->id)){
                    cluster.id_mapping[n->id] = idm;
                    ++idm;
                }
            }
        }

        cluster.matrix_width = dof*idm;
        cluster.solver->initialize_matrix(true, cluster.matrix_width);
        cluster.b.resize(cluster.matrix_width);
    }

    this->merged_nodes_mapping.clear();
    this->merged_nodes.clear();

    this->helmholtz_solver = std::make_unique<general_solver::MUMPSGeneral>();
    this->helmholtz_solver->initialize_matrix(true, this->optimized_nodes.size()*dof);
    //this->helmholtz_solver->initialize_matrix(true, this->optimized_nodes.size());

    this->generate_helmholtz();
}

void ShapeHandler::generate_helmholtz(){
    const size_t dof = this->mesh->elem_info->get_dof_per_node();
    const size_t bnum = this->mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dim =
        (this->mesh->elem_info->get_problem_type() == utils::ProblemType::PROBLEM_TYPE_2D)
        ? 2 : 3;

    const math::Matrix D(math::Matrix::identity(dim*dim));

    std::vector<long> bpos(bnum*dof);
    this->helmholtz_solver->make_zero();
    for(const auto& e:this->shape_mesh){
        for(size_t i = 0; i < bnum; ++i){
            const auto id = this->optimized_nodes_mapping[e->nodes[i]->id];
            for(size_t dim_i = 0; dim_i < dof; ++dim_i){
                bpos[i*dof + dim_i] = dof*id + dim_i;
            }
        }
        const auto a = e->absorption_Ndof();
        const auto k = (this->r*this->r)*e->diffusion_Ndof(D) + a;
        this->helmholtz_solver->add_element(k, bpos);
    }
    logger::quick_log("Smoothing...");
    this->helmholtz_solver->compute();
}

void ShapeHandler::filter_displacement(std::vector<double>& dx){
    logger::quick_log("Applying shape changes...");
    const size_t dof = this->mesh->elem_info->get_dof_per_node();
    const size_t bnum = this->mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dim =
        (this->mesh->elem_info->get_problem_type() == utils::ProblemType::PROBLEM_TYPE_2D)
        ? 2 : 3;

    // Update other nodes
    const math::Matrix D(math::Matrix::identity(dim*dim));
    //math::Matrix D(dim*dim, dim*dim, 0);
    //D(0,0) = 1;
    //D(4,4) = 1;
    //D(8,8) = 1;
    auto A = math::Matrix::identity(dim);
    /**
     * Smoothen displacements
     */
    math::Vector dx_be(bnum*dof);
    std::vector<double> dx_f(dx.size(), 0);
    for(const auto& e:this->shape_mesh){
        for(size_t i = 0; i < bnum; ++i){
            const auto id = this->optimized_nodes_mapping[e->nodes[i]->id];
            for(size_t dim_i = 0; dim_i < dof; ++dim_i){
                dx_be[i*dof + dim_i] = dx[dof*id + dim_i];
            }
        }
        const auto a = e->absorption_Ndof();
        const auto f = a*dx_be;
        for(size_t i = 0; i < bnum; ++i){
            const auto id = this->optimized_nodes_mapping[e->nodes[i]->id];
            for(size_t dim_i = 0; dim_i < dof; ++dim_i){
                dx_f[dof*id + dim_i] += f[i*dof + dim_i];
            }
        }
    }
    logger::quick_log("Smoothing...");
    this->helmholtz_solver->solve(dx_f);
    std::copy(dx_f.begin(), dx_f.end(), dx.begin());
}
    
void ShapeHandler::update_nodes(const std::vector<double>& dx){
    logger::quick_log("Applying shape changes...");
    const size_t dof = this->mesh->elem_info->get_dof_per_node();
    const size_t bnum = this->mesh->elem_info->get_boundary_nodes_per_element();
    const size_t num = this->mesh->elem_info->get_nodes_per_element();
    const size_t dim =
        (this->mesh->elem_info->get_problem_type() == utils::ProblemType::PROBLEM_TYPE_2D)
        ? 2 : 3;

    // Update other nodes
    const math::Matrix D(math::Matrix::identity(dim*dim));
    //math::Matrix D(dim*dim, dim*dim, 0);
    //D(0,0) = 1;
    //D(4,4) = 1;
    //D(8,8) = 1;
    auto A = math::Matrix::identity(dim);

    /**
     * Diffuse displacements to domain nodes
     */
    math::Vector bn_vals(dof*num, 0);
    math::Vector Fe(dof*num, 0);
    std::vector<gp_Pnt> pnts(bnum);
    std::vector<long> pos(num*dof);
    for(auto& cluster:this->clusters){
        std::fill(cluster.b.begin(), cluster.b.end(), 0);

        for(const auto g:cluster.geometries){
            // Generate global matrix
            for(const auto& e:g->mesh){
                for(size_t i = 0; i < num; ++i){
                    const auto id = cluster.id_mapping[e->nodes[i]->id];
                    for(size_t dim_i = 0; dim_i < dof; ++dim_i){
                        if(id < 0){
                            pos[i*dof + dim_i] = -1;
                        } else {
                            pos[i*dof + dim_i] = dof*id + dim_i;
                        }
                    }
                }
                const auto k = e->diffusion_dim_dof(this->mesh->thickness, D);
                cluster.solver->add_element(k, pos);

                if(this->elem_to_affected_node_mapping.contains(e.get())){
                    for(size_t i = 0; i < num; ++i){
                        if(this->optimized_nodes_mapping.contains(e->nodes[i]->id)){
                            const auto opt_id = this->optimized_nodes_mapping[e->nodes[i]->id];
                            for(size_t dim_i = 0; dim_i < dof; ++dim_i){
                                bn_vals[dof*i + dim_i] = dx[dof*opt_id + dim_i];
                            }
                        }
                    }
                    const math::Vector Fe = k*bn_vals;
                    for(size_t j = 0; j < num; ++j){
                        const long global_id = cluster.id_mapping[e->nodes[j]->id];
                        if(global_id > -1){
                            for(size_t dim_i = 0; dim_i < dof; ++dim_i){
                                cluster.b[dof*global_id + dim_i] -= Fe[dof*j + dim_i];
                            }
                        }
                    }
                    bn_vals.fill(0);
                }
            }
            // Apply outer boundary constraint
            for(const auto& b:g->boundary_mesh){
                // Check if:
                // - Not in optimization region
                // - Not in intergeometry boundary (merged boundaries
                //   within cluster are not subject to it, even if not
                //   being optimized)
                const auto found_pos = std::find(this->mesh->inter_geometry_boundary.begin(), this->mesh->inter_geometry_boundary.end(), b);
                const bool is_in_intergeom = (found_pos != this->mesh->inter_geometry_boundary.end());
                bool is_optimized = false;
                if(this->elem_to_affected_node_mapping.contains(b->parent)){
                    const auto& nodes = this->elem_to_affected_node_mapping[b->parent];
                    if(nodes.size() == bnum){
                        is_optimized = true;
                    }
                }
                if(!is_in_intergeom && !is_optimized){
                    for(size_t i = 0; i < bnum; ++i){
                        pnts[i] = b->nodes[i]->point;
                    }
                    for(size_t dim_i = 0; dim_i < dim; ++dim_i){
                        for(size_t dim_j = 0; dim_j < dim; ++dim_j){
                            A(dim_i,dim_j) = b->normal.Coord(1+dim_i)*b->normal.Coord(1+dim_j);
                        }
                    }
                    for(size_t i = 0; i < num; ++i){
                        const auto id = cluster.id_mapping[b->parent->nodes[i]->id];
                        for(size_t dim_i = 0; dim_i < dof; ++dim_i){
                            if(id < 0){
                                pos[i*dof + dim_i] = -1;
                            } else {
                                pos[i*dof + dim_i] = dof*id + dim_i;
                            }
                        }
                    }
                    const auto k = this->MU*b->parent->get_R(A, this->mesh->thickness, pnts);
                    cluster.solver->add_element(k, pos);
                }
            }
        }
        
        cluster.solver->compute();
        cluster.solver->solve(cluster.b);
        cluster.solver->make_zero();

        // Apply nodal displacements
        for(const auto g:cluster.geometries){
            for(const auto& ni:g->node_list){
                const long global_id = cluster.id_mapping[ni->id];
                // GT9 breaks this, but I'll fix it later if needed
                if(global_id > -1){
                    for(size_t dim_i = 0; dim_i < dim; ++dim_i){
                        ni->point.SetCoord(1+dim_i, ni->point.Coord(1+dim_i) + cluster.b[dim*global_id + dim_i]);
                    }
                }
            }
        }
    }

    // Update boundary coordinates
    for(size_t i = 0; i < this->optimized_nodes.size(); ++i){
        auto& nids = this->optimized_nodes[i].node_ids;
        for(auto& nid:nids){
            auto& n = this->mesh->node_list[nid];
            for(size_t d = 0; d < dim; ++d){
                n->point.SetCoord(1+d, n->point.Coord(1+d) + dx[i*dim + d]);
            }
        }
    }
    for(auto& b:this->mesh->boundary_elements){
        b.update_normal(bnum, this->mesh->proj_data->type);
    }

    // Update view
    for(size_t i = 0; i < this->mesh->node_list.size(); ++i){
        const auto& n = this->mesh->node_list[i];
        for(size_t j = 0; j < dim; ++j){
            this->shape_displacement[dim*i + j] = n->point.Coord(1+j) - this->original_points[dim*i + j];
        }
    }

    // Update all coefficients
    for(auto& g:this->mesh->geometries){
        #pragma omp parallel for
        for(auto& e:g->mesh){
            e->calculate_coefficients();
        }
    }

    this->generate_helmholtz();

    logger::quick_log("Done.");
}

void ShapeHandler::filter_gradient(std::vector<double>& df){
    const size_t dof = this->mesh->elem_info->get_dof_per_node();
    const size_t bnum = this->mesh->elem_info->get_boundary_nodes_per_element();

    math::Vector dx_be(bnum*dof);
    std::vector<double> dx_f(df.size(), 0);
    this->helmholtz_solver->solve(df);
    for(const auto& e:this->shape_mesh){
        for(size_t i = 0; i < bnum; ++i){
            const auto id = this->optimized_nodes_mapping[e->nodes[i]->id];
            for(size_t dim_i = 0; dim_i < dof; ++dim_i){
                dx_be[i*dof + dim_i] = df[dof*id + dim_i];
            }
        }
        const auto a = e->absorption_Ndof();
        const auto f = a*dx_be;
        for(size_t i = 0; i < bnum; ++i){
            const auto id = this->optimized_nodes_mapping[e->nodes[i]->id];
            for(size_t dim_i = 0; dim_i < dof; ++dim_i){
                dx_f[dof*id + dim_i] += f[i*dof + dim_i];
            }
        }
    }
    std::copy(dx_f.begin(), dx_f.end(), df.begin());
}

std::set<ShapeHandler::SuperimposedNodes*> ShapeHandler::apply_op(shape_op::ShapeOp* op) const{
    switch(op->get_type()){
        case shape_op::Code::GEOMETRY:{
            logger::log_assert(op->get_id() < this->geometries.size(), logger::ERROR,
                    "unknown geometry id: {}", op->get_id());
            const auto g = this->geometries[op->get_id()];
            std::set<SuperimposedNodes*> nodes;//(g->boundary_node_list.begin(), g->boundary_node_list.end());
            for(auto bn:g->boundary_node_list){
                nodes.insert(this->merged_nodes_mapping.at(bn->id));
            }

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
            std::vector<SuperimposedNodes*> result_vec;
            std::set<SuperimposedNodes*> result;
            std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(result_vec));
            s1.clear();
            s2.clear();
            result.insert(result_vec.begin(), result_vec.end());
            result_vec.clear();
            return result;
        }
        case shape_op::Code::DIFFERENCE:{
            auto s1 = this->apply_op(op->first());
            auto s2 = this->apply_op(op->second());
            std::vector<SuperimposedNodes*> result_vec;
            std::set<SuperimposedNodes*> result;
            std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(result_vec));
            s1.clear();
            s2.clear();
            result.insert(result_vec.begin(), result_vec.end());
            result_vec.clear();
            return result;
        }
        case shape_op::Code::SHELL:{
            logger::log_assert(false, logger::ERROR,
                    "SHELL shape operation is currently not implemented");
            return std::set<SuperimposedNodes*>();
        }
    }

    return std::set<SuperimposedNodes*>();
}
