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

#include "meshing.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include "project_data.hpp"
#include <BRepBuilderAPI_Transform.hxx>
#include <algorithm>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRep_Builder.hxx>
#include <BOPAlgo_Splitter.hxx>
#include <BOPAlgo_Builder.hxx>
#include <limits>
#include <memory>
#include <set>
#include <queue>
#include <BRepBuilderAPI_Copy.hxx>
#include <vector>


void Meshing::generate_elements(const TopoDS_Shape& shape,
                                const std::vector<size_t>& geom_elem_mapping, 
                                const std::vector<size_t>& elem_node_tags, 
                                const std::vector<size_t>& bound_elem_node_tags,
                                std::unordered_map<size_t, MeshNode*>& id_map,
                                std::unordered_map<size_t, size_t>& duplicate_map,
                                const bool deduplicate,
                                const bool boundary_condition_inside){

    this->orig_shape = shape;
    logger::quick_log("Generating elements and preparing for finite element analysis...");

    std::unordered_map<size_t, MeshNode*> original_map = id_map;
    const bool rigid = (this->proj_data->contact_data.contact_type == FiniteElement::ContactType::RIGID);

    const bool delete_duplicated = deduplicate && rigid;
    if(this->geometries.size() > 1){
        this->deduplicate(id_map, duplicate_map, delete_duplicated);
    }

    const size_t nodes_per_elem = this->elem_info->get_nodes_per_element();
    const size_t bound_nodes_per_elem = this->elem_info->get_boundary_nodes_per_element();

    {
        std::vector<ElementShape> list;
        if(rigid){
            // Generate element shapes from trimmed node list
            list = this->generate_element_shapes(elem_node_tags, nodes_per_elem, id_map);
        } else {
            // Generate element shapes from original node list
            // Currently no per-geometry boundary setting, applies to whole list
            list = this->generate_element_shapes(elem_node_tags, nodes_per_elem, original_map);
        }
        // Re-id nodes and remove unused nodes from map
        this->optimize(list, &id_map);
        if(!rigid){
            for(const auto& n:id_map){
                const size_t oid = n.first;
                MeshNode* dedup_n = n.second;
                MeshNode* prev_n = original_map[oid];
                const size_t curr_id = prev_n->id;
                this->to_rigid_map[curr_id] = dedup_n;
            }
        }
        auto elements = this->create_element_list(list, this->elem_info);
        list.clear();
        populate_inverse_mesh(elements);
        this->distribute_elements(geom_elem_mapping, elements);
        elements.clear();
        // Generate element ids
        size_t e_id = 0;
        for(auto& g:this->geometries){
            for(auto& e:g->mesh){
                e->id = e_id;
                ++e_id;
            }
        }
    }

    {
        std::vector<ElementShape> bound_list;
        if(rigid){
            // Generate element shapes from trimmed node list
            bound_list = this->generate_element_shapes(bound_elem_node_tags, bound_nodes_per_elem, id_map, true);
        } else {
            // Generate element shapes from original node list
            // Currently no per-geometry boundary setting, applies to whole list
            bound_list = this->generate_element_shapes(bound_elem_node_tags, bound_nodes_per_elem, original_map, true);
        }
        this->populate_boundary_elements(bound_list, boundary_condition_inside);
        this->distribute_boundary_elements();
    }
    this->distribute_node_pointers();
    logger::quick_log("Done.");
}

void Meshing::apply_boundary_conditions(const std::vector<Force>& forces, 
                                        const std::vector<Support>& supports,
                                        std::vector<Spring>& springs,
                                        std::vector<InternalLoads>& internal_loads,
                                        std::vector<SubProblem>& sub_problems){

    logger::quick_log("Applying boundary conditions...");
    const bool rigid = (this->proj_data->contact_data.contact_type != FiniteElement::ContactType::RIGID);
    if (!rigid){
        logger::log_assert(this->to_rigid_map.size() > 0,
                logger::ERROR,
                "ensure elements are generated before applying boundary conditions.");
    }

    const size_t dof = this->elem_info->get_dof_per_node();
    const size_t vec_size = this->node_list.size()*dof;
    this->max_dofs = vec_size;

    this->sub_problems = &sub_problems;

    this->node_positions.clear();
    this->node_positions.resize(sub_problems.size());
    this->dofs_per_subproblem.clear();
    this->dofs_per_subproblem.resize(sub_problems.size());
    this->load_vector.clear();
    this->load_vector.resize(sub_problems.size());
    for(auto& v:this->node_positions){
        v.resize(vec_size, 0);
    }
    this->global_load_vector.resize(vec_size);
    // TODO: revert back to old way of removing nodes?
    // Factorization is faster and uses less memory, but it would require
    // passing `node_positions` and handle used instances of node->u_pos

    // Number positioning
    #pragma omp parallel for
    for(size_t i = 0; i < this->node_list.size(); ++i){
        auto& n = this->node_list[i];
        for(size_t j = 0; j < dof; ++j){
            n->u_pos[j] = n->id*dof + j;
        }
    }

    if(supports.size() > 0){
        logger::quick_log("supports");
        this->apply_supports();
    }

    for(size_t it = 0; it < this->node_positions.size(); ++it){
        auto& v = this->node_positions[it];
        long current = 0;
        //if(rigid){
        for(auto& i:v){
            if(i >= 0){
                i = current;
                ++current;
            }
        }
        /*
        } else {
            std::vector<bool> pos_set(this->node_list.size(), false);
            for(size_t i = 0; i < this->node_list.size(); ++i){
                const auto& n = this->node_list[i];
                const size_t n_id = n->id;
                const size_t on_id = this->to_rigid_map[n->id]->id;
                if(pos_set[i]){
                    continue;
                }
                if(n_id == on_id){
                    pos_set[i] = true;
                    for(size_t d = 0; d < dof; ++d){
                        if(v[n_id*dof + d] >= 0){
                            v[n_id*dof + d] = current;
                            ++current;
                        }
                    }
                } else if(on_id < n_id){
                    pos_set[i] = true;
                    for(size_t d = 0; d < dof; ++d){
                        v[n_id*dof + d] = v[on_id*dof + d];
                    }
                } else {
                    pos_set[i] = true;
                    pos_set[on_id] = true;
                    for(size_t d = 0; d < dof; ++d){
                        if(v[n_id*dof + d] >= 0){
                            v[n_id*dof + d] = current;
                            v[on_id*dof + d] = current;
                            ++current;
                        }
                    }
                }
            }
        }*/
        this->dofs_per_subproblem[it] = current;
        this->load_vector[it].resize(current, 0);
    }

    logger::quick_log("loads");
    this->generate_load_vector(this->orig_shape);

    if(springs.size() > 0){
        logger::quick_log("springs");
        this->apply_springs(springs);
    }
    this->springs = &springs;

    if(internal_loads.size() > 0){
        logger::quick_log("internal_loads");
        this->apply_internal_loads(internal_loads);
    }
    this->internal_loads = &internal_loads;

    for(size_t i = 0; i < this->node_positions.size(); ++i){
        const auto& n = this->node_positions[i];
        const auto& l = this->load_vector[i];
        for(size_t j = 0; j < n.size(); ++j){
            if(n[j] >= 0){
                this->global_load_vector[j] += l[n[j]];
            }
        }
    }

    logger::quick_log("Done.");
}

std::vector<std::unique_ptr<MeshElement>> Meshing::create_element_list(
                               const std::vector<ElementShape>& base_mesh, 
                               const MeshElementFactory * const elem_info) const{
    std::vector<std::unique_ptr<MeshElement>> element_list(base_mesh.size());

    #pragma omp parallel for
    for(size_t i = 0; i < base_mesh.size(); ++i){
        element_list[i].reset(elem_info->make_element(base_mesh[i]));
    }

    return element_list;
}

void Meshing::populate_inverse_mesh(const std::vector<std::unique_ptr<MeshElement>>& element_list){
    const size_t N = this->elem_info->get_nodes_per_element();
    this->inverse_mesh.clear();
    for(const auto& e : element_list){
        for(size_t i = 0; i < N; ++i){
            this->inverse_mesh.emplace(e->nodes[i]->id, e.get());
        }
    }
}

void Meshing::populate_boundary_elements(const std::vector<ElementShape>& boundary_base_mesh,
                                         const bool boundary_condition_inside){
    (void) boundary_condition_inside;
    const size_t N = this->elem_info->get_boundary_nodes_per_element();
    this->boundary_elements.clear();
    this->boundary_elements.reserve(boundary_base_mesh.size());
    std::vector<size_t> inter_geom;
    // If contact type is rigid, overlapping nodes have been merged, so finding
    // the intergeometry boundary is a matter of looking for intersections
    // using the inverse mesh. Such boundary involves two boundary elements
    // sharing the same boundary nodes, but having opposite normals.
    for(size_t i = 0; i < boundary_base_mesh.size(); ++i){
        const auto& b = boundary_base_mesh[i];
        std::set<MeshElement*> common_nodes;
        const auto eq_range = this->inverse_mesh.equal_range(b.nodes[0]->id);
        for(auto k = eq_range.first; k != eq_range.second; ++k){
            common_nodes.insert(k->second);
        }
        for(size_t j = 1; j < N; ++j){
            std::set<MeshElement*> tmp_common_nodes;
            std::set<MeshElement*> tmp_comp;
            const auto eq_range = this->inverse_mesh.equal_range(b.nodes[j]->id);
            for(auto k = eq_range.first; k != eq_range.second; ++k){
                tmp_comp.insert(k->second);
            }
            std::set_intersection(common_nodes.begin(), common_nodes.end(),
                                  tmp_comp.begin(), tmp_comp.end(),
                                  std::inserter(tmp_common_nodes, tmp_common_nodes.begin()));
            common_nodes = std::move(tmp_common_nodes);
        }
        if(common_nodes.size() > 0){
            //if(boundary_condition_inside){
            //    this->boundary_elements.emplace_back(b.nodes, *common_nodes.begin(), b.normal);
            //} else 
            if(common_nodes.size() >= 1){
                double x = 0, y = 0, z = 0;
                for(auto& n:b.nodes){
                    x += n->point.X();
                    y += n->point.Y();
                    z += n->point.Z();
                }
                x /= b.nodes.size();
                y /= b.nodes.size();
                z /= b.nodes.size();
                gp_Pnt bc(x, y, z);
                gp_Pnt ec((*common_nodes.begin())->get_centroid());
                gp_Dir d(gp_Vec(bc, ec));
                gp_Dir n(b.normal);
                const double mult = (d.Dot(n) < 0) ? 1 : -1;
                if(common_nodes.size() == 1){
                    auto it = common_nodes.begin();
                    size_t geom_id = 0;
                    for(const auto& g:this->geometries){
                        bool found = false;
                        for(const auto& e:g->mesh){
                            if(e.get() == *it){
                                found = true;
                                break;
                            }
                        }
                        if(found){
                            break;
                        }
                        ++geom_id;
                    }
                    this->boundary_elements.emplace_back(b.nodes, *it, mult*n, geom_id);
                } else if(common_nodes.size() == 2){
                    auto it = common_nodes.begin();
                    bool found = false;
                    for(const auto& be:this->boundary_elements){
                        if(be.parent == *it && std::abs(be.normal.Dot(n)) >= 1.0 - Precision::Confusion()){
                            found = true;
                            break;
                        }
                    }
                    if(!found){
                        size_t geom_id = 0;
                        for(const auto& g:this->geometries){
                            bool found = false;
                            for(const auto& e:g->mesh){
                                if(e.get() == *it){
                                    found = true;
                                    break;
                                }
                            }
                            if(found){
                                break;
                            }
                            ++geom_id;
                        }
                        this->boundary_elements.emplace_back(b.nodes, *it, mult*n, geom_id);
                        inter_geom.push_back(this->boundary_elements.size()-1);
                        ++it;
                        geom_id = 0;
                        for(const auto& g:this->geometries){
                            bool found = false;
                            for(const auto& e:g->mesh){
                                if(e.get() == *it){
                                    found = true;
                                    break;
                                }
                            }
                            if(found){
                                break;
                            }
                            ++geom_id;
                        }
                        this->boundary_elements.emplace_back(b.nodes, *it, -mult*n, geom_id);
                        inter_geom.push_back(this->boundary_elements.size()-1);
                    }
                } else {
                    logger::log_assert(false, logger::ERROR, "More than two elements are connected to a single boundary (this shouldn't happen).");
                }
            }
        }
    }
    // If contact type is not rigid, nodes have not been merged, so one needs
    // to use to_rigid_map to find overlapping nodes and generate intergeometry
    // boundary metadata
    if(this->proj_data->contact_data.contact_type != FiniteElement::ContactType::RIGID){
        // Generate inverse mesh for boundary
        for(auto& b:this->boundary_elements){
            for(size_t i = 0; i < N; ++i){
                this->boundary_inverse_mesh.emplace(b.nodes[i]->id, &b);
            }
        }
        std::vector<size_t> bnodes(N);
        for(auto& b:this->boundary_elements){
            bool found_equivalent = true;
            for(size_t i = 0; i < N; ++i){
                const auto rn = this->to_rigid_map[b.nodes[i]->id];
                // All nodes must have equivalents
                // If at the edge of intergeometry boundary, not all nodes
                // will have equivalents, so check that all fit into the
                // intergeometry region
                if(rn->id != b.nodes[i]->id){
                    bnodes[i] = rn->id;
                } else {
                    found_equivalent = false;
                    break;
                }
            }
            if(found_equivalent){
                std::set<BoundaryElement*> common_nodes;
                const auto eq_range = this->boundary_inverse_mesh.equal_range(bnodes[0]);
                for(auto k = eq_range.first; k != eq_range.second; ++k){
                    common_nodes.insert(k->second);
                }
                for(size_t j = 1; j < N; ++j){
                    std::set<BoundaryElement*> tmp_common_nodes;
                    std::set<BoundaryElement*> tmp_comp;
                    const auto eq_range = this->boundary_inverse_mesh.equal_range(bnodes[j]);
                    for(auto k = eq_range.first; k != eq_range.second; ++k){
                        tmp_comp.insert(k->second);
                    }
                    std::set_intersection(common_nodes.begin(), common_nodes.end(),
                                          tmp_comp.begin(), tmp_comp.end(),
                                          std::inserter(tmp_common_nodes, tmp_common_nodes.begin()));
                    common_nodes = std::move(tmp_common_nodes);
                }
                if(common_nodes.size() == 1){
                    auto e1 = &b;
                    auto e2 = *common_nodes.begin();
                    this->inter_geometry_boundary.push_back(e1);
                    this->inter_geometry_boundary.push_back(e2);
                    std::vector<MeshNode*> nodes(N);
                    for(size_t i = 0; i < N; ++i){
                        nodes[i] = this->node_list[e1->nodes[i]->id].get();
                        logger::log_assert(nodes[i]->id == e1->nodes[i]->id, logger::ERROR, "assumption is wrong");
                        this->lag_node_map[nodes[i]->id] = 0;
                    }
                    ElementShape es{std::move(nodes), e1->normal};
                    std::unique_ptr<ContactMeshElement> bme(nullptr);
                    if(this->proj_data->contact_data.contact_type == FiniteElement::FRICTIONLESS_DISPL){
                        bme.reset(this->elem_info->get_contact_element_info()->make_element(es, e1->parent, e2->parent));
                    }
                    this->paired_boundary.emplace_back(e1, e2, std::move(bme));

                } else {
                    logger::quick_log("????????", common_nodes.size());
                }
            }
        }
        // std::vector<size_t> bnodes(N);
        // // Not optimal, but easier to implement.
        // // This is not working correctly and I don't know why.
        // for(size_t i = 0; i < this->boundary_elements.size(); ++i){
        //     const BoundaryElement* b = &this->boundary_elements[i];
        //     const gp_Pnt c(b->get_centroid(N));
        //     //double dist = 1e100;
        //     //gp_Pnt c_test;
        //     for(size_t j = i+1; j < this->boundary_elements.size(); ++j){
        //         const BoundaryElement* b2 = &this->boundary_elements[j];
        //         const gp_Pnt c2(b2->get_centroid(N));
        //         //if(c.Distance(c2) < dist){
        //         //    c_test = c2;
        //         //    dist = c.Distance(c2);
        //         //}
        //         if(c.IsEqual(c2, Precision::Confusion())){
        //             inter_geom.push_back(i);
        //             inter_geom.push_back(j);
        //             break;
        //         }
        //     }
        //     //if(std::abs(c.X()) < 11 && std::abs(c.Z()) < 11 && std::abs(c.Y()) < 19){
        //     //    logger::quick_log(c.X(), c.Y(), c.Z(), dist);
        //     //}
        // }
    }
    logger::quick_log("bound_elems", this->boundary_elements.size());
    if(inter_geom.size() > 0){
        this->inter_geometry_boundary.resize(inter_geom.size());
        for(size_t i = 0; i < inter_geom.size(); ++i){
            this->inter_geometry_boundary[i] = &this->boundary_elements[inter_geom[i]];
        }
    } else {
        this->inter_geometry_boundary.shrink_to_fit();
        this->paired_boundary.shrink_to_fit();
    }
    logger::quick_log("inter_geom", this->inter_geometry_boundary.size());
    if(this->paired_boundary.size() > 0){
        long lag_id = 0;
        for(auto& m:this->lag_node_map){
            m.second = lag_id;
            ++lag_id;
        }
    }
}

void Meshing::apply_supports(){
    const size_t dof = this->elem_info->get_dof_per_node();

    for(size_t spid = 0; spid < this->sub_problems->size(); ++spid){
        const auto& p = this->sub_problems->at(spid);
        auto& node_pos = this->node_positions[spid];
        for(auto& s:p.supports){
            if(s->S.get_area() > 0){
                #pragma omp parallel for
                for(size_t i = 0; i < this->boundary_node_list.size(); ++i){
                    const auto& n = this->boundary_node_list[i];

                    if(s->S.is_inside(n->point)){
                        std::vector<bool> sup_pos = this->get_support_dof(*s, this->elem_info);
                        for(size_t j = 0; j < dof; ++j){
                            if(sup_pos[j]){
                                const size_t p = n->id*dof + j;
                                node_pos[p] = -1;
                            }
                        }
                    }
                }
            } else {
                const gp_Pnt c = s->S.get_centroid();
                const MeshNode* curr = nullptr;
                double dist = 1e100;
                for(size_t i = 0; i < this->node_list.size(); ++i){
                    const auto& n = this->node_list[i];

                    double d = c.Distance(n->point);
                    if(d < dist){
                        dist = d;
                        curr = n.get();
                    }
                }
                logger::quick_log(curr->point.X(), curr->point.Y(), curr->point.Z());
                std::vector<bool> sup_pos = this->get_support_dof(*s, this->elem_info);
                for(size_t j = 0; j < dof; ++j){
                    if(sup_pos[j]){
                        const size_t p = curr->id*dof + j;
                        node_pos[p] = -1;
                    }
                }
            }
        }
    }
}

void Meshing::apply_springs(std::vector<Spring>& springs){
    for(auto& s:springs){
        s.generate_mesh(this->boundary_elements);
    }
}

void Meshing::apply_internal_loads(std::vector<InternalLoads>& loads){
    const auto problem_type = this->elem_info->get_problem_type();
    for(auto& s:loads){
        s.calculate_curvature(this->boundary_elements);
    }
    for(size_t spid = 0; spid < this->sub_problems->size(); ++spid){
        const auto& p = this->sub_problems->at(spid);
        const auto& n = this->node_positions[spid];
        auto& load_vec = this->load_vector[spid];
        if(problem_type == utils::PROBLEM_TYPE_2D){
            for(auto& s:p.internal_loads){
                s->apply_load_2D(n, load_vec);
            }
        } else if(problem_type == utils::PROBLEM_TYPE_3D){
            for(auto& s:p.internal_loads){
                s->apply_load_3D(n, load_vec);
            }
        }
    }
    for(auto& s:loads){
        s.clear_curvature_data();
    }
}


void Meshing::distribute_elements(const std::vector<size_t>& geom_elem_mapping, 
                                  std::vector<std::unique_ptr<MeshElement>>& element_list){
    if(geometries.size() > 1){
        for(size_t i = 0; i < geometries.size(); ++i){
            auto& g = geometries[i];
            g->mesh.clear();
            if(i == 0){
                g->mesh.resize(geom_elem_mapping[i]);
            } else {
                g->mesh.resize(geom_elem_mapping[i] - geom_elem_mapping[i-1]);
            }
        }
        size_t geom_num = 0;
        auto& g = geometries[geom_num];
        std::move(element_list.begin(), 
                  element_list.begin() + geom_elem_mapping[geom_num], 
                  g->mesh.begin());
        ++geom_num;
        while(geom_num < geometries.size()){
            auto& g = geometries[geom_num];
            std::move(element_list.begin() + geom_elem_mapping[geom_num-1], 
                      element_list.begin() + geom_elem_mapping[geom_num], 
                      g->mesh.begin());
            ++geom_num;
        }
    } else if(geometries.size() == 1){
        auto& g = geometries[0];
        g->mesh.clear();
        g->mesh.reserve(element_list.size());
        for(auto& e:element_list){
            g->mesh.push_back(std::move(e));
        }
    }
}

void Meshing::distribute_boundary_elements(){
    if(geometries.size() > 1){
        // Dirty and slow, but it was easier to implement
        // TODO: improve this
        for(auto& b:this->boundary_elements){
            for(auto& g:this->geometries){
                bool found = false;
                for(const auto& e:g->mesh){
                    if(e.get() == b.parent){
                        g->boundary_mesh.push_back(&b);
                        found = true;
                        break;
                    }
                }
                if(found){
                    break;
                }
            }
        }
        size_t sum = 0;
        for(auto& g:this->geometries){
            sum += g->boundary_mesh.size();
        }
        logger::log_assert(sum == this->boundary_elements.size(), logger::ERROR, "Not all boundary elements were distributed: {} != {}", sum, this->boundary_elements.size());
        for(auto& b:this->inter_geometry_boundary){
            for(auto& g:this->geometries){
                bool found = false;
                for(const auto& e:g->mesh){
                    if(e.get() == b->parent){
                        g->inter_geometry_boundary_mesh.push_back(b);
                        found = true;
                        break;
                    }
                }
                if(found){
                    break;
                }
            }
        }
        sum = 0;
        for(auto& g:this->geometries){
            sum += g->inter_geometry_boundary_mesh.size();
        }
        logger::log_assert(sum == this->inter_geometry_boundary.size(), logger::ERROR, "Not all inter geoemtry boundary elements were distributed: {} != {}", sum, this->inter_geometry_boundary.size());
        for(auto& g:this->geometries){
            g->boundary_mesh.shrink_to_fit();
            g->inter_geometry_boundary_mesh.shrink_to_fit();
        }
    } else if(geometries.size() == 1){
        auto& g = geometries[0];
        g->boundary_mesh.clear();
        g->boundary_mesh.reserve(boundary_elements.size());
        for(auto& e:boundary_elements){
            g->boundary_mesh.push_back(&e);
        }
    }
}

void Meshing::generate_load_vector(const TopoDS_Shape& shape){
    (void) shape;
    const size_t N = this->elem_info->get_nodes_per_element();
    const size_t Nb = this->elem_info->get_boundary_nodes_per_element();
    size_t dof = this->elem_info->get_dof_per_node();

    if(this->elem_info->get_problem_type() == utils::PROBLEM_TYPE_2D){
        auto is_between_points = [](gp_Pnt p1, gp_Pnt p2, gp_Pnt p)->bool{
            gp_Mat M(1, p1.X(), p1.Y(), 1, p2.X(), p2.Y(), 1, p.X(), p.Y());
            bool in_line = std::abs(M.Determinant()) < Precision::Confusion();
            bool within_bounds = p.Distance(p1) - p1.Distance(p2) < Precision::Confusion() &&
                                 p.Distance(p2) - p1.Distance(p2) < Precision::Confusion();

            return in_line && within_bounds;
        };
        double F = 0;
        /**
         * This is a mess, but I didn't have too much of a choice. The clean
         * approach is the one used in the 3D, but it doesn't work well for
         * regular square meshes for some reason (it's probably due to how
         * the meshing is implemented in gmsh or something).
         *
         * So, this is what thing does:
         * 1. Gets info on each load and caches it
         * 2. Goes through the elements and checks if the nodes are completely
         *    inside it
         * 3. If they are, create local load vector and apply it to the global
         *    one
         * 4. Otherwise, check if it *partially* overlaps it. If it does, get
         *    the end point and the point inside it to do the calculations. The
         *    resulting local load vector is a bit weird, but it works. Apply
         *    it to the global vector and continue
         * 5. Otherwise, just proceed to the next element
         * 6. Do so for every load applied
         */
        for(size_t spid = 0; spid < this->sub_problems->size(); ++spid){
            const auto& p = this->sub_problems->at(spid);
            const auto& node_pos = this->node_positions[spid];
            auto& load_vec = this->load_vector[spid];
            for(auto& f : p.forces){
                gp_Vec vec(f->vec/(thickness*f->S.get_dimension()));

                gp_Dir Snormal = f->S.get_normal();
                double Ssize = f->S.get_dimension();
                gp_Pnt center = f->S.get_centroid();
                gp_Dir line_dir = Snormal.Rotated(gp_Ax1(center, gp_Dir(0,0,1)), M_PI/2);
                gp_Pnt p1 = center.Translated( 0.5*Ssize*line_dir);
                gp_Pnt p2 = center.Translated(-0.5*Ssize*line_dir);

                TopoDS_Edge line = BRepBuilderAPI_MakeEdge(p1, p2);

                for(const auto& e : this->boundary_elements){
                    std::vector<gp_Pnt> list;
                    list.reserve(Nb);
                    for(size_t i = 0; i < Nb; ++i){
                        auto n = e.nodes[i];
                        if(is_between_points(p1, p2, n->point)){
                            list.push_back(n->point);
                        }
                    }
                    std::vector<double> fe;
                    if(list.size() >= 2){
                        fe = e.parent->get_f(thickness, vec, list);
                    } else if(list.size() > 1) {
                        for(size_t i = 0; i < Nb; ++i){
                            size_t j = (i+1)%Nb;
                            auto n1 = e.nodes[i]->point;
                            auto n2 = e.nodes[j]->point;
                            if(p1.IsEqual(n1, Precision::Confusion()) || p1.IsEqual(n2, Precision::Confusion()) ||
                               p2.IsEqual(n1, Precision::Confusion()) || p2.IsEqual(n2, Precision::Confusion())){
                                continue;
                            }
                            if(is_between_points(n1, n2, p1)){
                                list.push_back(p1);
                                break;
                            } else if(is_between_points(n1, n2, p2)){
                                list.push_back(p2);
                                break;
                            }
                        }
                        fe = e.parent->get_f(thickness, vec, list);
                    }
                    if(fe.size() > 0){
                        //logger::quick_log(fe);
                        for(auto& ff:fe){
                            F += ff;
                        }
                        for(size_t i = 0; i < N; ++i){
                            for(size_t j = 0; j < dof; ++j){
                                auto n = e.parent->nodes[i];
                                const long p = node_pos[n->id*dof + j];
                                if(p >= 0){
                                    load_vec[p] += fe[i*dof+j];
                                }
                            }
                        }
                    }
                }
                logger::quick_log("Force sum ", F);
            }
        }
    } else if(this->elem_info->get_problem_type() == utils::PROBLEM_TYPE_3D){
        const auto A_tri = [](const MeshNode** const nodes)->double{
            const gp_Vec v1(nodes[0]->point, nodes[1]->point);
            const gp_Vec v2(nodes[0]->point, nodes[2]->point);

            return 0.5*v1.Crossed(v2).Magnitude();
        };
        const auto A_quad = [](const MeshNode** const nodes)->double{
            const gp_Vec v1(nodes[0]->point, nodes[1]->point);
            const gp_Vec v2(nodes[0]->point, nodes[3]->point);

            const gp_Vec v3(nodes[2]->point, nodes[1]->point);
            const gp_Vec v4(nodes[2]->point, nodes[3]->point);

            return 0.5*(v1.Crossed(v2).Magnitude() + v3.Crossed(v4).Magnitude());
        };
        for(size_t spid = 0; spid < this->sub_problems->size(); ++spid){
            double F = 0;
            const auto& p = this->sub_problems->at(spid);
            const auto& node_pos = this->node_positions[spid];
            auto& load_vec = this->load_vector[spid];
            for(auto& f : p.forces){
                if(f->vec.Magnitude() < Precision::Confusion()){
                    continue;
                }

                std::vector<bool> apply_force(this->boundary_elements.size());
                double A = 0;
                if(this->elem_info->get_shape_type() == Element::Shape::TRI){
                    #pragma omp parallel for reduction(+:A)
                    for(size_t i = 0; i < apply_force.size(); ++i){
                        const auto& e = this->boundary_elements[i];
                        apply_force[i] = f->S.is_inside(e.get_centroid(Nb));
                        if(apply_force[i]){
                            A += A_tri(e.nodes);
                        }
                    }
                } else if(this->elem_info->get_shape_type() == Element::Shape::QUAD){
                    #pragma omp parallel for reduction(+:A)
                    for(size_t i = 0; i < apply_force.size(); ++i){
                        const auto& e = this->boundary_elements[i];
                        apply_force[i] = f->S.is_inside(e.get_centroid(Nb));
                        if(apply_force[i]){
                            A += A_quad(e.nodes);
                        }
                    }
                }

                gp_Vec vec(f->vec/A);

                for(size_t i = 0; i < apply_force.size(); ++i){
                    if(apply_force[i]){
                        const auto& e = this->boundary_elements[i];
                        std::vector<gp_Pnt> points(Nb);
                        for(size_t i = 0; i < Nb; ++i){
                            points[i] = e.nodes[i]->point;
                        }
                        const auto fe = e.parent->get_f(1, vec, points);
                        //logger::quick_log(fe);
                        for(auto& ff:fe){
                            F += ff;
                        }
                        for(size_t i = 0; i < N; ++i){
                            for(size_t j = 0; j < dof; ++j){
                                const auto n = e.parent->nodes[i];
                                const long p = node_pos[n->id*dof + j];
                                if(p >= 0){
                                    load_vec[p] += fe[i*dof+j];
                                }
                            }
                        }
                    }
                }
            }
            logger::quick_log("Force sum:", F);
        }
    }
}

void Meshing::prepare_for_FEM(const TopoDS_Shape& shape,
                              const std::vector<size_t>& geom_elem_mapping,
                              const std::vector<ElementShape>& base_mesh,
                              const std::vector<Force>& forces, 
                              const std::vector<Support>& supports){

    (void) shape;
    std::vector<std::unique_ptr<MeshElement>> element_list;
    element_list.reserve(base_mesh.size());

    auto comp = [](const std::unique_ptr<MeshNode>& a, const std::unique_ptr<MeshNode>& b){ return a->id < b->id;};
    std::sort(this->node_list.begin(), this->node_list.end(), comp);

    size_t dof = this->elem_info->get_dof_per_node();

    size_t current = 0;
    for(size_t i = 0; i < this->node_list.size(); ++i){
        auto& n = this->node_list[i];

        for(size_t j = 0; j < dof; ++j){
            n->u_pos[j] = 0;
        }

        for(auto& s : supports){
            if(s.S.is_inside(n->point)){
                std::vector<bool> sup_pos = this->get_support_dof(s, this->elem_info);
                for(size_t j = 0; j < dof; ++j){
                    if(sup_pos[j]){
                        n->u_pos[j] = -1;
                    }
                }
            }
        }
        for(size_t j = 0; j < dof; ++j){
            if(n->u_pos[j] != -1){
                n->u_pos[j] = current;
                ++current;
            }
        }
    }

    this->load_vector.clear();
    this->load_vector.resize(current);

    for(auto& e : base_mesh){
        element_list.emplace_back(this->elem_info->make_element(e));
    }


    const size_t N = this->elem_info->get_nodes_per_element();
    this->inverse_mesh.clear();
    for(const auto& e : element_list){
        for(size_t i = 0; i < N; ++i){
            this->inverse_mesh.emplace(e->nodes[i]->id, e.get());
        }
    }


    if(this->elem_info->get_problem_type() == utils::PROBLEM_TYPE_2D){
        for(auto& f : forces){
            gp_Vec vec(f.vec/(thickness*f.S.get_dimension()));

            gp_Dir Snormal = f.S.get_normal();
            double Ssize = f.S.get_dimension();
            gp_Pnt center = f.S.get_centroid();
            gp_Dir line_dir = Snormal.Rotated(gp_Ax1(center, gp_Dir(0,0,1)), M_PI/2);
            gp_Pnt p1 = center.Translated( 0.5*Ssize*line_dir);
            gp_Pnt p2 = center.Translated(-0.5*Ssize*line_dir);

            TopoDS_Edge line = BRepBuilderAPI_MakeEdge(p1, p2);

            auto is_inside = [&](gp_Pnt p)->bool{
                if(p.IsEqual(center, Precision::Confusion())){
                    return true;
                } else {
                    bool in_line = std::abs(Snormal.Dot(gp_Vec(p, center))) < Precision::Confusion();
                    bool within_bounds = p.Distance(p1) - p1.Distance(p2) < Precision::Confusion() &&
                                         p.Distance(p2) - p1.Distance(p2) < Precision::Confusion();
                    return in_line && within_bounds;
                }
            };

            auto is_between_points = [](gp_Pnt p1, gp_Pnt p2, gp_Pnt p)->bool{
                gp_Mat M(1, p1.X(), p1.Y(), 1, p2.X(), p2.Y(), 1, p.X(), p.Y());
                bool in_line = std::abs(M.Determinant()) < Precision::Confusion();
                bool within_bounds = p.Distance(p1) - p1.Distance(p2) < Precision::Confusion() &&
                                     p.Distance(p2) - p1.Distance(p2) < Precision::Confusion();

                return in_line && within_bounds;
            };

            for(auto& e : element_list){
                std::vector<Node*> list;
                int last = -1;
                for(size_t i = 0; i < N; ++i){
                    auto n = e->nodes[i];
                    if(is_inside(n->point)){
                        // maintain ordering
                        if(last == -1){
                            list.push_back(n);
                            last = i;
                        } else if(i - last == 1){
                            list.push_back(n);
                        } else {
                            list.insert(list.begin(), n);
                        }
                    }
                }
                std::vector<double> fe;
                if(list.size() == 2){
                    fe = e->get_f(thickness, vec, {list[0]->point, list[1]->point});
                } else if(list.size() == 1){
                    for(size_t i = 0; i < N; ++i){
                        size_t j = (i+1)%N;
                        auto n1 = e->nodes[i]->point;
                        auto n2 = e->nodes[j]->point;
                        if(p1.IsEqual(n1, Precision::Confusion()) || p1.IsEqual(n2, Precision::Confusion()) ||
                           p2.IsEqual(n1, Precision::Confusion()) || p2.IsEqual(n2, Precision::Confusion())){
                            break;
                        }
                        if(is_between_points(n1, n2, p1)){
                            fe = e->get_f(thickness, vec, {n1, p1});
                            // logger::quick_log(n1.X(), n1.Y(), p1.X(), p1.Y(), n2.X(), n2.Y());
                            // logger::quick_log(fe);
                            break;
                        } else if(is_between_points(n1, n2, p2)){
                            fe = e->get_f(thickness, vec, {n1, p2});
                            // logger::quick_log(n1.X(), n1.Y(), p2.X(), p2.Y(), n2.X(), n2.Y());
                            // logger::quick_log(fe);
                            break;
                        }
                    }
                }
                // if(fe.size() == 0){
                //     std::vector<gp_Pnt> points = e->get_intersection_points(line);
                //     if(points.size() == 2){
                //         fe = e->get_f(dir, norm, {points[0], points[1]});
                //     }
                // }
                if(fe.size() > 0){
                    logger::quick_log(fe);
                    for(size_t i = 0; i < N; ++i){
                        for(size_t j = 0; j < dof; ++j){
                            auto n = e->nodes[i];
                            if(n->u_pos[j] >= 0){
                                this->load_vector[0][n->u_pos[j]] += fe[i*dof+j];
                            }
                        }
                    }
                }
            }
        }
    } else if(this->elem_info->get_problem_type() == utils::PROBLEM_TYPE_3D){
        for(auto& f : forces){
            gp_Vec vec(f.vec/f.S.get_dimension());

            // TODO: generalize this and make it faster
            for(auto& e : element_list){
                std::vector<gp_Pnt> points;
                points.reserve(N);
                for(size_t i = 0; i < N; ++i){
                    const auto& p = e->nodes[i]->point;
                    if(f.S.is_inside(p)){
                        points.push_back(p);
                    }
                }
                if(points.size() == 3){
                    auto fe = e->get_f(1, vec, points);
                    logger::quick_log(fe);
                    for(size_t i = 0; i < N; ++i){
                        for(size_t j = 0; j < dof; ++j){
                            auto n = e->nodes[i];
                            if(n->u_pos[j] >= 0){
                                this->load_vector[0][n->u_pos[j]] += fe[i*dof+j];
                            }
                        }
                    }
                }
            }
        }
    }

    if(geometries.size() > 1){
        for(size_t i = 0; i < geometries.size(); ++i){
            auto& g = geometries[i];
            g->mesh.clear();
            if(i == 0){
                g->mesh.resize(geom_elem_mapping[i]);
            } else {
                g->mesh.resize(geom_elem_mapping[i] - geom_elem_mapping[i-1]);
            }
        }
        size_t geom_num = 0;
        auto& g = geometries[geom_num];
        std::move(element_list.begin(), 
                  element_list.begin() + geom_elem_mapping[geom_num], 
                  g->mesh.begin());
        ++geom_num;
        while(geom_num < geometries.size()){
            auto& g = geometries[geom_num];
            std::move(element_list.begin() + geom_elem_mapping[geom_num-1], 
                      element_list.begin() + geom_elem_mapping[geom_num], 
                      g->mesh.begin());
            ++geom_num;
        }
    } else if(geometries.size() == 1){
        auto& g = geometries[0];
        g->mesh.clear();
        g->mesh.reserve(element_list.size());
        for(auto& e:element_list){
            g->mesh.push_back(std::move(e));
        }
    }
}

void Meshing::prune(const std::vector<Force>& forces, 
                    const std::vector<Support>& supports,
                    const std::vector<double>& rho, double threshold){
    auto shape = this->make_compound(this->geometries);
    std::vector<ElementShape> list;
    std::vector<size_t> geom_elem_mapping;
    geom_elem_mapping.reserve(this->geometries.size());
    { // Remove elements
        size_t N = this->elem_info->get_nodes_per_element();
        auto r = rho.begin();
        for(const auto& g:this->geometries){
            if(g->do_topopt){
                for(const auto& e:g->mesh){
                    if(*r >= threshold){
                        list.emplace_back();
                        for(size_t i = 0; i < N; ++i){
                            auto& n = e->nodes[i];
                            list.back().nodes.push_back(static_cast<MeshNode*>(n));
                        }
                    }
                    ++r;
                }
            } else {
                for(const auto& e:g->mesh){
                    list.emplace_back();
                    for(size_t i = 0; i < N; ++i){
                        auto& n = e->nodes[i];
                        list.back().nodes.push_back(static_cast<MeshNode*>(n));
                    }
                }
            }
            geom_elem_mapping.push_back(list.size());
        }
    }
    // Redo boundary mesh? Hopefully it's not necessary

    this->optimize(list);

    this->prepare_for_FEM(shape, geom_elem_mapping, list, forces, supports);
}

std::vector<bool> Meshing::get_support_dof(const Support& support, const MeshElementFactory* elem_info) const{
    size_t size = elem_info->get_dof_per_node();
    utils::ProblemType prob_type = elem_info->get_problem_type();
    std::vector<bool> pos(size);
    switch(size){
        case 6:
            pos[5] = support.MZ;
            pos[4] = support.MY;
            pos[3] = support.MX;
            [[fallthrough]];
        case 3:
            if(prob_type == utils::PROBLEM_TYPE_2D){
                pos[2] = support.MZ;
            } else {
                pos[2] = support.Z;
            }
            [[fallthrough]];
        case 2:
            pos[1] = support.Y;
            pos[0] = support.X;
    }

    return pos;
}

std::vector<double> Meshing::get_force_dof(const Force& force, const MeshElementFactory* elem_info) const{
    size_t size = elem_info->get_dof_per_node();
    utils::ProblemType prob_type = elem_info->get_problem_type();
    std::vector<double> f(size);
    switch(size){
        case 6:
            [[fallthrough]];
        case 3:
            if(prob_type == utils::PROBLEM_TYPE_3D){
                f[2] = -force.vec.Z();
            }
            [[fallthrough]];
        case 2:
            f[1] = -force.vec.Y();
            f[0] = -force.vec.X();
    }

    return f;
}


void Meshing::reverse_cuthill_mckee(const std::vector<ElementShape>& elem_list){
    // Create adjacency "matrix"
    std::vector<std::set<size_t>> adjacents(this->node_list.size());
    for(auto& e:elem_list){
        for(size_t i = 0; i < e.nodes.size(); ++i){
            for(size_t j = 1; j < e.nodes.size(); ++j){
                adjacents[e.nodes[i]->id].insert(e.nodes[j]->id);
                adjacents[e.nodes[j]->id].insert(e.nodes[i]->id);
            }
        }
    }
    std::vector<bool> added(this->node_list.size(), false);

    // Get node with least degree
    size_t min_degree = adjacents[0].size();
    size_t min_node = 0;
    for(size_t i = 1; i < adjacents.size(); ++i){
        if(adjacents[i].size() < min_degree){
            min_node = i;
            min_degree = adjacents[i].size();
        }
    }

    // Generate Cuthill-McKee
    
    auto comp = [&](size_t n1, size_t n2){
        return adjacents[n1].size() < adjacents[n2].size();
    };
    std::queue<size_t> queue;
    queue.push(min_node);
    added[min_node] = true;
    std::vector<size_t> result;
    result.reserve(this->node_list.size());
    if(this->proj_data->contact_data.contact_type == FiniteElement::ContactType::RIGID){
        while(!queue.empty()){
            size_t node = queue.front();
            auto& adj = adjacents[node];

            std::vector<size_t> new_nodes;
            new_nodes.reserve(adj.size());
            for(auto& n:adj){
                if(!added[n]){
                    added[n] = true;
                    auto upper = std::upper_bound(new_nodes.begin(), new_nodes.end(), n, comp);
                    new_nodes.insert(upper, n);
                }
            }
            for(auto& n:new_nodes){
                queue.push(n);
            }

            result.push_back(node);
            queue.pop();
        }
        logger::log_assert(result.size() == this->node_list.size(), logger::ERROR, "Mesh contains disconnected nodes. result: {}, node_list: {}.", result.size(), node_list.size());
    } else {
        size_t loop_count = 0;
        do{
            bool all_added = true;
            while(!queue.empty()){
                size_t node = queue.front();
                auto& adj = adjacents[node];

                std::vector<size_t> new_nodes;
                new_nodes.reserve(adj.size());
                for(auto& n:adj){
                    if(!added[n]){
                        added[n] = true;
                        auto upper = std::upper_bound(new_nodes.begin(), new_nodes.end(), n, comp);
                        new_nodes.insert(upper, n);
                    }
                }
                for(auto& n:new_nodes){
                    queue.push(n);
                }

                result.push_back(node);
                queue.pop();
            }
            for(size_t i = 0; i < added.size(); ++i){
                if(!added[i]){
                    queue.push(i);
                    added[i] = true;
                    all_added = false;
                    break;
                }
            }
            ++loop_count;
            if(all_added){
                break;
            }
        } while(result.size() != this->node_list.size());
        logger::log_assert(loop_count == this->geometries.size(), logger::WARNING, "Number of separate meshes is different from number of geometries: {} and {}. Disconnected nodes?.", loop_count, geometries.size());
        logger::log_assert(result.size() == this->node_list.size(), logger::ERROR, "Inconsistency in the number of nodes after RCM. result: {}, node_list: {}.", result.size(), node_list.size());
    }

    // Reorder node list
    std::vector<std::unique_ptr<MeshNode>> new_node_list(this->node_list.size());
    size_t pos = 0;
    for(auto i = result.rbegin(); i < result.rend(); ++i){
        new_node_list[pos] = std::move(this->node_list[*i]);
        new_node_list[pos]->id = pos;
        ++pos;
    }

    this->node_list = std::move(new_node_list);
}

bool Meshing::is_strictly_inside2D(gp_Pnt p, TopoDS_Shape s) const{
    BRepClass3d_SolidClassifier insider(s);
    insider.Perform(p, Precision::Confusion());
    bool on_edge = false;
    TopoDS_Vertex v = BRepBuilderAPI_MakeVertex(p);
    for(TopExp_Explorer exp(s, TopAbs_EDGE); exp.More(); exp.Next()){
        BRepExtrema_DistShapeShape edge_checker(v, exp.Current());
        edge_checker.Perform();
        double dist = edge_checker.Value();
        if(dist < Precision::Confusion()){
            on_edge = true;
            break;
        }
    }
    return insider.State() == TopAbs_ON && !on_edge;
}

bool Meshing::is_strictly_inside3D(gp_Pnt p, TopoDS_Shape s) const{
    BRepClass3d_SolidClassifier insider(s);
    insider.Perform(p, Precision::Confusion());
    bool on_bounds = false;
    TopoDS_Vertex v = BRepBuilderAPI_MakeVertex(p);
    for(TopExp_Explorer exp(s, TopAbs_FACE); exp.More(); exp.Next()){
        BRepClass3d_SolidClassifier insider(exp.Current());
        insider.Perform(p, Precision::Confusion());
        if(insider.State() == TopAbs_ON){
            on_bounds = true;
            break;
        }
    }
    return insider.State() == TopAbs_IN && !on_bounds;
}

void Meshing::optimize(std::vector<ElementShape>& list, std::unordered_map<size_t, MeshNode*>* id_map){
    // Prune unused nodes
    size_t new_id = 0;
    for(auto& n : this->node_list){
        n->id = new_id;
        ++new_id;
    }

    std::vector<bool> used(this->node_list.size(), false);
    for(size_t i = 0; i < list.size(); ++i){
        for(auto& n:list[i].nodes){
            used[n->id] = true;
        }   
    }

    bool renumber = false;
    std::vector<size_t> ndel;
    std::vector<size_t> bndel;
    auto it = this->node_list.begin();
    while(it < this->node_list.end()){
        if(!used[(*it)->id]){
            ndel.push_back(it - this->node_list.begin());
            renumber = true;
            if(id_map != nullptr){
                auto mit = id_map->begin();
                while(mit != id_map->end()){
                    if(mit->second == it->get()){
                        mit = id_map->erase(mit);
                    } else {
                        ++mit;
                    }
                }
            }
        }
        ++it;
    }
    auto itb = this->boundary_node_list.begin();
    while(itb < this->boundary_node_list.end()){
        if(!used[(*itb)->id]){
            bndel.push_back(itb - this->boundary_node_list.begin());
            renumber = true;
        }
        ++itb;
    }
    it = this->node_list.begin();
    size_t offset = 0;
    for(auto n:ndel){
        this->node_list.erase(it+n - offset);
        ++offset;
    }
    itb = this->boundary_node_list.begin();
    offset = 0;
    for(auto n:bndel){
        this->boundary_node_list.erase(itb+n - offset);
        ++offset;
    }
    used.clear();

    if(renumber){
        new_id = 0;
        for(auto& n : this->node_list){
            n->id = new_id;
            ++new_id;
        }
    }

    this->reverse_cuthill_mckee(list);
}

std::vector<ElementShape> Meshing::generate_element_shapes(
            const std::vector<size_t>& elem_node_tags, 
            size_t nodes_per_elem,
            const std::unordered_map<size_t, MeshNode*>& id_map,
            bool calc_normals){
    std::vector<size_t> filtered_tags;
    filtered_tags.reserve(elem_node_tags.size());
    std::vector<size_t> elem_tmp(nodes_per_elem);
    const size_t N0 = elem_node_tags.size()/nodes_per_elem;
    for(size_t i = 0; i < N0; ++i){
        bool found_all_nodes = true;
        bool all_nodes_equal = true;
        for(size_t j = 0; j < nodes_per_elem; ++j){
            size_t n = elem_node_tags[i*nodes_per_elem + j];
            auto pos = id_map.find(n);
            if(pos != id_map.end()){
                elem_tmp[j] = n;
            } else {
                found_all_nodes = false;
                break;
            }
            if(j > 0){
                all_nodes_equal &= id_map.at(n) == id_map.at(elem_tmp[j-1]);
            }
        }
        if(found_all_nodes && !all_nodes_equal){
            filtered_tags.insert(filtered_tags.end(), elem_tmp.begin(), elem_tmp.end());
        }
    }
    const size_t N = filtered_tags.size()/nodes_per_elem;
    std::vector<ElementShape> list(N);
    std::vector<gp_Dir> normals(N);
    if(calc_normals){
        if(this->elem_info->get_problem_type() == utils::PROBLEM_TYPE_2D){
            for(size_t i = 0; i < N; ++i){
                gp_Pnt p0 = id_map.at(filtered_tags[i*nodes_per_elem + 0])->point;
                gp_Pnt p1 = id_map.at(filtered_tags[i*nodes_per_elem + 1])->point;
                gp_Dir n(gp_Vec(p0, p1));
                n.Cross(gp_Dir(0,0,1));
                normals[i] = std::move(n);
            }
        } else {
            for(size_t i = 0; i < N; ++i){
                gp_Pnt p0 = id_map.at(filtered_tags[i*nodes_per_elem + 0])->point;
                gp_Pnt p1 = id_map.at(filtered_tags[i*nodes_per_elem + 1])->point;
                gp_Pnt p2 = id_map.at(filtered_tags[i*nodes_per_elem + 2])->point;
                gp_Dir n0(gp_Vec(p1, p0));
                gp_Dir n1(gp_Vec(p1, p2));
                normals[i] = n0.Crossed(n1);
            }
        }
    }
    #pragma omp parallel for
    for(size_t j = 0; j < N; ++j){
        list[j].nodes.resize(nodes_per_elem);
        for(size_t i = 0; i < nodes_per_elem; ++i){
            size_t n = filtered_tags[j*nodes_per_elem + i];

            MeshNode* node = id_map.at(n);
            list[j].nodes[i] = node;
        }
        list[j].normal = normals[j];
    }

    return list;
}

void Meshing::deduplicate(std::unordered_map<size_t, MeshNode*>& id_map, const std::unordered_map<size_t, size_t>& duplicate_map, const bool delete_dups){
    for(const auto& dn:duplicate_map){
        auto n = this->boundary_node_list.begin();
        while(n < this->boundary_node_list.end()){
            if((*n)->id == dn.first){
                break;
            }
            ++n;
        }
        logger::log_assert(n != this->boundary_node_list.end(),
                logger::ERROR,
                "WHAT");
        id_map[dn.second] = *n;
        if(delete_dups){
            this->boundary_node_list.erase(n);

            auto nn = this->node_list.begin();
            while(nn < this->node_list.end()){
                if((*nn)->id == dn.second){
                    break;
                }
                ++nn;
            }
            this->node_list.erase(nn);
        }
    }
}

std::vector<size_t> Meshing::find_duplicates(std::unordered_map<size_t, MeshNode*>& id_map, const bool delete_dups){
    //const auto point_sort = 
    //    [](const std::unique_ptr<MeshNode>& n1, const std::unique_ptr<MeshNode>& n2)->bool{
    //        if(n1->point.X() < n2->point.X()){
    //            return true;
    //        } else if(n1->point.X() == n2->point.X() && n1->point.Y() < n2->point.Y()){
    //            return true;
    //        } else if(n1->point.X() == n2->point.X() && n1->point.Y() == n2->point.Y() && n1->point.Z() < n2->point.Z()){
    //            return true;
    //        }
    //        return false;
    //    };

    std::vector<size_t> duplicates;
    this->number_duplicated_nodes = 0;

    // Assume there are only duplicates on the boundary (as it should be)
    for(auto i = this->boundary_node_list.begin(); i < this->boundary_node_list.end(); ++i){
        if(!delete_dups){
            auto dup_it = std::find(duplicates.begin(), duplicates.end(), (*i)->id);
            if(dup_it != duplicates.end()){
                continue;
            }
        }
        auto j = i + 1;
        while(j < this->boundary_node_list.end()){
            if((*i)->point.IsEqual((*j)->point, Precision::Confusion())){

                id_map[(*j)->id] = id_map[(*i)->id];
                ++this->number_duplicated_nodes;
                //gp_Pnt test = (*i)->point;
                if(delete_dups){
                    auto n = this->node_list.begin();
                    while(n < this->node_list.end()){
                        if(n->get() == *j){
                            break;
                        }
                        ++n;
                    }
                    logger::log_assert(n != this->node_list.end(), logger::ERROR, "boundary node not present in main node list");
                    this->node_list.erase(n);
                    //const size_t diff = i - this->boundary_node_list.begin();
                    j = this->boundary_node_list.erase(j);
                    //i = this->boundary_node_list.begin() + diff;
                } else {
                    duplicates.push_back((*j)->id);
                    ++j;
                }

                // There may be more than one duplicate, e.g., multiple geometry
                // intersection, do don't break
            } else {
                ++j;
            }
        }
        ++i;
    }

    // // The usual method to find duplicates is to test each one against the
    // // other, which is O(N^2). However, std::sort() is Nlog2N, then,
    // // figuring out the duplicates becomes faster, as you can limit search
    // // to nodes with equal X. It should be O(M^2) with M << N.
    // std::sort(this->node_list.begin(), this->node_list.end(), point_sort);

    // double eps = Precision::Confusion();
    // logger::quick_log("Original number of nodes:", node_list.size());
    // auto i = this->node_list.begin();
    // while(i < this->node_list.end()-1){
    //     auto j = i + 1;
    //     if(std::abs((*i)->point.X()) < 10 &&
    //        std::abs((*i)->point.Y()) < 10 &&
    //        std::abs((*i)->point.Z()) < 10){
    //         logger::quick_log((*i)->point.X(), (*i)->point.Y(), (*i)->point.Z());
    //     }
    //     //while(j < this->node_list.end()){
    //     while(j < this->node_list.end() &&
    //           (*i)->point.X() == (*j)->point.X() &&
    //           (*i)->point.Y() == (*j)->point.Y() &&
    //           (*i)->point.Z() == (*j)->point.Z()){

    //         id_map[(*j)->id] = id_map[(*i)->id];
    //         ++this->number_duplicated_nodes;
    //         //gp_Pnt test = (*i)->point;
    //         if(delete_dups){
    //             auto bn = std::find(this->boundary_node_list.begin(), this->boundary_node_list.end(), j->get());
    //             if(bn != this->boundary_node_list.end()){
    //                 this->boundary_node_list.erase(bn);
    //             }
    //             const size_t diff = i - this->node_list.begin();
    //             j = this->node_list.erase(j);
    //             i = this->node_list.begin() + diff;
    //         } else {
    //             ++j;
    //         }
    //     }
    //     if(delete_dups){
    //         ++i;
    //     } else {
    //         i += j - i;
    //     }
    // }
    logger::quick_log("Final number of nodes:", node_list.size());
    logger::quick_log("Duplicates:", this->number_duplicated_nodes);
    
    return duplicates;
}

TopoDS_Shape Meshing::make_compound(const std::vector<Geometry*>& geometries, const double scale) const{
    // Assuming linear with a single element type
    BRep_Builder builder;
    TopoDS_CompSolid comp;
    builder.MakeCompSolid(comp);
    if(scale == 1.0){
        for(const auto& geom:geometries){
            builder.Add(comp, geom->shape);
        }
    } else {
        for(const auto& geom:geometries){
            gp_Trsf t;
            t.SetScale(gp_Pnt(0, 0, 0), scale);
            BRepBuilderAPI_Transform transf(geom->shape, t, true);
            builder.Add(comp, transf.Shape());
        }
    }
    return comp;
}
bool Meshing::adapt_for_boundary_condition_inside(TopoDS_Shape& shape, const std::vector<Force>& forces, const std::vector<Support>& supports){
    bool has_condition_inside = false;
    // TODO: adapt for 3D, maybe (methods to insert loads/supports into the
    // geometry if needed).
    TopoDS_Shape sh = BRepBuilderAPI_Copy(shape);
    for(auto& f:forces){
        if(this->is_strictly_inside2D(f.S.get_centroid(), shape)){
            has_condition_inside = true;

            BOPAlgo_Splitter splitter;
            splitter.SetNonDestructive(true);
            splitter.AddArgument(shape);
            splitter.AddTool(f.S.get_shape());
            splitter.Perform();
            sh = splitter.Shape();
        }
    }
    for(auto& s:supports){
        if(this->is_strictly_inside2D(s.S.get_centroid(), shape)){
            has_condition_inside = true;

            BOPAlgo_Splitter splitter;
            splitter.SetNonDestructive(true);
            splitter.AddArgument(shape);
            splitter.AddTool(s.S.get_shape());
            splitter.Perform();
            sh = splitter.Shape();
        }
    }
    shape = sh;

    return has_condition_inside;
}

void Meshing::fix_node_point_precision(){
    const std::function<bool(const std::unique_ptr<MeshNode>&,const std::unique_ptr<MeshNode>&)> coord_sort[3] ={
        [](const std::unique_ptr<MeshNode>& n1, const std::unique_ptr<MeshNode>& n2)->bool{
            if(n1->point.X() < n2->point.X()){
                return true;
            }
            return false;
        },
        [](const std::unique_ptr<MeshNode>& n1, const std::unique_ptr<MeshNode>& n2)->bool{
            if(n1->point.Y() < n2->point.Y()){
                return true;
            }
            return false;
        },
        [](const std::unique_ptr<MeshNode>& n1, const std::unique_ptr<MeshNode>& n2)->bool{
            if(n1->point.Z() < n2->point.Z()){
                return true;
            }
            return false;
        }
    };

    const double eps = Precision::Confusion();

    for(size_t i = 0; i < 3; ++i){
        std::sort(this->node_list.begin(), this->node_list.end(), coord_sort[i]);
        const size_t C = 1+i;
        double curr_val = this->node_list[0]->point.Coord(C);
        for(size_t it = 1; it < this->node_list.size(); ++it){
            if(std::abs(this->node_list[it]->point.Coord(C) - curr_val) < eps){
                this->node_list[it]->point.SetCoord(C, curr_val);
            } else {
                curr_val = this->node_list[it]->point.Coord(C);
            }
        }
    }

}


void Meshing::extend_vector(const size_t subproblem, const std::vector<double>& v, std::vector<double>& v_ext) const{
    const size_t dof = this->elem_info->get_dof_per_node();
    if(v_ext.size() < this->node_list.size()*dof){
        v_ext.resize(this->node_list.size()*dof, 0);
    }
    const auto& npos = this->node_positions[subproblem];
    for(size_t i = 0; i < this->node_list.size(); ++i){
        const auto& n = this->node_list[i];
        for(size_t j = 0; j < dof; ++j){
            const long p = npos[dof*n->id + j];
            if(p > -1){
                v_ext[n->u_pos[j]] = v[p];
            } else {
                v_ext[n->u_pos[j]] = 0;
            }
        }
    }
}

void Meshing::de_extend_vector(const size_t subproblem, const std::vector<double>& v_ext, std::vector<double>& v) const{
    const size_t dof = this->elem_info->get_dof_per_node();
    if(v.size() < this->load_vector[subproblem].size()){
        v.resize(this->load_vector[subproblem].size(), 0);
    }
    const auto& npos = this->node_positions[subproblem];
    for(size_t i = 0; i < this->node_list.size(); ++i){
        const auto& n = this->node_list[i];
        for(size_t j = 0; j < dof; ++j){
            const long p = npos[dof*n->id + j];
            if(p > -1){
                v[p] = v_ext[n->u_pos[j]];
            }
        }
    }
}

void Meshing::distribute_node_pointers(){
    const size_t bnodes_per_elem = this->elem_info->get_boundary_nodes_per_element();
    for(auto& g:this->geometries){
        std::set<Node*> nodes;
        for(auto& e:g->mesh){
            for(size_t i = 0; i < bnodes_per_elem; ++i){
                nodes.insert(e->nodes[i]);
            }
        }
        g->node_list.insert(g->node_list.begin(), nodes.begin(), nodes.end());
        std::set<const MeshNode*> bnodes;
        for(auto& e:g->boundary_mesh){
            for(size_t i = 0; i < bnodes_per_elem; ++i){
                bnodes.insert(e->nodes[i]);
            }
        }
        g->boundary_node_list.insert(g->boundary_node_list.begin(), bnodes.begin(), bnodes.end());
    }
}

void Meshing::generate_initial_u_contact(std::vector<double>& u) const{
    const size_t bnum = this->elem_info->get_boundary_nodes_per_element();
    const size_t dof = this->elem_info->get_dof_per_node();

    //const double EPS = 1e-7;
    //for(const auto& e:this->paired_boundary){
    //    for(size_t i = 0; i < bnum; ++i){
    //        const auto n1 = e.b1->nodes[i];
    //        const auto n2 = e.b2->nodes[i];
    //        const auto normal = e.b1->normal;
    //        for(size_t j = 0; j < dof; ++j){
    //            auto ni1 = node_positions[0][n1->u_pos[j]];
    //            auto ni2 = node_positions[0][n2->u_pos[j]];
    //            if(ni1 > -1){
    //                u[ni1] -= normal.Coord(1+j)*EPS;
    //            }
    //            if(ni2 > -1){
    //                u[ni2] += normal.Coord(1+j)*EPS;
    //            }
    //        }
    //    }
    //}
    /*
    utils::ProblemType prob_type = this->elem_info->get_problem_type();
    size_t dim = 3;
    if(prob_type == utils::PROBLEM_TYPE_2D){
        dim = 2;
    } else if(prob_type == utils::PROBLEM_TYPE_3D){
        dim = 3;
    }
    for(const auto& g:this->geometries){
        std::vector<double> F(dof, 0);
        bool supported = false;
        for(const auto& e:g->boundary_node_list){
            for(size_t i = 0; i < dof; ++i){
                if(this->node_positions[0][e->u_pos[i]] > -1){
                    F[i] += this->global_load_vector[e->u_pos[i]];
                } else {
                    supported = true;
                    break;
                }
            }
            if(supported){
                break;
            }
        }
        if(supported){
            continue;
        }
        double Fnorm = 0;
        for(const auto& fi:F){
            Fnorm += fi*fi;
        }
        Fnorm = std::sqrt(Fnorm);
        for(auto& fi:F){
            fi *= 1e-7/Fnorm;
        }
        for(const auto& e:g->node_list){
            for(size_t i = 0; i < dof; ++i){
                const size_t ni = this->node_positions[0][e->u_pos[i]];
                u[ni] = F[i];
            }
            if(supported){
                break;
            }
        }
    }
    */
}
