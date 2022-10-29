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


void Meshing::generate_elements(const TopoDS_Shape& shape,
                                const std::vector<size_t>& geom_elem_mapping, 
                                const std::vector<size_t>& elem_node_tags, 
                                const std::vector<size_t>& bound_elem_node_tags,
                                std::unordered_map<size_t, MeshNode*>& id_map,
                                const std::vector<Force>& forces, 
                                const std::vector<Support>& supports,
                                const bool deduplicate,
                                const bool prune){

    logger::quick_log("Generating elements and preparing for finite element analysis...");

    if(deduplicate){
       this->find_duplicates(id_map);
    }

    const size_t nodes_per_elem = this->elem_info->get_nodes_per_element();
    const size_t bound_nodes_per_elem = this->elem_info->get_boundary_nodes_per_element();

    auto list = this->generate_element_shapes(elem_node_tags, nodes_per_elem, id_map);
    this->optimize(list, prune);
    auto elements = this->create_element_list(list, this->elem_info);
    list.clear();
    populate_inverse_mesh(elements);

    {
        auto bound_list = this->generate_element_shapes(bound_elem_node_tags, bound_nodes_per_elem, id_map);
        this->populate_boundary_elements(bound_list);
    }

    this->apply_supports(supports);

    this->generate_load_vector(shape, forces);

    this->distribute_elements(geom_elem_mapping, elements);
    elements.clear();

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

void Meshing::populate_boundary_elements(const std::vector<ElementShape>& boundary_base_mesh){
    const size_t N = this->elem_info->get_boundary_nodes_per_element();
    this->boundary_elements.clear();
    this->boundary_elements.reserve(boundary_base_mesh.size());
    for(size_t i = 0; i < boundary_base_mesh.size(); ++i){
        const auto& b = boundary_base_mesh[i];
        std::set<MeshElement*> common_nodes;
        auto eq_range = this->inverse_mesh.equal_range(b.nodes[0]->id);
        for(auto k = eq_range.first; k != eq_range.second; ++k){
            common_nodes.insert(k->second);
        }
        for(size_t j = 1; j < N; ++j){
            std::set<MeshElement*> tmp_common_nodes;
            std::set<MeshElement*> tmp_comp;
            eq_range = this->inverse_mesh.equal_range(b.nodes[j]->id);
            for(auto k = eq_range.first; k != eq_range.second; ++k){
                tmp_comp.insert(k->second);
            }
            std::set_intersection(common_nodes.begin(), common_nodes.end(),
                                  tmp_comp.begin(), tmp_comp.end(),
                                  std::inserter(tmp_common_nodes, tmp_common_nodes.begin()));
            common_nodes = std::move(tmp_common_nodes);
        }
        logger::log_assert(common_nodes.size() == 1, logger::WARNING, "boundary element associated to more than one element.");
        this->boundary_elements.emplace_back(b.nodes, *common_nodes.begin());
    }
}

void Meshing::apply_supports(const std::vector<Support>& supports){
    size_t dof = this->elem_info->get_dof_per_node();

    #pragma omp parallel for
    for(size_t i = 0; i < this->boundary_node_list.size(); ++i){
        auto& n = this->boundary_node_list[i];

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
    }
    size_t current = 0;
    for(size_t i = 0; i < this->node_list.size(); ++i){
        auto& n = this->node_list[i];
        for(size_t j = 0; j < dof; ++j){
            if(n->u_pos[j] > -1){
                n->u_pos[j] = current;
                ++current;
            }
        }
    }
    this->load_vector.clear();
    this->load_vector.resize(current);
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

void Meshing::generate_load_vector(const TopoDS_Shape& shape,
                                   const std::vector<Force>& forces){
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
        for(auto& f : forces){
            double norm = f.vec.Magnitude()/(thickness*f.S.get_dimension());
            gp_Dir dir(f.vec);

            if(this->is_strictly_inside2D(f.S.get_centroid(), shape)){
                norm = norm/2;
            }

            gp_Dir Snormal = f.S.get_normal();
            double Ssize = f.S.get_dimension();
            gp_Pnt center = f.S.get_centroid();
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
                    fe = e.parent->get_f(thickness, dir, norm, list);
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
                    fe = e.parent->get_f(thickness, dir, norm, list);
                }
                if(fe.size() > 0){
                    logger::quick_log(fe);
                    for(size_t i = 0; i < N; ++i){
                        for(size_t j = 0; j < dof; ++j){
                            auto n = e.parent->nodes[i];
                            if(n->u_pos[j] >= 0){
                                this->load_vector[n->u_pos[j]] += fe[i*dof+j];
                            }
                        }
                    }
                }
            }
        }
    } else if(this->elem_info->get_problem_type() == utils::PROBLEM_TYPE_3D){
        for(auto& f : forces){
            double norm = f.vec.Magnitude()/f.S.get_dimension();
            gp_Dir dir(f.vec);

            for(const auto& e : this->boundary_elements){
                if(f.S.is_inside(e.get_centroid(Nb))){
                    std::vector<gp_Pnt> points(Nb);
                    for(size_t i = 0; i < Nb; ++i){
                        points[i] = e.nodes[i]->point;
                    }
                    const auto fe = e.parent->get_f(1, dir, norm, points);
                    logger::quick_log(fe);
                    for(size_t i = 0; i < N; ++i){
                        for(size_t j = 0; j < dof; ++j){
                            const auto n = e.parent->nodes[i];
                            if(n->u_pos[j] >= 0){
                                this->load_vector[n->u_pos[j]] += fe[i*dof+j];
                            }
                        }
                    }
                }
            }
        }
    }

}

void Meshing::prepare_for_FEM(const TopoDS_Shape& shape,
                              const std::vector<size_t>& geom_elem_mapping,
                              const std::vector<ElementShape>& base_mesh,
                              const std::vector<Force>& forces, 
                              const std::vector<Support>& supports){

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
            double norm = f.vec.Magnitude()/(thickness*f.S.get_dimension());
            gp_Dir dir(f.vec);

            if(this->is_strictly_inside2D(f.S.get_centroid(), shape)){
                norm = norm/2;
            }

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
                    fe = e->get_f(thickness, dir, norm, {list[0]->point, list[1]->point});
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
                            fe = e->get_f(thickness, dir, norm, {n1, p1});
                            // logger::quick_log(n1.X(), n1.Y(), p1.X(), p1.Y(), n2.X(), n2.Y());
                            // logger::quick_log(fe);
                            break;
                        } else if(is_between_points(n1, n2, p2)){
                            fe = e->get_f(thickness, dir, norm, {n1, p2});
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
                                this->load_vector[n->u_pos[j]] += fe[i*dof+j];
                            }
                        }
                    }
                }
            }
        }
    } else if(this->elem_info->get_problem_type() == utils::PROBLEM_TYPE_3D){
        for(auto& f : forces){
            double norm = f.vec.Magnitude()/f.S.get_dimension();
            gp_Dir dir(f.vec);

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
                    auto fe = e->get_f(1, dir, norm, points);
                    logger::quick_log(fe);
                    for(size_t i = 0; i < N; ++i){
                        for(size_t j = 0; j < dof; ++j){
                            auto n = e->nodes[i];
                            if(n->u_pos[j] >= 0){
                                this->load_vector[n->u_pos[j]] += fe[i*dof+j];
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
            geom_elem_mapping.push_back(list.size());
        }
    }
    // Redo boundary mesh? Hopefully it's not necessary

    this->optimize(list, true);

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
                size_t k = (i+j) % e.nodes.size();
                adjacents[e.nodes[i]->id].insert(e.nodes[k]->id);
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
    logger::log_assert(result.size() == this->node_list.size(), logger::ERROR, "Mesh contains disconnected nodes.");

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

void Meshing::optimize(std::vector<ElementShape>& list, const bool prune){
    // Prune unused nodes
    size_t new_id = 0;
    for(auto& n : this->node_list){
        n->id = new_id;
        ++new_id;
    }

    if(prune){
        std::vector<bool> used(this->node_list.size(), false);
        for(size_t i = 0; i < list.size(); ++i){
            for(auto& n:list[i].nodes){
                used[n->id] = true;
            }   
        }

        bool renumber = false;
        auto it = this->node_list.begin();
        while(it < this->node_list.end()){
            if(!used[(*it)->id]){
                it = this->node_list.erase(it);
                renumber = true;
            } else {
                ++it;
            }
        }
        used.clear();

        if(renumber){
            new_id = 0;
            for(auto& n : this->node_list){
                n->id = new_id;
                ++new_id;
            }
        }
    }

    this->reverse_cuthill_mckee(list);
}

std::vector<ElementShape> Meshing::generate_element_shapes(
            const std::vector<size_t>& elem_node_tags, 
            size_t nodes_per_elem,
            const std::unordered_map<size_t, MeshNode*>& id_map){
    const size_t N = elem_node_tags.size()/nodes_per_elem;
    std::vector<ElementShape> list(N);
    #pragma omp parallel for
    for(size_t j = 0; j < N; ++j){
        list[j].nodes.resize(nodes_per_elem);
        for(size_t i = 0; i < nodes_per_elem; ++i){
            size_t n = elem_node_tags[j*nodes_per_elem + i];

            MeshNode* node = id_map.at(n);
            list[j].nodes[i] = node;
        }
    }

    return list;
}

std::vector<size_t> Meshing::find_duplicates(std::unordered_map<size_t, MeshNode*>& id_map){
    const auto point_sort = 
        [](const std::unique_ptr<MeshNode>& n1, const std::unique_ptr<MeshNode>& n2)->bool{
            if(n1->point.X() < n2->point.X()){
                return true;
            } else if(n1->point.Y() < n2->point.Y()){
                return true;
            } else if(n1->point.Z() < n2->point.Z()){
                return true;
            }
            return false;
        };

    std::vector<size_t> duplicates;
    // The usual method to find duplicates is to test each one against the
    // other, which is O(N^2). However, std::sort() is Nlog2N, then,
    // figuring out the duplicates is O(N), so it's a much better
    // alternative.
    std::sort(this->node_list.begin(), this->node_list.end(), point_sort);

    auto i = this->node_list.begin();
    while(i < this->node_list.end()-1){
        if((*i)->point.IsEqual((*(i+1))->point, Precision::Confusion())){
            id_map[(*(i+1))->id] = id_map[(*i)->id];
            this->node_list.erase(i+1);
        } else {
            ++i;
        }
    }
    
    return duplicates;
}

TopoDS_Shape Meshing::make_compound(const std::vector<Geometry*>& geometries) const{
    // Assuming linear with a single element type
    BRep_Builder builder;
    TopoDS_Compound comp;
    builder.MakeCompound(comp);
    for(const auto& geom:geometries){
        builder.Add(comp,geom->shape);
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
