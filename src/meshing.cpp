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


void Meshing::prepare_for_FEM(const TopoDS_Shape& shape,
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
        if(n->u_pos != nullptr){
            delete[] n->u_pos;
        }
        n->u_pos = new long[dof]();
        bool supported = false;
        for(auto& s : supports){
            if(s.S.is_inside(n->point)){
                size_t offset = 0;
                std::vector<long> sup_pos = this->get_support_dof(offset, 0, s, this->elem_info);
                for(size_t j = 0; j < dof; ++j){
                    if(sup_pos[j] >= 0){
                        if(n->u_pos[j] == 0){
                            n->u_pos[j] = current;
                            ++current;
                        }
                    } else {
                        n->u_pos[j] = -1;
                    }
                }
                supported = true;
            }
        }
        if(!supported){
            for(size_t j = 0; j < dof; ++j){
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
                size_t N = e->nodes.size();
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
        // TODO
    }

    if(geometries.size() > 1){
        // Slower, but handles memory much better, especially considering it
        // may involve large amounts of memory.
        // Fortunately, it should be done only once.
        std::vector<size_t> mesh_size_per_geom(geometries.size(), 0);
        for(auto& e:element_list){
            for(size_t i = 0; i < geometries.size(); ++i){
                auto& g = geometries[i];
                if(g->is_inside(e->get_centroid())){
                    ++mesh_size_per_geom[i];
                    break;
                }
            }
        }
        for(size_t i = 0; i < geometries.size(); ++i){
            auto& g = geometries[i];
            g->mesh.clear();
            g->mesh.reserve(mesh_size_per_geom[i]);
        }
        for(auto& e:element_list){
            for(auto& g:geometries){
                if(g->is_inside(e->get_centroid())){
                    g->mesh.push_back(std::move(e));
                    break;
                }
            }
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
    { // Remove elements
        auto r = rho.begin();
        for(const auto& g:this->geometries){
            for(const auto& e:g->mesh){
                if(*r >= threshold){
                    list.emplace_back();
                    for(const auto& n:e->nodes){
                        list.back().nodes.push_back(static_cast<MeshNode*>(n));
                    }
                }
                ++r;
            }
        }
    }

    this->prune(list);

    this->prepare_for_FEM(shape, list, forces, supports);
}

std::vector<long> Meshing::get_support_dof(size_t& offset, size_t id, const Support& support, const MeshElementFactory* elem_info) const{
    size_t size = elem_info->get_dof_per_node();
    utils::ProblemType prob_type = elem_info->get_problem_type();
    std::vector<long> pos(size);
    id *= size;
    switch(size){
        case 6:
            pos[5] = support.MZ ? -1 : (id + offset++);
            pos[4] = support.MY ? -1 : (id + offset++);
            pos[3] = support.MX ? -1 : (id + offset++);
            [[fallthrough]];
        case 3:
            if(prob_type == utils::PROBLEM_TYPE_2D){
                pos[2] = support.MZ ? -1 : (id + offset++);
            } else {
                pos[2] = support.Z ? -1 : (id + offset++);
            }
            [[fallthrough]];
        case 2:
            pos[1] = support.Y ? -1 : (id + offset++);
            pos[0] = support.X ? -1 : (id + offset++);
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
    std::vector<bool> added(elem_list.size(), false);

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

void Meshing::prune(const std::vector<ElementShape>& list){
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
    this->reverse_cuthill_mckee(list);
}

std::vector<ElementShape> Meshing::generate_element_shapes(const std::vector<size_t>& elem_tags, const std::vector<size_t>& elem_node_tags, size_t nodes_per_elem, const std::unordered_map<size_t, size_t>& duplicate_map){
    std::vector<ElementShape> list;
    list.reserve(elem_tags.size());
    list.emplace_back();
    size_t i = 0;
    for(auto n:elem_node_tags){
        // gp_Pnt p(nodeCoords[n*3], nodeCoords[n*3+1], nodeCoords[n*3+2]);
        // auto get_id = [&p](const std::unique_ptr<MeshNode>& m)->bool{ return p.IsEqual(m->point, Precision::Confusion()); };
        
        // If there are duplicates, redirect the node tag to the correct,
        // deduplicated node.
        auto map_find(duplicate_map.find(n));
        if(map_find != duplicate_map.end()){
            n = map_find->second;
        }

        // Find node with id `n`
        auto get_id = [n](const std::unique_ptr<MeshNode>& m)->bool{ return n == m->id; };
        MeshNode* node = std::find_if(this->node_list.begin(), this->node_list.end(), get_id)->get();
        list.back().nodes.push_back(node);
        ++i;

        // Begin next element
        if(i == nodes_per_elem){
            // Checking for colinear points. Implemented when there some
            // problems with mesh generation, but were due to a bad imple-
            // mentation of mine. Probably not needed anymore.
            //
            // auto& nodes = list.back().nodes;
            // double Delta = 0;
            // if(node_per_elem % 3 == 0){
            //     gp_Pnt p[3] = {nodes[0]->point, nodes[1]->point, nodes[2]->point};
            //     gp_Mat deltaM(1, p[0].X(), p[0].Y(), 1, p[1].X(), p[1].Y(), 1, p[2].X(), p[2].Y());
            //     Delta = std::abs(deltaM.Determinant());
            // } else {
            //     // TODO: Checking for 3D
            // }
            // if(Delta < 1e-3){
            //     list.pop_back();
            // }

            list.emplace_back();
            i = 0;
        }
    }
    list.pop_back();

    return list;
}

std::unordered_map<size_t, size_t> Meshing::find_duplicates(){
    std::unordered_map<size_t, size_t> duplicate_map;
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
    // The usual method to find duplicates is to test each one against the
    // other, which is O(N^2). However, std::sort() is Nlog2N, then,
    // figuring out the duplicates is O(N), so it's a much better
    // alternative.
    std::sort(this->node_list.begin(), this->node_list.end(), point_sort);

    auto i = this->node_list.begin();
    while(i < this->node_list.end()-1){
        if((*i)->point.IsEqual((*(i+1))->point, Precision::Confusion())){
            duplicate_map.emplace((*(i+1))->id, (*i)->id);
            this->node_list.erase(i+1);
        } else {
            ++i;
        }
    }
    
    return duplicate_map;
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

    return has_condition_inside;
}
