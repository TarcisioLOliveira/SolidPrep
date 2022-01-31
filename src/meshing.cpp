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
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>


void Meshing::prepare_for_FEM(const std::vector<ElementShape>& base_mesh,
                              MeshElementFactory::MeshElementType element_type,
                              ProjectData* data, bool force_only){
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

    this->type = element_type;
    this->element_list.clear();
    this->element_list.reserve(base_mesh.size());

    auto comp = [](const std::unique_ptr<MeshNode>& a, const std::unique_ptr<MeshNode>& b){ return a->id < b->id;};
    std::sort(this->node_list.begin(), this->node_list.end(), comp);

    size_t dof = MeshElementFactory::get_dof_per_node(element_type);
    size_t k_size = MeshElementFactory::get_k_dimension(element_type);

    size_t current = 0;
    for(size_t i = 0; i < this->node_list.size(); ++i){
        auto& n = this->node_list[i];
        if(n->u_pos != nullptr){
            delete[] n->u_pos;
        }
        n->u_pos = new long[dof]();
        bool supported = false;
        for(auto& s : data->supports){
            if(s.S.is_inside(n->point)){//s.S.get_distance(n->point) <= this->size/4){//
                size_t offset = 0;
                std::vector<long> sup_pos = this->get_support_dof(offset, 0, s, element_type);
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

    this->load_vector.resize(current);

    for(auto& e : base_mesh){
        this->element_list.emplace_back(MeshElementFactory::make_element(element_type, e, data));
    }


    if(data->type == utils::PROBLEM_TYPE_2D){
        if(force_only){
            for(auto& f : data->forces){
                std::vector<MeshNode*> node_list;
                for(auto& n : this->node_list){
                    if(f.S.is_inside(n->point)){//f.S.get_distance(n->point) < this->size/3){//
                        node_list.push_back(n.get());
                    }
                }
                for(auto& n : node_list){
                    for(size_t i = 0; i < dof; ++i){
                        std::vector<double> f_vec = this->get_force_dof(f, element_type);
                        if(n->u_pos[i] >= 0){
                            this->load_vector[n->u_pos[i]] += f_vec[i]/node_list.size();
                        }
                    }
                }
            }
        } else {
            // Distributed loading
            for(auto& f : data->forces){
                double norm = f.vec.Magnitude()/(data->thickness*f.S.get_dimension());
                gp_Dir dir(f.vec);

                if(this->is_strictly_inside2D(f.S.get_centroid(), this->shape)){
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

                for(auto& e : this->element_list){
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
                        fe = e->get_f(dir, norm, {list[0]->point, list[1]->point});
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
                                fe = e->get_f(dir, norm, {n1, p1});
                                // logger::quick_log(n1.X(), n1.Y(), p1.X(), p1.Y(), n2.X(), n2.Y());
                                // logger::quick_log(fe);
                                break;
                            } else if(is_between_points(n1, n2, p2)){
                                fe = e->get_f(dir, norm, {n1, p2});
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
        }
    } else if(data->type == utils::PROBLEM_TYPE_3D){
        // TODO
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


void Meshing::clear_results(){
    for(auto& n:this->node_list){
        for(size_t i = 0; i < n->get_result_size(); ++i){
            n->results[i] = 0;
        }
    }
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

