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

#include "pathfinding/visibility_graph.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BOPTools_AlgoTools.hxx>
#include <gp.hxx>
#include <gp_Circ.hxx>
#include <TopExp.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BOPTools_AlgoTools.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS.hxx>
#include <IntCurvesFace_ShapeIntersector.hxx>
#include <IntTools_Context.hxx>
#include <TopoDS_Shape.hxx>
#include <BRep_Tool.hxx>
#include <IntTools_EdgeEdge.hxx>
#include <vector>
#include <algorithm>
#include <memory>
#include <GeomAPI_ExtremaCurveCurve.hxx>
#include <lapacke.h>

namespace pathfinding{

VisibilityGraph::VisibilityGraph(const std::unique_ptr<Geometry>& topology, double step, double turn_angle, double restriction, utils::ProblemType type):
    step(step), angle(turn_angle*M_PI/180), restriction(restriction), topology(topology.get()), type(type){}

std::vector<gp_Pnt> VisibilityGraph::find_path(const CrossSection& begin, const CrossSection& end){
    struct Vertex{
        gp_Pnt v;
        std::vector<size_t> linked_vertices;
    };

    double restr = this->restriction + begin.get_dimension()/2;

    // Get edges and vertices.
    std::vector<Edge> edges;
    std::vector<Vertex> vertices;
    for (TopExp_Explorer exp(this->topology->shape, TopAbs_EDGE); exp.More(); exp.Next()){
        TopoDS_Edge edge = TopoDS::Edge(exp.Current());
        double c1 = 0;
        double c2 = 0;
        Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, c1, c2);
        Edge e{curve, c1, c2};
        edges.push_back(e);

        gp_Pnt p1 = curve->Value(c1);
        gp_Pnt p2 = curve->Value(c2);
        if(utils::equal(p1, p2, 1e-3)){
            const int SUBDIV = 8;
            double cstep = c2-c1/SUBDIV;
            size_t pos = vertices.size();
            for(size_t i = SUBDIV; i < 2*SUBDIV-1; ++i){
                Vertex v;
                v.v = curve->Value(c1 + (i%SUBDIV)*cstep);
                v.linked_vertices.push_back(pos + (i-1)%SUBDIV);
                v.linked_vertices.push_back(pos + (i+1)%SUBDIV);
                vertices.push_back(v);
            }

        } else {
            bool found_v1 = false;
            size_t v1_pos = vertices.size();
            for(size_t i = 0; i < vertices.size(); ++i){
                if(utils::equal(p1, vertices[i].v, 1e-3)){
                    found_v1 = true;
                    v1_pos = i;
                    break;
                }
            }
            bool found_v2 = false;
            size_t v2_pos = vertices.size();
            if(!found_v1){
                ++v2_pos;
            }
            for(size_t i = 0; i < vertices.size(); ++i){
                if(utils::equal(p2, vertices[i].v, 1e-3)){
                    found_v2 = true;
                    v2_pos = i;
                    break;
                }
            }
            if(found_v1){
                vertices[v1_pos].linked_vertices.push_back(v2_pos);
            } else {
                Vertex v1;
                v1.v = p1;
                v1.linked_vertices.push_back(v2_pos);
                vertices.push_back(v1);
            }
            if(found_v2){
                vertices[v2_pos].linked_vertices.push_back(v1_pos);
            } else {
                Vertex v2;
                v2.v = p2;
                v2.linked_vertices.push_back(v1_pos);
                vertices.push_back(v2);
            }
        }
    }

    // Path node generation
    std::vector<gp_Pnt> node_list;
    for(auto& v:vertices){
        gp_Vec sum(0,0,0);
        for(auto& i:v.linked_vertices){
            sum += gp_Vec(v.v, vertices[i].v).Normalized();
        }
        if(sum.SquareMagnitude() == 0){
            continue;
        }
        gp_Vec dir(-sum.Normalized());
        if(this->topology->is_inside(v.v.Translated(10*dir))){
            node_list.push_back(v.v.Translated(restr*dir));
        }
    }

    std::vector<TopoDS_Face> faces;
    if(this->type == utils::PROBLEM_TYPE_3D){
        for (TopExp_Explorer exp(this->topology->shape, TopAbs_FACE); exp.More(); exp.Next()){
            TopoDS_Face face = TopoDS::Face(exp.Current());
            faces.push_back(face);
        }
    }

    // Get node path
    std::vector<std::unique_ptr<PathPoint>> full_list;
    node_list.push_back(begin.get_centroid());
    full_list.emplace_back(new PathPoint(node_list.size()-1, end.get_distance(begin.get_centroid()), nullptr));
    for(size_t i = 0; i < node_list.size()-1; ++i){
        full_list[0]->remaining.push_back(i);
    }
    PathPoint* current = full_list.front().get();
    auto queue = &current->visible;
    bool end_visible = this->visible_from_here(node_list[current->point], end.get_centroid(), true, edges, faces);
    while(!end_visible){
        std::vector<size_t> non_visible;
        std::vector<PathPoint*> visible;
        for(const auto& p:current->remaining){
            if(this->visible_from_here(node_list[current->point], node_list[p], false, edges, faces)){
                full_list.emplace_back(new PathPoint(p, node_list[p].Distance(node_list[current->point])+end.get_distance(node_list[p]), current));
                visible.push_back(full_list.back().get());
                current->visible.push(visible.back());
            } else {
                non_visible.push_back(p);
            }
        }
        if(visible.size() > 0){
            for(auto& p:visible){
                p->remaining = non_visible;
            }
            queue = &current->visible;
            current = queue->top();
            queue->pop();
        } else {
            if(!queue->empty()){
                current = queue->top();
                queue->pop();
            } else if(current->prev != nullptr){
                queue = &current->prev->visible;
                current = queue->top();
                queue->pop();
            } else {
                break;
            }
        }
        end_visible = this->visible_from_here(node_list[current->point], end.get_centroid(), false, edges, faces);
    }
    std::vector<size_t> node_path;
    if(end_visible){
        do{
            node_path.push_back(current->point);
            current = current->prev;
        } while(current != nullptr);
    } else {
        logger::log_assert(false, logger::ERROR, "pathfinder was unable to find path.");
        return std::vector<gp_Pnt>();
    }

    // Get final path
    if(node_path.size() == 1){
        std::vector<gp_Pnt> final_list(this->path_section(begin, end));
        std::reverse(final_list.begin(), final_list.end());
        return final_list;
    } else if(node_path.size() > 1){
        gp_Pnt cur_point = node_list[*(node_path.end()-2)];
        std::vector<gp_Pnt> final_list(this->path_section(begin, CrossSection(cur_point)));
        std::reverse(final_list.begin(), final_list.end());
        for(auto i = node_path.rbegin()+2; i < node_path.rend(); ++i){
            cur_point = node_list[*i];
            gp_Dir n(gp_Vec(*(final_list.begin()+2), *(final_list.begin()+1)));
            gp_Pnt prev_point = final_list.front();
            CrossSection cs(prev_point);
            cs.set_normal(n);
            std::vector<gp_Pnt> cur_list(this->path_section(cs, CrossSection(cur_point)));
            final_list.insert(final_list.begin(), cur_list.rbegin(), cur_list.rend()-1);
        }
        gp_Dir n(gp_Vec(*(final_list.begin()+2), *(final_list.begin()+1)));
        gp_Pnt prev_point = final_list.front();
        CrossSection cs(prev_point);
        cs.set_normal(n);
        std::vector<gp_Pnt> cur_list(this->path_section(cs, end));
        final_list.insert(final_list.begin(), cur_list.rbegin(), cur_list.rend()-1);

        return final_list;
    }

    return std::vector<gp_Pnt>();
}
std::vector<gp_Pnt> VisibilityGraph::path_section(const CrossSection& begin, const CrossSection& end){
    std::vector<gp_Pnt> list;
    list.push_back(begin.get_centroid());
    gp_Pnt b1 = begin.get_centroid().Translated(this->step*begin.get_normal());
    gp_Pnt b2 = begin.get_centroid().Translated(-this->step*begin.get_normal());

    if(begin.get_dimension() == 0){
        b2 = b1;
    }

    gp_Pnt closest1 = this->get_closest_point(b1, end.get_shape());
    gp_Pnt closest2 = this->get_closest_point(b2, end.get_shape());
    double dist1 = b1.Distance(closest1);
    double dist2 = b2.Distance(closest2);

    double prev_dist = 0;
    double curr_dist = 0;
    if(dist1 < dist2){
        list.push_back(b1);
        curr_dist = dist1;
    } else {
        list.push_back(b2);
        curr_dist = dist2;
    }
    do{
        prev_dist = curr_dist;
        gp_Pnt current = list.back();
        gp_Pnt prev = *(list.end()-2);
        gp_Pnt closest = this->get_closest_point(current, end.get_shape());
        gp_Dir current_dir(gp_Vec(prev, current));
        gp_Dir path_dir(gp_Vec(current, closest));
        curr_dist = current.Distance(closest);
        if(std::abs(current_dir.AngleWithRef(path_dir, gp_Dir(0,0,1))) <= this->angle){
            if(curr_dist <= this->step + 1e-3){
                list.push_back(closest);
                break;
            } else {
                list.push_back(current.Translated(this->step*path_dir));
            }
        } else {
            if(this->type == utils::PROBLEM_TYPE_2D){
                double needed_ang = current_dir.AngleWithRef(path_dir, gp_Dir(0,0,1));
                double ang = this->angle*(needed_ang/std::abs(needed_ang));
                gp_Ax1 axis(current, gp_Dir(0,0,1));
                list.push_back(current.Translated(this->step*current_dir.Rotated(axis, ang)));
            } else if(this->type == utils::PROBLEM_TYPE_3D){
                // TODO
            }
        }
        if(!this->topology->is_inside(list.back())){
            logger::log_assert(curr_dist <= prev_dist, logger::ERROR, "pathfinding algorithm did not converge.");
        }
    } while(true);//curr_dist <= prev_dist + Precision::Confusion());

    logger::log_assert(curr_dist <= prev_dist, logger::ERROR, "pathfinding algorithm did not converge.");

    return list;
}

gp_Pnt VisibilityGraph::get_closest_point(const gp_Pnt& p, const TopoDS_Shape& t) const{
    if(t.ShapeType() == TopAbs_VERTEX){
        return BRep_Tool::Pnt(TopoDS::Vertex(t));
    } else {
        TopoDS_Vertex v = BRepBuilderAPI_MakeVertex(p);
        BRepExtrema_DistShapeShape d(v, t, 0.0001, Extrema_ExtFlag_MINMAX, Extrema_ExtAlgo_Grad);
        gp_Pnt p(d.PointOnShape2(1));

        // Workaround
        if(this->type == utils::PROBLEM_TYPE_2D){
            p.SetZ(0);
        }
        return p;
    }
}

bool VisibilityGraph::visible_from_here(const gp_Pnt& p1, const gp_Pnt& p2, bool starter, const std::vector<Edge>& edges, const std::vector<TopoDS_Face>& faces) const{
    if(this->type == utils::PROBLEM_TYPE_2D){
        return !this->intersects_edge_2D(p1, p2, starter, edges);
    } else if(this->type == utils::PROBLEM_TYPE_3D){
        return !this->intersects_edge_3D(p1, p2, starter, faces);
    }
    return false;
}

bool VisibilityGraph::intersects_edge_2D(const gp_Pnt& p1, const gp_Pnt& p2, bool starter, const std::vector<Edge>& edges) const{
    // Check if Vertex on Edge counts as intersection.

    size_t allowance = 1; // p2 is always in an edge
    if(starter) ++allowance;
    TopoDS_Edge line = BRepBuilderAPI_MakeEdge(p1, p2);
    size_t int_count = 0;
    double a = 0;
    double b = 0;
    auto curve = BRep_Tool::Curve(line, a, b);
    for(auto& e:edges){
        GeomAPI_ExtremaCurveCurve cc(curve, e.curve, a-0.01, b+0.01, e.a-0.01, e.b+0.01);
        if(cc.NbExtrema() > 0 && cc.LowerDistance() < Precision::Confusion()){
            ++int_count;
            if(int_count > allowance){
                return true;
            }
        }
    }

    return false;
}

bool VisibilityGraph::intersects_edge_3D(const gp_Pnt& p1, const gp_Pnt& p2, bool starter, const std::vector<TopoDS_Face>& faces) const{
    // Check if Vertex on Face counts as intersection.

    (void)starter;
    for(auto& f:faces){
        gp_Lin line(p1, gp_Vec(p1, p2));
        IntCurvesFace_ShapeIntersector intersector;
        intersector.Load(f, 0.01);
        intersector.Perform(line, 0, p1.Distance(p2));
        if(intersector.NbPnt() > 0){
            return true;
        }
    }
    return false;
}


}

/*
    // Attempt at creating the final path using splines. Hopefully it will work
    // one day. The tricky part is to get the function to pick the best point
    // in the the final cross section in order to obtain the shortest path.
    // Other than that, it works really well.

    std::vector<double> dists(node_path.size()+1);
    double total_dist = end.get_distance(node_list[node_path.front()]);
    dists[0] = total_dist;
    for(size_t i = 1; i < node_path.size(); ++i){
        double d = node_list[node_path[i]].Distance(node_list[node_path[i-1]]);
        dists[i] = d;
        total_dist += d;
    }
    std::vector<double> t(dists.size());
    t[0] = 1;
    for(size_t i = 0; i < t.size(); ++i){
        t[i+1] = t[i] - dists[i]/total_dist;
        std::cout << t[i] << std::endl;
    }
    t.back() = 0;

    if(this->type == utils::PROBLEM_TYPE_2D){
        size_t W = (node_path.size()+1+2)*2; // n+1 points (includes unknown end point) + 2 derivatives, 2D
        size_t H = W-1+node_path.size(); // 1 equation for every coefficient + a restriction to obtain end point
        std::vector<double> A(W*H);
        std::vector<double> B(W);
        size_t yoffsetA = W*(W/2-1)+W/2;
        size_t yoffsetB = W/2-1;
        // Line equation from end point (t = 1)
        gp_Vec end_normal = end.get_normal();
        B[0] = end_normal.Dot(gp_Vec(gp_Pnt(0,0,0), end.get_centroid()));
        for(size_t i = 0; i < W/2; ++i){
            A[i] = end_normal.X();
            A[i+yoffsetB+1] = end_normal.Y();
        }
        // Derivative from start point (t = 0)
        B[1] = begin.get_normal().X();
        B[1+yoffsetB] = begin.get_normal().Y();
        A[W+1] = 1;
        A[W+1+yoffsetA] = 1;
        // Derivative from end point (t = 1)
        B[2] = end.get_normal().X();
        B[2+yoffsetB] = end.get_normal().Y();

        for(size_t i = 2*W+1; i < 2*W+W/2; ++i){
            A[i] = i-2*W;
            A[i+yoffsetA] = i-2*W;
        }

        for(size_t i = 1; i < t.size(); ++i){
            B[3+i-1] = node_list[node_path[i-1]].X();
            B[3+i-1+yoffsetB] = node_list[node_path[i-1]].Y();
            for(size_t j = 0; j < W/2; ++j){
                A[3*W+(i-1)*W+j] = std::pow(t[i], j);
                A[3*W+(i-1)*W+j+yoffsetA] = std::pow(t[i], j);
            }
        }

        size_t offset = W*(H - node_path.size());
        for(size_t i = 0; i < node_path.size()-1; ++i){
            gp_Dir n = node_dirs[node_path[i]];
            for(size_t j = 1; j < W/2; ++j){
                A[offset+i*W+j] = j*n.X()*std::pow(t[i+1], j-1);
                A[offset+i*W+j+yoffsetB+1] = j*n.Y()*std::pow(t[i+1], j-1);
            }
        }
        offset = (W-1)*H+1;
        gp_Dir n = end.get_normal();
        for(size_t j = 1; j < W/2; ++j){
            A[offset+j] = j*n.X();
            A[offset+j+yoffsetB+1] = j*n.Y();
        }

        for(size_t i = 0; i < W*H; ++i){
            std::cout << A[i] << " ";
            if((i+1)%W == 0){
                std::cout << std::endl;
            }
        }

        std::cout << std::endl;
        for(auto& i:B){
            std::cout << i << " ";
        }

        std::cout << std::endl;
        std::cout << std::endl;
        LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', H, W, 1, A.data(), W, B.data(), 1);
        for(auto& i:B){
            std::cout << i << " ";
        }
        std::cout << std::endl;
        std::cout << std::endl;
        double X = 0;
        double Y = 0;
        for(size_t i = 0; i < W/2; ++i){
            X += B[i];
            Y += B[i+yoffsetB+1];
        }
        std::cout << X << " " << Y << std::endl;

    }
*/
