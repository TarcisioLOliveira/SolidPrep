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

#include "pathfinding/meshless_astar.hpp"
#include "logger.hpp"
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
#include <BOPTools_AlgoTools.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS.hxx>
#include <IntCurvesFace_ShapeIntersector.hxx>
#include <vector>
#include <queue>
#include <algorithm>

namespace pathfinding{

MeshlessAStar::MeshlessAStar(TopoDS_Shape topology, double step, double angle, int choices, double restriction, ProbType type){
    this->step = step;
    this->restriction = restriction;
    this->type = type;
    if(this->type == TYPE_2D){
        double angle_step = (M_PI/360)*2*angle/(choices-1);
        for(int i = 0; i < choices; ++i){
            this->angles2D.push_back(-angle + i*angle_step);
        }
    } else if(this->type == TYPE_3D){
        double angle_step = (M_PI/360)*2*angle/(choices-1);
        double spin_step = 2*M_PI/(choices-1);
        bool odd = (choices % 2) == 1;
        if(odd){
            std::array<double, 2> arr = {0, 0};
            this->angles3D.push_back(arr);
            for(int j = 0; j < choices - 1; ++j){
                for(int i = 0; i < choices; ++i){
                    arr[0] = -angle + i*angle_step;
                    if(arr[0] != 0){
                        arr[1] = j*spin_step;
                        this->angles3D.push_back(arr);
                    }
                }
            }
        } else {
            std::array<double, 2> arr = {0, 0};
            for(int j = 0; j < choices+1; ++j){
                for(int i = 0; i < choices; ++i){
                    arr[0] = -angle + i*angle_step;
                    arr[1] = j*spin_step;
                    this->angles3D.push_back(arr);
                }
            }
        }
    }
}

std::vector<gp_Pnt> MeshlessAStar::find_path(gp_Pnt p, const TopoDS_Shape& dest, gp_Dir initial_direction){
    std::priority_queue<PathPoint*, std::vector<PathPoint*>, PathPointCompare> point_queue;
    std::vector<PathPoint> point_list;
    point_list.emplace_back(p, this->get_distance(p, dest), nullptr);
    point_queue.push(&point_list[0]);

    // Forward direction
    gp_Pnt point = p.Translated(step*initial_direction);
    point_list.emplace_back(point, this->get_distance(point, dest), &point_list[0]);
    point_queue.push(&point_list[point_list.size()-1]);
    // Backward direction
    point = p.Translated(step*(-initial_direction));
    point_list.emplace_back(point, this->get_distance(point, dest), &point_list[0]);
    point_queue.push(&point_list[point_list.size()-1]);

    PathPoint* current = point_queue.top();
    gp_Dir direction = initial_direction; // GET DIRECTION FROM SUBTRACTING POINTS
    point_queue.pop();

    if(this->type == TYPE_2D){
        gp_Ax1 axis(p, gp_Dir(0.0,0.0,1.0));
        TopTools_IndexedMapOfShape edges;
        TopExp::MapShapes(this->topology, TopAbs_EDGE, edges);
        Handle(IntTools_Context) context;
        gp_Circ circle(gp::XOY(), this->restriction/2);
        TopoDS_Edge cedge = BRepBuilderAPI_MakeEdge(circle);
        TopoDS_Wire cwire = BRepBuilderAPI_MakeWire(cedge);
        TopoDS_Face cface = BRepBuilderAPI_MakeFace(cwire);
        TopoDS_Solid solid_topo = TopoDS::Solid(this->topology);

        gp_Trsf translation;
        bool reached_obj = false;
        while(!reached_obj && !point_queue.empty()){
            if(this->is_inside(current->point, this->topology)){
                translation.SetTranslation(gp::Origin(), current->point);
                BRepBuilderAPI_Transform transf(cface, translation, true);
                TopoDS_Face cface_t = TopoDS::Face(transf.Shape());

                TopAbs_State state = BOPTools_AlgoTools::ComputeState(cface_t, solid_topo, 0.001, edges, context);
                if(state == TopAbs_IN || current->point.Distance(point_list[0].point) <= this->restriction){
                    for(double a:this->angles2D){
                        gp_Dir dir = direction.Rotated(axis, a);
                        std::pair<bool, gp_Pnt> intersec = this->get_intersection_point(current->point, dir, dest);
                        if(intersec.first == false){
                            gp_Pnt point = p.Translated(step*dir);
                            point_list.emplace_back(point, this->get_distance(point, dest), &point_list[0]);
                            point_queue.push(&point_list[point_list.size()-1]);
                        } else {
                            gp_Pnt point = intersec.second;
                            point_list.emplace_back(point, this->get_distance(point, dest), &point_list[0]);
                            point_queue.push(&point_list[point_list.size()-1]);
                            reached_obj = true;
                            break;
                        }
                    }
                }
            }
            current = point_queue.top();
            direction = gp_Vec(current->prev->point, current->point);
            point_queue.pop();
        }
    } else if(this->type == TYPE_3D){
        while(this->get_distance(current->point, dest) > this->step && !point_queue.empty()){
            // TODO
        }
    }

    logger::log_assert(!point_queue.empty(), logger::ERROR, "Failed to find path. Initial point was ({},{},{})", p.X(), p.Y(), p.Z());
    std::vector<gp_Pnt> path;
    while(current->prev != nullptr){
        path.push_back(current->point);
        current = current->prev;
    }

    return path;
}

bool MeshlessAStar::is_inside(gp_Pnt p, const TopoDS_Shape& t){
    BRepClass3d_SolidClassifier insider(t);
    insider.Perform(p, 0);
    return insider.State() == TopAbs_ON || insider.State() == TopAbs_IN;
}

double MeshlessAStar::get_distance(gp_Pnt p, const TopoDS_Shape& t){
    TopoDS_Vertex v = BRepBuilderAPI_MakeVertex(p);

    BRepExtrema_DistShapeShape d(v, t);
    d.Perform();
    return d.Value();

}

std::pair<bool, gp_Pnt> MeshlessAStar::get_intersection_point(gp_Pnt p, gp_Dir dir, const TopoDS_Shape& t){
    gp_Lin line(p, dir);
    gp_Pnt itsc(0,0,0);
    IntCurvesFace_ShapeIntersector intersector;
    intersector.Load(t, 0.001);
    intersector.PerformNearest(line, 0.0, this->step);
    if(intersector.State(0) == TopAbs_IN && intersector.State(0) == TopAbs_ON){
        itsc = intersector.Pnt(0);
        return std::pair<bool, gp_Pnt>(true, itsc);
    } else {
        return std::pair<bool, gp_Pnt>(false, itsc);
    }
}

}
