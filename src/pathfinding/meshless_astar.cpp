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
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BOPTools_AlgoTools.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS.hxx>
#include <IntCurvesFace_ShapeIntersector.hxx>
#include <IntTools_Context.hxx>
#include <vector>
#include <queue>
#include <algorithm>
#include <memory>

namespace pathfinding{

MeshlessAStar::MeshlessAStar(TopoDS_Shape topology, double step, double angle, int choices, double restriction, ProbType type){
    this->step = step;
    this->restriction = restriction;
    this->type = type;
    angle = (M_PI/180)*angle;
    if(type == TYPE_2D){
        this->topology = BRepBuilderAPI_MakeSolid(TopoDS::Shell(topology));
    } else {
        this->topology = TopoDS::Solid(topology);
    }
    if(this->type == TYPE_2D){
        double angle_step = 2*angle/(choices-1);
        for(int i = 0; i < choices; ++i){
            this->angles2D.push_back(-angle + i*angle_step);
        }
    } else if(this->type == TYPE_3D){
        double angle_step = 2*angle/(choices-1);
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

std::vector<gp_Pnt> MeshlessAStar::find_path(const Force& f, const TopoDS_Shape& dest){
    std::priority_queue<PathPoint*, std::vector<PathPoint*>, PathPointCompare> point_queue;
    std::vector<std::unique_ptr<PathPoint>> point_list;
    gp_Pnt p = f.get_centroid();
    gp_Dir initial_direction = f.get_normal();

    point_list.emplace_back(new PathPoint(p, this->get_distance(p, dest), nullptr));

    // Forward direction
    gp_Pnt point = p.Translated(this->step*initial_direction);
    point_list.emplace_back(new PathPoint(point, this->get_distance(point, dest), point_list[0].get()));
    point_queue.push(point_list[point_list.size()-1].get());
    // Backward direction
    point = p.Translated(this->step*(-initial_direction));
    point_list.emplace_back(new PathPoint(point, this->get_distance(point, dest), point_list[0].get()));
    point_queue.push(point_list[point_list.size()-1].get());

    PathPoint* current = point_queue.top();
    gp_Dir direction = gp_Vec(current->prev->point, current->point);
    point_queue.pop();

    if(this->type == TYPE_2D){
        bool reached_obj = false;
        while(!reached_obj && !point_queue.empty()){
            std::cout << current->point.X() << " " << current->point.Y() << std::endl;
            bool fully_inside_topology = this->shape_inside_2D(current->point, direction, this->restriction + f.get_dimension()/2, this->step, this->topology);
            bool center_inside = this->is_inside(current->point, this->topology);
            bool close_to_start = current->point.Distance(point_list[0]->point) <= this->restriction + f.get_dimension()/2;
            bool close_to_end = this->get_distance(current->point, dest) <= this->restriction + f.get_dimension()/2;
            if(center_inside && (fully_inside_topology || close_to_start || close_to_end)){
                std::cout << "inside" << std::endl;
                gp_Ax1 axis(current->point, gp_Dir(0.0,0.0,1.0));
                for(double a:this->angles2D){
                    gp_Dir dir = direction.Rotated(axis, a);
                    std::pair<bool, gp_Pnt> intersec = this->get_intersection_point(current->point, dir, this->step, dest);
                    if(intersec.first == false){
                        gp_Pnt point = current->point.Translated(this->step*dir);
                        point_list.emplace_back(new PathPoint(point, this->get_distance(point, dest), current));
                        point_queue.push(point_list[point_list.size()-1].get());
                    } else {
                        gp_Pnt point = intersec.second;
                        point_list.emplace_back(new PathPoint(point, this->get_distance(point, dest), current));
                        point_queue.push(point_list[point_list.size()-1].get());
                        reached_obj = true;
                        break;
                    }
                }
            }
            else {
                std::cout << center_inside << " " << fully_inside_topology << " " << close_to_start << " " << close_to_end << std::endl;
            }
            current = point_queue.top();
            direction = gp_Vec(current->prev->point, current->point);
            point_queue.pop();
        }
    } else if(this->type == TYPE_3D){
        // TODO
    }

    logger::log_assert(!point_queue.empty(), logger::ERROR, "Failed to find path. Initial point was ({},{},{})", p.X(), p.Y(), p.Z());
    std::vector<gp_Pnt> path;
    while(current != nullptr){
        path.push_back(current->point);
        current = current->prev;
    }

    return path;
}

bool MeshlessAStar::shape_inside_2D(gp_Pnt center, gp_Dir dir, double restr, double step, const TopoDS_Shape& t){
    gp_Pnt p[4] = {gp_Pnt(center.X() + step/2, center.Y() + restr, 0), 
                   gp_Pnt(center.X() - step/2, center.Y() + restr, 0), 
                   gp_Pnt(center.X() - step/2, center.Y() - restr, 0), 
                   gp_Pnt(center.X() + step/2, center.Y() - restr, 0)};
    gp_Ax1 axis(center, gp_Dir(0, 0, 1));
    double ang = dir.AngleWithRef(gp_Dir(1,0,0), gp_Dir(0,0,1));
    bool inside = true;
    for(int i = 0; i < 4; ++i){
        p[i].Rotate(axis, ang);
        inside = inside && this->is_inside(p[i], t);
        if(!inside){
            break;
        }
    }
    return inside;
}

bool MeshlessAStar::is_inside(gp_Pnt p, const TopoDS_Shape& t){
    BRepClass3d_SolidClassifier insider(t, p, 0.001);
    return insider.State() == TopAbs_ON || insider.State() == TopAbs_IN;
}

double MeshlessAStar::get_distance(gp_Pnt p, const TopoDS_Shape& t){
    TopoDS_Vertex v = BRepBuilderAPI_MakeVertex(p);

    BRepExtrema_DistShapeShape d(v, t, Extrema_ExtFlag_MIN, Extrema_ExtAlgo_Grad);
    return d.Value();

}

std::pair<bool, gp_Pnt> MeshlessAStar::get_intersection_point(gp_Pnt p, gp_Dir dir, double step, const TopoDS_Shape& t){
    gp_Lin line(p, dir);
    gp_Pnt itsc(0,0,0);
    IntCurvesFace_ShapeIntersector intersector;
    intersector.Load(t, 0.001);
    intersector.PerformNearest(line, 0, step);
    try{
        itsc = intersector.Pnt(1);
        return std::pair<bool, gp_Pnt>(true, itsc);
    } catch(const Standard_OutOfRange& e){
        return std::pair<bool, gp_Pnt>(false, itsc);
    }
    return std::pair<bool, gp_Pnt>(false, itsc);
}

}
