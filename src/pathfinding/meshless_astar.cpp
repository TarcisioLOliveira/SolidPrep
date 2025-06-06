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

#include "pathfinding/meshless_astar.hpp"
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
#include <BRep_Tool.hxx>
#include <vector>
#include <algorithm>
#include <memory>
#include "project_data.hpp"

namespace pathfinding{

void MeshlessAStar::PriorityQueue::push(PathPoint* p){
    if(this->empty()){
        this->c.push_back(p);
    }
    auto i = this->c.crbegin();
    while(i < this->c.crend()){
        if((*i)->cost > p->cost){
            break;
        }
        ++i;
    }
    auto base = i.base();
    if(base == this->c.end()){
        this->c.push_back(p);
    } else if(this->equal((*base)->point, p->point)){
        return;
    } else {
        this->c.insert(base, p);
    }
}

bool MeshlessAStar::PriorityQueue::equal(gp_Pnt p1, gp_Pnt p2, double eps){
    return ((p1.X() - p2.X()) < eps) && ((p2.X() - p1.X()) < eps) &&
           ((p1.Y() - p2.Y()) < eps) && ((p2.Y() - p1.Y()) < eps) &&
           ((p1.Z() - p2.Z()) < eps) && ((p2.Z() - p1.Z()) < eps);
}

MeshlessAStar::MeshlessAStar(const projspec::DataMap& data):
    step(data.get_double("step")),
    turn_angle(data.get_double("max_turn_angle")),
    angles2D(),
    angles3D(),
    restriction(data.get_double("restriction_size")),
    type(data.proj->type),
    topology(TopoDS::Solid(data.proj->geometries[0]->shape)){

    const auto choices = data.get_int("turn_options");

    const double angle = (M_PI/180.0)*this->turn_angle;
    if(type == utils::PROBLEM_TYPE_2D){
        this->topology = BRepBuilderAPI_MakeSolid(TopoDS::Shell(topology));
    } else {
        this->topology = TopoDS::Solid(topology);
    }
    if(this->type == utils::PROBLEM_TYPE_2D){
        double angle_step = 2*angle/(choices-1);
        for(int i = 0; i < choices; ++i){
            this->angles2D.push_back(-angle + i*angle_step);
        }
    } else if(this->type == utils::PROBLEM_TYPE_3D){
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

std::vector<gp_Pnt> MeshlessAStar::find_path(const CrossSection& begin, const CrossSection& end){
    TopoDS_Shape dest = end.get_shape();
    PriorityQueue point_queue;
    std::vector<std::unique_ptr<PathPoint>> point_list;
    gp_Pnt p = begin.get_centroid();
    gp_Dir initial_direction = begin.get_normal();

    point_list.emplace_back(new PathPoint(p, end.get_distance(p), nullptr));

    // Forward direction
    gp_Pnt point = p.Translated(this->step*initial_direction);
    point_list.emplace_back(new PathPoint(point, end.get_distance(p), point_list[0].get()));
    point_queue.push(point_list[point_list.size()-1].get());
    // Backward direction
    point = p.Translated(this->step*(-initial_direction));
    point_list.emplace_back(new PathPoint(point, end.get_distance(p), point_list[0].get()));
    point_queue.push(point_list[point_list.size()-1].get());

    PathPoint* current = point_queue.top();
    gp_Dir direction = gp_Vec(current->prev->point, current->point);
    point_queue.pop();

    double f_dim = begin.get_dimension()*1e3;

    bool reached_obj = false;
    if(this->type == utils::PROBLEM_TYPE_2D){
        while(!reached_obj){
            bool fully_inside_topology = this->shape_inside_2D(current->point, direction, this->restriction + f_dim/2, this->step, this->topology);
            bool center_inside = this->is_inside_2D(current->point, this->topology);
            bool close_to_start = current->point.Distance(point_list[0]->point) <= (this->restriction + f_dim/2);
            bool close_to_end = end.get_distance(current->point) <= this->step;
            if(center_inside && (fully_inside_topology || close_to_start || close_to_end)){
                gp_Ax1 axis(current->point, gp_Dir(0.0,0.0,1.0));
                for(double a:this->angles2D){
                    gp_Dir dir = direction.Rotated(axis, a);
                    std::pair<bool, gp_Pnt> intersec = this->get_intersection_point(current->point, dir, this->step, dest);
                    if(intersec.first == false){
                        gp_Pnt point = current->point.Translated(this->step*dir);
                        point_list.emplace_back(new PathPoint(point, end.get_distance(point), current));
                        point_queue.push(point_list[point_list.size()-1].get());
                    } else {
                        gp_Pnt point = intersec.second;
                        point_list.emplace_back(new PathPoint(point, current->point.Distance(point), current));
                        point_queue.push(point_list[point_list.size()-1].get());
                        reached_obj = true;
                    }
                }
            }
            if(point_queue.empty()){
                break;
            }
            current = point_queue.top();
            direction = gp_Vec(current->prev->point, current->point);
            point_queue.pop();
        }
    } else if(this->type == utils::PROBLEM_TYPE_3D){
        // TODO
    }

    logger::log_assert(reached_obj, logger::ERROR, "Failed to find path. Initial point was ({},{},{})", p.X(), p.Y(), p.Z());
    std::vector<gp_Pnt> path;
    while(current != nullptr){
        path.push_back(current->point);
        current = current->prev;
    }

    return path;
}

bool MeshlessAStar::shape_inside_2D(gp_Pnt center, gp_Dir dir, double restr, double step, const TopoDS_Shape& t){
    std::vector<gp_Pnt> p({gp_Pnt(center.X() + step/2, center.Y() + restr, 0), 
                           gp_Pnt(center.X() - step/2, center.Y() + restr, 0), 
                           gp_Pnt(center.X() - step/2, center.Y() - restr, 0), 
                           gp_Pnt(center.X() + step/2, center.Y() - restr, 0),
                           gp_Pnt(center.X() - step/2, center.Y(), 0), 
                           gp_Pnt(center.X() + step/2, center.Y(), 0), 
                           gp_Pnt(center.X(), center.Y() + restr, 0), 
                           gp_Pnt(center.X(), center.Y() - restr, 0)}); 
    gp_Ax1 axis(center, gp_Dir(0, 0, 1));
    double ang = dir.AngleWithRef(gp_Dir(1,0,0), gp_Dir(0,0,1));
    bool inside = true;
    for(size_t i = 0; i < p.size(); ++i){
        p[i].Rotate(axis, ang);
        inside = inside && this->is_inside_2D(p[i], t);
        if(!inside){
            break;
        }
    }

    return inside;
}

bool MeshlessAStar::is_inside_2D(gp_Pnt p, const TopoDS_Shape& t){
    BRepClass3d_SolidClassifier insider(t, p, 0.01);
    return insider.State() == TopAbs_ON;
}

bool MeshlessAStar::is_inside_3D(gp_Pnt p, const TopoDS_Shape& t){
    BRepClass3d_SolidClassifier insider(t, p, 0.01);
    return insider.State() == TopAbs_IN;
}

std::pair<bool, gp_Pnt> MeshlessAStar::get_intersection_point(gp_Pnt p, gp_Dir dir, double step, const TopoDS_Shape& t){
    if(t.ShapeType() == TopAbs_VERTEX){
        gp_Pnt v = BRep_Tool::Pnt(TopoDS::Vertex(t));
        if(v.Distance(p) < step){
            gp_Vec vec(p, v);
            double ang = dir.Angle(vec);
            if(ang <= this->turn_angle){
                return std::pair<bool, gp_Pnt>(true, v);
            }
        }
        return std::pair<bool, gp_Pnt>(false, v);
    } else {
        gp_Pnt maxp = p.Translated(step*dir);
        if(this->is_inside_2D(maxp, t)){
            return std::pair<bool, gp_Pnt>(true, maxp);
        }
        gp_Lin line(p, dir);
        gp_Pnt itsc(0,0,0);
        IntCurvesFace_ShapeIntersector intersector;
        intersector.Load(t, 0.01);
        intersector.PerformNearest(line, 0, step);
        if(intersector.NbPnt() > 0){
            itsc = intersector.Pnt(1);
            return std::pair<bool, gp_Pnt>(true, itsc);
        } else {
            return std::pair<bool, gp_Pnt>(false, itsc);
        }
    }
}

using namespace projspec;
const bool MeshlessAStar::reg = Factory<Pathfinding>::add(
    [](const DataMap& data){
        return std::make_unique<MeshlessAStar>(data);
    },
    ObjectRequirements{
        "meshless_astar",
        {
            DataEntry{.name = "step", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "max_turn_angle", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "turn_options", .type = TYPE_INT, .required = true},
            DataEntry{.name = "restriction_size", .type = TYPE_DOUBLE, .required = false},
        }
    }
);

}
