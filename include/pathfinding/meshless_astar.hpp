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

#ifndef MESHLESS_ASTAR_HPP
#define MESHLESS_ASTAR_HPP

#include "pathfinding.hpp"

namespace pathfinding{

struct PathPoint{
    PathPoint(gp_Pnt point, double cost, PathPoint* prev):
        point(point), cost(cost), prev(prev){}
    gp_Pnt point;
    double cost;
    PathPoint* prev;
};

class PathPointCompare{
    public:
    inline bool operator() (PathPoint* p1, PathPoint* p2) const{
        return p1->cost > p2->cost;
    }
};

class MeshlessAStar : public Pathfinding{
    public:
    enum ProbType{
        TYPE_2D,
        TYPE_3D
    };

    MeshlessAStar(TopoDS_Shape topology, double step, double turn_angle, int choices, double restriction, ProbType type);
    ~MeshlessAStar() = default;

    virtual std::vector<gp_Pnt> find_path(gp_Pnt p, const TopoDS_Shape& dest, gp_Dir initial_direction) override;

    private:
    double step;
    std::vector<double> angles2D;
    std::vector<std::array<double,2>> angles3D;
    double restriction;
    ProbType type;
    TopoDS_Shape topology;

    bool is_inside(gp_Pnt p, const TopoDS_Shape& t);
    double get_distance(gp_Pnt p, const TopoDS_Shape& t);
    std::pair<bool, gp_Pnt> get_intersection_point(gp_Pnt p, gp_Dir dir, const TopoDS_Shape& t);

};

}

#endif
