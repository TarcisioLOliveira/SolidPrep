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

#ifndef VISIBILITY_GRAPH_HPP
#define VISIBILITY_GRAPH_HPP

#include "TopoDS_Solid.hxx"
#include "pathfinding.hpp"
#include "utils.hpp"
#include "geometry.hpp"
#include "Geom_Curve.hxx"
#include "TopoDS_Face.hxx"
#include <memory>
#include <queue>
#include <vector>

namespace pathfinding{

class VisibilityGraph : public Pathfinding{
    public:
    struct PathPoint;

    class PathPointCompare{
        public:
        inline bool operator() (PathPoint* p1, PathPoint* p2) const{
            return p1->cost > p2->cost;
        }
    };

    struct PathPoint{
        PathPoint(size_t point, double cost, PathPoint* prev):
            point(point), cost(cost), prev(prev){}
        size_t point;
        double cost;
        PathPoint* prev;
        std::vector<size_t> remaining;
        std::priority_queue<PathPoint*, std::vector<PathPoint*>, PathPointCompare> visible;
    };

    struct Edge{
        Handle(Geom_Curve) curve;
        double a;
        double b;
    };

    VisibilityGraph(const std::unique_ptr<Geometry>& topology, double step, double turn_angle, double restriction, utils::ProblemType type);
    ~VisibilityGraph() = default;

    virtual std::vector<gp_Pnt> find_path(const CrossSection& begin, const CrossSection& end) override;

    private:
    double step;
    double angle;
    double restriction;
    const Geometry* topology;
    utils::ProblemType type;

    gp_Pnt get_closest_point(const gp_Pnt& p, const TopoDS_Shape& t) const;
    std::vector<gp_Pnt> path_section(const CrossSection& begin, const CrossSection& end);

    bool visible_from_here(const gp_Pnt& p1, const gp_Pnt& p2, bool starter, const std::vector<Edge>& edges, const std::vector<TopoDS_Face>& faces) const;
    bool intersects_edge_2D(const gp_Pnt& p1, const gp_Pnt& p2, bool starter, const std::vector<Edge>& edges) const;
    bool intersects_edge_3D(const gp_Pnt& p1, const gp_Pnt& p2, bool starter, const std::vector<TopoDS_Face>& faces) const;
};

}

#endif
