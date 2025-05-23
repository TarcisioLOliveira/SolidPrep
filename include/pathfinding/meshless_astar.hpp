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

#ifndef MESHLESS_ASTAR_HPP
#define MESHLESS_ASTAR_HPP

#include <array>
#include "TopoDS_Solid.hxx"
#include "pathfinding.hpp"
#include "project_specification/data_map.hpp"
#include "utils.hpp"

namespace pathfinding{

class MeshlessAStar : public Pathfinding{
    public:
    struct PathPoint{
        PathPoint(gp_Pnt point, double cost, PathPoint* prev):
            point(point), cost(cost), prev(prev){}
        gp_Pnt point;
        double cost;
        PathPoint* prev;
    };
    /**
     * Custom priority queue made to avoid inserting duplicates into the container.
     * While it may look like it slows down the process, the number of duplicates
     * actually tends to become considerably big the harder it is to find a path.
     */

    class PriorityQueue{
        public:

        void push(PathPoint* p);
        void pop(){
            this->c.pop_back();
        }
        inline PathPoint* top(){
            return *this->c.rbegin();
        }
        inline bool empty(){
            return this->c.empty();
        }

        private:

        bool equal(gp_Pnt p1, gp_Pnt p2, double eps = 0.01);

        std::vector<PathPoint*> c;
    };

    MeshlessAStar(const projspec::DataMap& data);
    ~MeshlessAStar() = default;

    virtual std::vector<gp_Pnt> find_path(const CrossSection& begin, const CrossSection& end) override;

    private:
    static const bool reg;
    double step;
    double turn_angle;
    std::vector<double> angles2D;
    std::vector<std::array<double,2>> angles3D;
    double restriction;
    utils::ProblemType type;
    TopoDS_Solid topology;

    bool shape_inside_2D(gp_Pnt center, gp_Dir dir, double restr, double step, const TopoDS_Shape& t);
    bool is_inside_2D(gp_Pnt p, const TopoDS_Shape& t);
    bool is_inside_3D(gp_Pnt p, const TopoDS_Shape& t);
    std::pair<bool, gp_Pnt> get_intersection_point(gp_Pnt p, gp_Dir dir, double step, const TopoDS_Shape& t);
};

}

#endif
