/*
 *   Copyright (C) 2023 Tarcísio Ladeia de Oliveira.
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

#ifndef UTILS_BOUNDARY_NULLIFIER_HPP
#define UTILS_BOUNDARY_NULLIFIER_HPP

#include <cstddef>
#include <vector>
#include <gp_Pnt.hxx>
#include "LIN3B.hpp"
#include "element.hpp"
#include "logger.hpp"

namespace utils{

class LineBoundary{
    public:
    std::array<const Node*, 2> edges;
    bool inner;
    BoundaryMeshElement* parent;

    bool is_connected_to(const LineBoundary& l) const{
        return (this->edges[0]->point.IsEqual(l.edges[0]->point, Precision::Confusion()) || this->edges[1]->point.IsEqual(l.edges[1]->point, Precision::Confusion())) ||
               (this->edges[1]->point.IsEqual(l.edges[0]->point, Precision::Confusion()) || this->edges[0]->point.IsEqual(l.edges[1]->point, Precision::Confusion()));
    }

    bool operator==(const LineBoundary& l) const{
        return (this->edges[0]->point.IsEqual(l.edges[0]->point, Precision::Confusion()) && this->edges[1]->point.IsEqual(l.edges[1]->point, Precision::Confusion())) ||
               (this->edges[1]->point.IsEqual(l.edges[0]->point, Precision::Confusion()) && this->edges[0]->point.IsEqual(l.edges[1]->point, Precision::Confusion()));
    }
};

// Coordinates numbers for "height" (Y) and "width" (Z)
//
// If h == y and w == z: Y=1, Z=2
// If h == z and w == y: Y=2, Z=1
template<size_t Y, size_t Z>
class BoundaryNullifier{
    public:
    BoundaryNullifier(std::vector<LineBoundary> bounds, const Material* mat, size_t mesh_size, gp_Pnt center);

    inline double height(const gp_Pnt& p) const{
        return p.Coord(1+Y) - center.Coord(1+Y);
    }
    inline double width(const gp_Pnt& p) const{
        return p.Coord(1+Z) - center.Coord(1+Z);
    }

    private:
    const gp_Pnt center;
    class F{
        public:
        F() = default;
        F(const std::vector<gp_Pnt>& points, const Material* mat, size_t mesh_size);

        double operator()(double y) const;

        private:
        std::vector<double> nodal_values;
        std::vector<LIN3B> mesh;
    } Fp, Fm;

    inline bool less(const double a, const double b) const{
        return (a - b) < -Precision::Confusion();
    }
    inline bool leq(const double a, const double b) const{
        return (a - b) <= Precision::Confusion();
    }
    inline bool greater(const double a, const double b) const{
        return (a - b) > Precision::Confusion();
    }
    inline bool geq(const double a, const double b) const{
        return (a - b) >= -Precision::Confusion();
    }
    inline bool equal(const double a, const double b) const{
        return std::abs(a - b) < Precision::Confusion();
    }

    // Assumes there is a single, external boundary!
    void order_boundary(std::vector<LineBoundary>& bounds) const;
};

template<size_t Y, size_t Z>
BoundaryNullifier<Y, Z>::BoundaryNullifier(std::vector<LineBoundary> bounds, const Material* mat, size_t mesh_size, gp_Pnt center):
    center(center){

    const auto Y_comp_pos = [&](const gp_Pnt& p1, const gp_Pnt p2)->bool{
        if(less(p1.Y(), p2.Y())){
            return true;
        } else if(equal(p1.Y(), p2.Y())){
            if(p1.Y() <= 0){
                return less(p1.Z(), p2.Z());
            } else {
                return greater(p1.Z(), p2.Z());
            }
        }
        return false;
    };
    const auto Y_comp_neg = [&](const gp_Pnt& p1, const gp_Pnt p2)->bool{
        if(less(p1.Y(), p2.Y())){
            return true;
        } else if(equal(p1.Y(), p2.Y())){
            if(p1.Y() <= 0){
                return greater(p1.Z(), p2.Z());
            } else {
                return less(p1.Z(), p2.Z());
            }
        }
        return false;
    };

    this->order_boundary(bounds);

    std::vector<gp_Pnt> points_pos, points_neg;
    points_pos.reserve(bounds.size());
    points_neg.reserve(bounds.size());

    for(const auto& l:bounds){
        gp_Pnt p1(0, this->height(l.edges[0]->point), this->width(l.edges[0]->point));
        gp_Pnt p2(0, this->height(l.edges[1]->point), this->width(l.edges[1]->point));
        if(p1.Z() >= 0 && p2.Z() >= 0){
            points_pos.push_back(p1);
            if(p1.Z() == 0){
                points_neg.push_back(p1);
            }
        } else if(p1.Z() < 0 && p2.Z() < 0){
            points_neg.push_back(p1);
        } else {
            if(p1.Z() == 0){
                points_pos.push_back(p1);
                points_neg.push_back(p1);
            } else {
                if(p1.Z() > 0){
                    points_pos.push_back(p1);
                } else if(p1.Z() < 0){
                    points_neg.push_back(p1);
                }
                const double a = -p1.Z()/(p2.Z() - p1.Z());
                const double y = a*(p2.Y() - p1.Y()) + p1.Y();
                gp_Pnt lim(0, y, 0);
                points_pos.push_back(lim);
                points_neg.push_back(lim);
            }
        }
    }

    std::sort(points_pos.begin(), points_pos.end(), Y_comp_pos);
    std::sort(points_neg.begin(), points_neg.end(), Y_comp_neg);
    for(const auto& p:points_pos){
        logger::quick_log(p.X(), p.Y(), p.Z());
    }
    logger::quick_log("");
    for(const auto& p:points_neg){
        logger::quick_log(p.X(), p.Y(), p.Z());
    }
    logger::quick_log("");
}

template<size_t Y, size_t Z>
void BoundaryNullifier<Y, Z>::order_boundary(std::vector<LineBoundary>& bounds) const{
    for(size_t i = 0; i < bounds.size()-1; ++i){
        auto li = bounds[i];
        for(size_t j = i+1; j < bounds.size(); ++j){
            auto lj = bounds[j];
            if(li.is_connected_to(lj)){
                if(li.edges[1] != lj.edges[0]){
                    if(li.edges[1] == lj.edges[1]){
                        std::swap(lj.edges[0], lj.edges[1]);
                    } else {
                        std::swap(li.edges[0], li.edges[1]);
                        if(li.edges[1] != lj.edges[0]){
                            std::swap(lj.edges[0], lj.edges[1]);
                        }
                    }
                }
                std::swap(bounds[i+1], bounds[j]);
                break;
            }
        }
    }
}

}

#endif
