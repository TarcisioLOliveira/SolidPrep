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

#ifndef CROSS_SECTION_HPP
#define CROSS_SECTION_HPP

#include "utils.hpp"
#include <vector>
#include <gp_Pnt.hxx>
#include <TopoDS_Shape.hxx>

class CrossSection{
    public:
    struct Rectangle{
        double w, h;
        gp_Pnt center;
        gp_Vec rotation;
    };
    struct Circle{
        double r;
        gp_Pnt center;
        gp_Vec rotation;
    };

    /**
     * Creates a cross section for 2D problems.
     *
     * @param vertices Vertices of the application geometry.
     * @param thickness Thickness of the plate.
     */
    CrossSection(std::vector<gp_Pnt> vertices, double thickness);

    /**
     * Creates a general cross section for 3D problems.
     *
     * @param vertices Vertices of the application geometry.
     */
    CrossSection(std::vector<gp_Pnt> vertices);

    /**
     * Creates a rectangular cross section for 3D problems.
     *
     * @param r Rectangle object.
     */
    CrossSection(Rectangle r);

    /**
     * Creates a circular cross section for 3D problems.
     *
     * @param c Circular object.
     */
    CrossSection(Circle c);

    /**
     * Creates a dummy point cross section, for pathfinding for example.
     *
     * @param p Point.
     */
    CrossSection(gp_Pnt p);

    /**
     * Creates a dummy cross section with size, for pathfinding for example.
     *
     * @param p Point.
     */
    CrossSection(gp_Pnt p, utils::ProblemType type, double radius = 10);

    bool is_inside(gp_Pnt p) const;
    double get_distance(gp_Pnt p) const;

    inline double get_moment_of_inertia(int i, int j) const{
        return this->inertia(i, j);
    }
    inline gp_Mat get_moment_of_inertia() const{
        return this->inertia;
    }
    inline gp_Pnt get_centroid() const{
        return this->centroid;
    }
    inline gp_Dir get_normal() const{
        return this->normal;
    }
    inline double get_dimension() const{
        return this->max_dim;
    }
    inline double get_area() const{
        return this->area;
    }
    inline TopoDS_Shape get_shape() const{
        return this->shape;
    }

    inline void set_normal(gp_Dir n){
        this->normal = n;
    }
    void set_centroid(gp_Pnt p);

    private:
    gp_Pnt centroid;
    gp_Mat inertia;
    gp_Dir normal;
    Standard_Real max_dim;
    TopoDS_Shape shape;
    double area;

};

#endif
