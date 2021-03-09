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

#ifndef FORCE_HPP
#define FORCE_HPP
#include "BRepClass3d_SolidClassifier.hxx"
#include <vector>

class Force{
    public:

    // MUST CHECK IF SHAPES ARE BEING CREATED PROPERLY.

    /**
     * Creates a Force object for 2D problems.
     *
     * @param vertices Vertices of the application geometry.
     * @param thickness Thickness of the plate.
     * @param force Force vector.
     */
    Force(std::vector<std::array<double, 2>> vertices, double thickness, std::array<double, 2> force);

    /**
     * Creates a Force object for 3D problems.
     *
     * @param vertices Vertices of the application geometry.
     * @param force Force vector.
     */
    Force(std::vector<std::array<double, 3>> vertices,  std::array<double, 3> force);

    bool is_inside(gp_Pnt p) const;
    double get_moment_of_inertia(int i, int j);


    private:
    gp_Pnt centroid;
    gp_Mat inertia;
    gp_Dir normal;
    gp_Vec force;
    Standard_Real max_dim;
    TopoDS_Shape shape;
};

#endif
