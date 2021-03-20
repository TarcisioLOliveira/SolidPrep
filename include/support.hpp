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
#ifndef SUPPORT_HPP
#define SUPPORT_HPP

#include <memory>
#include <vector>
#include "BRepClass3d_SolidClassifier.hxx"

class Support{
    public:

    /**
     * Creates a Support object for 2D problems.
     *
     * @param X Whether it supports the X axis.
     * @param Y Whether it supports the Y axis.
     * @param vertices Vertices of the support geometry.
     */
    Support(bool X, bool Y, double thickness, std::vector<std::array<double, 2>> vertices);
    /**
     * Creates a Support object for 3D problems.
     *
     * @param X Whether it supports the X axis.
     * @param Y Whether it supports the Y axis.
     * @param Z Whether it supports the Z axis.
     * @param vertices Vertices of the support geometry.
     */
    Support(bool X, bool Y, bool Z, std::vector<std::array<double, 3>> vertices);

    bool is_inside(gp_Pnt p) const;
    double get_distance(gp_Pnt p) const;
    inline TopoDS_Shape get_shape(){
        return this->shape;
    }

    private:
    TopoDS_Shape shape;
    bool X, Y, Z;
};

#endif
