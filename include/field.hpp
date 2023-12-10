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

#ifndef FIELD_HPP
#define FIELD_HPP

#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include "element.hpp"
#include "visualization.hpp"

class Field{
    public:
    enum class Type{
        SCALAR,
        VECTOR,
        DIRECTION,
        COORDINATE
    };
    enum class SubType{
        DOMAIN,       // Limited to the domain of the union of specified geometries
        PROJECTION,   // Limited to the volumetric projection of a cross-section
        BOUNDARY      // Limited to the boundary of the union of geometries
    };
    virtual void generate() = 0;
    virtual void initialize_views(Visualization* viz) = 0;
    virtual void display_views() const = 0;
    virtual Type get_type() const = 0;
    virtual SubType get_sub_type() const = 0;
};

class ScalarField : public Field{
    public:
    inline virtual Type get_type() const override{
        return Type::SCALAR;
    };
    virtual double get(const MeshElement* e, const gp_Pnt& p) const = 0;
};

class VectorField : public Field{
    public:
    inline virtual Type get_type() const override{
        return Type::VECTOR;
    };
    virtual gp_Vec get(const MeshElement* e, const gp_Pnt& p) const = 0;
};

class DirectionField : public Field{
    public:
    inline virtual Type get_type() const override{
        return Type::DIRECTION;
    };
    virtual gp_Dir get(const MeshElement* e, const gp_Pnt& p) const = 0;
};

class CoordinateField : public Field{
    public:
    enum class Class{
        ORTHOTROPIC_FLOW
    };
    inline virtual Type get_type() const override{
        return Type::COORDINATE;
    };
    virtual std::array<gp_Dir, 3> get_array(const MeshElement* e, const gp_Pnt& p) const = 0;
    virtual Eigen::Matrix<double, 3, 3> get_matrix(const MeshElement* e, const gp_Pnt& p) const = 0;
    virtual Class get_class() const = 0;
};

#endif
