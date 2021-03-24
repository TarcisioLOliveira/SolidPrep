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

#include "internal_force.hpp"


InternalForce::InternalForce(const Force& f):
    point(f.get_centroid()), section_normal(f.get_normal()), force(f.get_force()), mass(f.get_mass()),
    variable(false), before(), after(){
}
InternalForce::InternalForce(double mass, gp_Pnt point, gp_Dir normal):
    point(point), section_normal(normal), force(), mass(mass),
    variable(true), before(), after(){
}

InternalForce InternalForce::merge(const InternalForce& f1, const InternalForce& f2){
    double mass = f1.mass + f2.mass;
    gp_Pnt point = f1.point;
    point.BaryCenter(f1.mass, f2.point, f2.mass);
    gp_Vec vnormal(0, 0, 0);
    for(auto i:f1.before){
        vnormal += gp_Vec(i->point, point);
    }
    for(auto i:f2.before){
        vnormal += gp_Vec(i->point, point);
    }
    for(auto i:f1.after){
        vnormal += gp_Vec(point, i->point);
    }
    for(auto i:f2.after){
        vnormal += gp_Vec(point, i->point);
    }
    gp_Dir normal(vnormal);

    InternalForce f(mass, point, normal);
    for(auto i:f1.before){
        f.add_before(i);
    }
    for(auto i:f2.before){
        f.add_before(i);
    }
    for(auto i:f1.after){
        f.add_after(i);
    }
    for(auto i:f2.after){
        f.add_after(i);
    }

    return f;
}
