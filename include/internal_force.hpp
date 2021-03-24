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

#ifndef INTERNAL_FORCE_HPP
#define INTERNAL_FORCE_HPP

#include <gp_Pnt.hxx>
#include <vector>
#include "force.hpp"

class InternalForce{
    public:
    InternalForce(const Force& f);
    InternalForce(double mass, gp_Pnt point, gp_Dir normal);

    InternalForce merge(const InternalForce& f1, const InternalForce& f2);

    inline void add_before(InternalForce* f){
        this->before.push_back(f);
    }
    inline void add_after(InternalForce* f){
        this->after.push_back(f);
    }

    inline void set_force(gp_Vec force){
        this->force = force;
    }
    inline gp_Vec get_force(){
        return this->force;
    }
    inline gp_Dir get_normal(){
        return this->section_normal;
    }
    inline gp_Pnt get_point(){
        return this->point;
    }
    inline double get_mass(){
        return this->mass;
    }
    inline bool is_variable(){
        return this->variable;
    }

    private:
    gp_Pnt point;
    gp_Dir section_normal;
    gp_Vec force;
    double mass;
    bool variable;
    std::vector<InternalForce*> before;
    std::vector<InternalForce*> after;
};

#endif
