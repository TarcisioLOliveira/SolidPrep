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

#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <vector>
#include <gp_Dir.hxx>

class Material{
    public:
    enum Type{
        NONE,
        LINEAR_ELASTIC_ISOTROPIC,
        LINEAR_ELASTIC_ORTHOTROPIC
    };


    /**
     * Initializes a basic material function with only maximum design stress
     * values. Accounts for anisotropy (both for positive and negative stress
     * values) if so desired.
     *
     * @param Smax Maximum normal stresses. Input either 1 or 3 values.
     * @param Tmax Maximum shear stresses. Input either 1 or 3 values.
     */
    Material(std::vector<double> Smax, std::vector<double> Tmax);

    virtual std::vector<double> stiffness_2D() const = 0;
    virtual std::vector<double> stiffness_3D() const = 0;

    virtual double beam_E_2D(gp_Dir d = gp_Dir(1,0,0)) const = 0;
    virtual double beam_E_3D(gp_Dir d = gp_Dir(1,0,0)) const = 0;

    virtual Type get_type() const{ return this->NONE; }

    /**
     * Considers each beam node as representing maximum stresses at its cross
     * section. Therefore, returns stresses acting on the plane with normal
     * d.
     *
     * @param d Cross section normal
     * @returns Maximum stresses at the cross section plane
     */
    virtual gp_Vec get_max_stresses(gp_Dir d) const = 0;

    virtual double get_max_Von_Mises_2D() const;
    virtual double get_max_Von_Mises_3D() const;

    protected:
    gp_Mat max_stress;
    
};

#endif
