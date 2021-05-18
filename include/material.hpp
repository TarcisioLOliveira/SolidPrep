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
    /**
     * Initializes a basic material function with only maximum design stress
     * values. Accounts for anisotropy (both for positive and negative stress
     * values) if so desired.
     *
     * @param Smax Maximum normal stresses. Input either 1, 3, or 6 values.
     * @param Tmax Maximum shear stresses. Input either 1, 3, or 6 values.
     */
    Material(std::vector<double> Smax, std::vector<double> Tmax);

    virtual std::vector<double> stiffness_2D() const = 0;
    virtual std::vector<double> stiffness_3D() const = 0;

    virtual double beam_E_2D(gp_Dir d) const = 0;
    virtual double beam_E_3D(gp_Dir d) const = 0;

    inline double get_Smax(size_t i = 0) const{return this->Smax[i];}
    inline double get_Tmax(size_t i = 0) const{return this->Tmax[i];}

    private:
    std::vector<double> Smax;
    std::vector<double> Tmax;
};

class LinearElasticOrthotropic : public Material{

};

#endif
