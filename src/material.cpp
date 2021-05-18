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

#include "material.hpp"
#include <cmath>

Material::Material(std::vector<double> Smax, std::vector<double> Tmax){
    if(Smax.size() == 1){
        Smax = std::vector<double>(3, Smax[0]);
    }
    if(Tmax.size() == 1){
        Tmax = std::vector<double>(3, Tmax[0]);
    }
    this->max_stress = gp_Mat(Smax[0], Tmax[0], Tmax[1],
                              Tmax[0], Smax[1], Tmax[2],
                              Tmax[1], Tmax[2], Smax[2]);
}

gp_Mat Material::get_max_stresses_2D(gp_Dir d) const{
    gp_Dir z(0,0,1);
    gp_Dir x(1,0,0);
    double a = d.AngleWithRef(x, z);
    gp_Mat rot;
    rot.SetRotation(z.XYZ(), a);
    gp_Mat result = rot.Transposed()*this->max_stress*rot;

    return result;
}

gp_Mat Material::get_max_stresses_3D(gp_Dir d) const{
    gp_Dir z(0,0,1);
    gp_Dir x(1,0,0);
    double a = d.AngleWithRef(x, z);
    gp_Dir cross(d.Crossed(z));
    double b = M_PI/2 + d.AngleWithRef(z, cross);

    gp_Mat rot1;
    rot1.SetRotation(z.XYZ(), a);
    gp_Mat rot2;
    rot1.SetRotation(cross.XYZ(), b);
    gp_Mat rot(rot1*rot2);
    gp_Mat result = rot.Transposed()*this->max_stress*rot;

    return result;
}

double Material::get_max_Von_Mises_2D() const{
    return std::sqrt(0.5*(std::pow(max_stress(1,1)-max_stress(2,2), 2)
                          + std::pow(max_stress(1,1), 2) + std::pow(max_stress(2,2),2)
                          + 6*std::pow(max_stress(1,2), 2)));
}
double Material::get_max_Von_Mises_3D() const{
    return std::sqrt(0.5*(std::pow(max_stress(1,1)-max_stress(2,2), 2)
                          + std::pow(max_stress(2,2)-max_stress(3,3), 2)
                          + std::pow(max_stress(1,1)-max_stress(3,3), 2)
                          + 6*(std::pow(max_stress(1,2), 2)
                               + std::pow(max_stress(1,3), 2)
                               + std::pow(max_stress(2,3), 2))));
}
