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

#include "material/linear_elastic_isotropic.hpp"
#include <cmath>

namespace material{

LinearElasticIsotropic::LinearElasticIsotropic(float E, float nu, float Smax, float Tmax, bool plane_stress):
    Material(std::vector<float>(Smax), std::vector<float>(Tmax)), E(E), nu(nu){
    
    this->D_2D.resize(9);
    if(plane_stress){
        this->D_2D[0] = E/(1-nu*nu);
        this->D_2D[1] = (E/(1-nu*nu))*nu;
        this->D_2D[3] = (E/(1-nu*nu))*nu;
        this->D_2D[4] = E/(1-nu*nu);
        this->D_2D[8] = (E/(1-nu*nu))*(1-nu)/2;
    } else {
        this->D_2D[0] = (E/((1+nu)*(1-2*nu)))*(1-nu);
        this->D_2D[1] = (E/((1+nu)*(1-2*nu)))*(nu);
        this->D_2D[3] = (E/((1+nu)*(1-2*nu)))*(nu);
        this->D_2D[4] = (E/((1+nu)*(1-2*nu)))*(1-nu);
        this->D_2D[8] = E/(2*(1+nu));
    }

    this->D_3D.resize(36);
    this->D_3D[ 0] = ((E*nu)/((1+nu)*(1-2*nu)))*(1-nu);
    this->D_3D[ 1] = ((E*nu)/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D[ 2] = ((E*nu)/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D[ 6] = ((E*nu)/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D[ 7] = ((E*nu)/((1+nu)*(1-2*nu)))*(1-nu);
    this->D_3D[ 8] = ((E*nu)/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D[12] = ((E*nu)/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D[13] = ((E*nu)/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D[14] = ((E*nu)/((1+nu)*(1-2*nu)))*(1-nu);
    this->D_3D[21] = (E*nu)/(2*(1+nu));
    this->D_3D[28] = (E*nu)/(2*(1+nu));
    this->D_3D[35] = (E*nu)/(2*(1+nu));
}

std::vector<float> LinearElasticIsotropic::stiffness_2D() const{
    return this->D_2D;
}
std::vector<float> LinearElasticIsotropic::stiffness_3D() const{
    return this->D_3D;
}

float LinearElasticIsotropic::beam_E_2D(gp_Dir d) const{
    (void)d;
    return this->E;
}
float LinearElasticIsotropic::beam_E_3D(gp_Dir d) const{
    (void)d;
    return this->E;
}

std::vector<float> LinearElasticIsotropic::get_max_stresses(gp_Dir d) const{
    (void)d;
    return {this->Smax[0], this->Smax[3], this->Tmax[0]};
}

}
