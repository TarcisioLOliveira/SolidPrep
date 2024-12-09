/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
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

#include "material/linear_elastic_isotropic.hpp"
#include "logger.hpp"
#include <lapacke.h>
#include <cblas.h>

namespace material{

LinearElasticIsotropic::LinearElasticIsotropic(const std::string& name, const double density, double E, double nu, double Smax, double Tmax, bool plane_stress):
    Material(name, {Smax}, {Tmax}), E(E), G(E/(2*(1 + nu))), nu(nu), density(density),
    D_2D(3,3), D_3D(6,6), S_2D(), S_3D(){

    if(plane_stress){
        this->D_2D.data()[0] = E/(1-nu*nu);
        this->D_2D.data()[1] = (E/(1-nu*nu))*nu;
        this->D_2D.data()[3] = (E/(1-nu*nu))*nu;
        this->D_2D.data()[4] = E/(1-nu*nu);
        this->D_2D.data()[8] = (E/(1-nu*nu))*(1-nu)/2;
    } else {
        this->D_2D.data()[0] = (E/((1+nu)*(1-2*nu)))*(1-nu);
        this->D_2D.data()[1] = (E/((1+nu)*(1-2*nu)))*(nu);
        this->D_2D.data()[3] = (E/((1+nu)*(1-2*nu)))*(nu);
        this->D_2D.data()[4] = (E/((1+nu)*(1-2*nu)))*(1-nu);
        this->D_2D.data()[8] = E/(2*(1+nu));
    }
    this->S_2D = D_2D.get_inverted_cholesky();

    this->D_3D.data()[ 0] = (E/((1+nu)*(1-2*nu)))*(1-nu);
    this->D_3D.data()[ 1] = (E/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D.data()[ 2] = (E/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D.data()[ 6] = (E/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D.data()[ 7] = (E/((1+nu)*(1-2*nu)))*(1-nu);
    this->D_3D.data()[ 8] = (E/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D.data()[12] = (E/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D.data()[13] = (E/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D.data()[14] = (E/((1+nu)*(1-2*nu)))*(1-nu);
    this->D_3D.data()[21] = E/(2*(1+nu));
    this->D_3D.data()[28] = E/(2*(1+nu));
    this->D_3D.data()[35] = E/(2*(1+nu));

    this->S_3D = D_3D.get_inverted_cholesky();
}

std::vector<double> LinearElasticIsotropic::get_max_stresses(gp_Dir d) const{
    (void)d;
    return {this->Smax[0], this->Smax[0], this->Tmax[0]};
}

}
