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
#include <cmath>
#include <lapacke.h>
#include <cblas.h>

namespace material{

LinearElasticIsotropic::LinearElasticIsotropic(const std::string& name, const double density, double E, double nu, double Smax, double Tmax, bool plane_stress):
    Material(name, {Smax}, {Tmax}), E(E), G(E/(2*(1 + nu))), nu(nu), density(density){

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
    std::vector<int> ipiv(9);
    std::vector<double> D_2D_tmp = D_2D;
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 3, 3, D_2D_tmp.data(), 3, ipiv.data());
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, 3, D_2D_tmp.data(), 3, ipiv.data());
    this->S_2D = std::move(D_2D_tmp);

    this->D_3D.resize(36);
    this->D_3D[ 0] = (E/((1+nu)*(1-2*nu)))*(1-nu);
    this->D_3D[ 1] = (E/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D[ 2] = (E/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D[ 6] = (E/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D[ 7] = (E/((1+nu)*(1-2*nu)))*(1-nu);
    this->D_3D[ 8] = (E/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D[12] = (E/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D[13] = (E/((1+nu)*(1-2*nu)))*(nu);
    this->D_3D[14] = (E/((1+nu)*(1-2*nu)))*(1-nu);
    this->D_3D[21] = E/(2*(1+nu));
    this->D_3D[28] = E/(2*(1+nu));
    this->D_3D[35] = E/(2*(1+nu));

    ipiv.resize(36);
    std::vector<double> D_3D_tmp = D_3D;
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 6, 6, D_3D_tmp.data(), 6, ipiv.data());
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, 6, D_3D_tmp.data(), 6, ipiv.data());
    this->S_3D = std::move(D_3D_tmp);
}

double LinearElasticIsotropic::beam_E_2D(const gp_Pnt& p, gp_Dir d) const{
    (void)d;
    (void)p;
    return this->E;
}
double LinearElasticIsotropic::beam_E_3D(const gp_Pnt& p, gp_Dir d) const{
    (void)d;
    (void)p;
    return this->E;
}
std::array<double, 2> LinearElasticIsotropic::beam_EG_2D(const gp_Pnt& p, gp_Dir d) const{
    (void)d;
    (void)p;
    return {E, G};
}
std::array<double, 4> LinearElasticIsotropic::beam_EG_3D(const gp_Pnt& p, gp_Dir d) const{
    (void)d;
    (void)p;
    return {E, G, G, G};
}

std::vector<double> LinearElasticIsotropic::get_max_stresses(gp_Dir d) const{
    (void)d;
    return {this->Smax[0], this->Smax[0], this->Tmax[0]};
}

}
