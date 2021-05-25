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

#include "material/linear_elastic_orthotropic.hpp"
#include "logger.hpp"
#include <cmath>
#include <lapacke.h>
#include <cblas.h>

namespace material{

LinearElasticOrthotropic::LinearElasticOrthotropic(std::vector<double> E, std::vector<double> nu, std::vector<double> G, std::vector<double> Smax, std::vector<double> Tmax):
    Material(std::move(Smax), std::move(Tmax)){
   
    std::vector<double> S_2D(9); 
    S_2D[0] = 1/E[0];
    S_2D[1] = -nu[0]/E[1];
    S_2D[3] = -nu[0]/E[1];
    S_2D[4] = 1/E[1];
    S_2D[8] = 1/G[0];
    std::vector<int> ipiv(9);
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 3, 3, S_2D.data(), 3, ipiv.data());
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, 3, S_2D.data(), 3, ipiv.data());
    this->D_2D = std::move(S_2D);

    std::vector<double> S_3D(36); 
    S_3D[ 0] = 1/E[0];
    S_3D[ 1] = -nu[0]/E[1];
    S_3D[ 2] = -nu[1]/E[2];
    S_3D[ 6] = -nu[0]/E[0];
    S_3D[ 7] = 1/E[1];
    S_3D[ 8] = -nu[2]/E[2];
    S_3D[12] = -nu[1]/E[0];
    S_3D[13] = -nu[2]/E[1];
    S_3D[14] = 1/E[2];
    S_3D[21] = 1/G[0];
    S_3D[28] = 1/G[1];
    S_3D[35] = 1/G[2];
    ipiv.resize(36);
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 6, 6, S_3D.data(), 6, ipiv.data());
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, 6, S_3D.data(), 6, ipiv.data());
    this->D_3D = std::move(S_3D);
}

std::vector<double> LinearElasticOrthotropic::stiffness_2D() const{
    return this->D_2D;
}
std::vector<double> LinearElasticOrthotropic::stiffness_3D() const{
    return this->D_3D;
}

double LinearElasticOrthotropic::beam_E_2D(gp_Dir d) const{
    // (Hyer, 2009), Stress analysis of fiber-reinforced composites
    return std::pow(d.X(),4)*D_2D[0] + (D_2D[1] + 2*D_2D[8])*std::pow(d.X(),2)*std::pow(d.Y(),2) + std::pow(d.Y(), 4)*D_2D[4];
}
double LinearElasticOrthotropic::beam_E_3D(gp_Dir d) const{
    // (Zhao, Song, Liu, 2016)
    gp_Dir z(0,0,1);
    gp_Dir x(1,0,0);
    if(d.IsEqual(x, 0.001)){
        return this->D_3D[0];
    }
    double a = d.AngleWithRef(x, z);
    gp_Dir cross(d.Crossed(z));
    double b = M_PI/2 + d.AngleWithRef(z, cross);
    double cz = std::cos(b);
    double sz = std::sin(b);
    double cx = std::cos(a);
    double sx = std::sin(a);
    double r_x[6] = {cz*cz, cx*cx*sz*sz, sz*sz*sx*sx, -2*sx*cx*sz*sz, 2*sz*cz*sx, -2*sz*cz*cx};
    std::vector<double> out(6);

    cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 6, 1, this->D_3D.data(), 6, r_x, 1, 0, out.data(), 1);
    double Ex = cblas_ddot(6, r_x, 1, out.data(), 1);
    
    return Ex;
}

gp_Vec LinearElasticOrthotropic::get_max_stresses(gp_Dir d) const{
    gp_Vec result(0,0,0);
    for(int i = 1; i < 4; ++i){
        double s = 0;
        for(int j = 1; j < 4; ++j){
            s += this->max_stress(i,j)*std::abs(d.Coord(j));
        }
        result.SetCoord(i, std::abs(s));
    }

    return result;
}

}
