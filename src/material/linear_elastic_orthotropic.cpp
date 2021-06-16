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

LinearElasticOrthotropic::LinearElasticOrthotropic(std::vector<float> E, std::vector<float> nu, std::vector<float> G, std::vector<float> Smax, std::vector<float> Tmax):
    Material(std::move(Smax), std::move(Tmax)){
   
    std::vector<float> S_2D(9); 
    S_2D[0] = 1/E[0];
    S_2D[1] = -nu[0]/E[1];
    S_2D[3] = -nu[0]/E[1];
    S_2D[4] = 1/E[1];
    S_2D[8] = 1/G[0];
    std::vector<int> ipiv(9);
    LAPACKE_sgetrf(LAPACK_ROW_MAJOR, 3, 3, S_2D.data(), 3, ipiv.data());
    LAPACKE_sgetri(LAPACK_ROW_MAJOR, 3, S_2D.data(), 3, ipiv.data());
    this->D_2D = std::move(S_2D);

    std::vector<float> S_3D(36); 
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
    LAPACKE_sgetrf(LAPACK_ROW_MAJOR, 6, 6, S_3D.data(), 6, ipiv.data());
    LAPACKE_sgetri(LAPACK_ROW_MAJOR, 6, S_3D.data(), 6, ipiv.data());
    this->D_3D = std::move(S_3D);
}

std::vector<float> LinearElasticOrthotropic::stiffness_2D() const{
    return this->D_2D;
}
std::vector<float> LinearElasticOrthotropic::stiffness_3D() const{
    return this->D_3D;
}

float LinearElasticOrthotropic::beam_E_2D(gp_Dir d) const{
    // (Hyer, 2009), Stress analysis of fiber-reinforced composites
    return std::pow(d.X(),4)*D_2D[0] + (D_2D[1] + 2*D_2D[8])*std::pow(d.X(),2)*std::pow(d.Y(),2) + std::pow(d.Y(), 4)*D_2D[4];
}
float LinearElasticOrthotropic::beam_E_3D(gp_Dir d) const{
    // (Zhao, Song, Liu, 2016)
    gp_Dir z(0,0,1);
    gp_Dir x(1,0,0);
    if(d.IsEqual(x, 0.001)){
        return this->D_3D[0];
    }
    float a = d.AngleWithRef(x, z);
    gp_Dir cross(d.Crossed(z));
    float b = M_PI/2 + d.AngleWithRef(z, cross);
    float cz = std::cos(b);
    float sz = std::sin(b);
    float cx = std::cos(a);
    float sx = std::sin(a);
    float r_x[6] = {cz*cz, cx*cx*sz*sz, sz*sz*sx*sx, -2*sx*cx*sz*sz, 2*sz*cz*sx, -2*sz*cz*cx};
    std::vector<float> out(6);

    cblas_sgemv(CblasRowMajor, CblasNoTrans, 6, 6, 1, this->D_3D.data(), 6, r_x, 1, 0, out.data(), 1);
    float Ex = cblas_sdot(6, r_x, 1, out.data(), 1);
    
    return Ex;
}

std::vector<float> LinearElasticOrthotropic::get_max_stresses(gp_Dir d) const{
    // Principles of Composite Material Mechanics (Gibson, 2016)
    // (adapted from Tsai-Wu criteria)
    // May have problems, as it is adapted from the calculation of principal
    // stress for normal and shear stresses separately. Still, beam sizing was
    // made assuming orthotropic materials, so it is not surprising it would
    // be a bit hard to adapt one to the other.
    //
    // As sizing method looks for the greatest dimensions among the ones
    // calculated, it *should* be fine most of the time (that is, most of the
    // geometry). If it is to fail somewhere, it would probably be somewhere
    // where shear and normal stress are close (that is, close to the place
    // that the load is applied) and the geometry is not aligned to the axes.
    // If it proves to be a problem, summing the greatest of the normal and
    // shear dimensions obtained should probably be enough. This problem may
    // be less noticeable in 3D problems though, as the geometry is not as 
    // sensitive to the loads in those cases.

    if(d.Z() == 0){
        gp_Dir z(0,0,1);
        gp_Dir x(1,0,0);
        float alpha = d.AngleWithRef(x, z);
        float s1 = std::pow(std::cos(alpha), 2);
        float s2 = std::pow(std::sin(alpha), 2);
        float s6 = -std::cos(alpha)*std::sin(alpha);

        float a = std::pow(s1, 2)/(this->Smax[0]*this->Smax[3])
                   + std::pow(s2, 2)/(this->Smax[1]*this->Smax[4])
                   + std::pow(s6, 2)/(this->Tmax[0]*this->Tmax[0])
                   + s1*s2/(-2*std::sqrt(this->Smax[0]*this->Smax[3]*this->Smax[1]*this->Smax[4]));
        float b = s1*(1/this->Smax[0] - 1/this->Smax[3])
                   + s2*(1/this->Smax[1] - 1/this->Smax[4]);
        float c = -1;

        float St = std::abs((-b + std::sqrt(b*b - 4*a*c))/(2*a));
        float Sc = std::abs((-b - std::sqrt(b*b - 4*a*c))/(2*a));

        s1 = 2*std::cos(alpha)*std::sin(alpha);
        s2 = -2*std::cos(alpha)*std::sin(alpha);
        s6 = std::pow(std::cos(alpha), 2) - std::pow(std::sin(alpha), 2);

        a = std::pow(s1, 2)/(this->Smax[0]*this->Smax[3])
            + std::pow(s2, 2)/(this->Smax[1]*this->Smax[4])
            + std::pow(s6, 2)/(this->Tmax[0]*this->Tmax[0])
            + s1*s2/(-2*std::sqrt(this->Smax[0]*this->Smax[3]*this->Smax[1]*this->Smax[4]));
        b = s1*(1/this->Smax[0] - 1/this->Smax[3])
            + s2*(1/this->Smax[1] - 1/this->Smax[4]);
        c = -1;

        float T1 = std::abs((-b + std::sqrt(b*b - 4*a*c))/(2*a));
        float T2 = std::abs((-b - std::sqrt(b*b - 4*a*c))/(2*a));
        float T = std::min(T1, T2);

        return {St, Sc, T};
    } else {
        // TODO 3D implementation
        // Need to find how to make that stress to principal stress transfor-
        // mation for 3D.
    }

    return std::vector<float>();
}

}
