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

#include "material/linear_elastic_orthotropic.hpp"
#include "logger.hpp"
#include <cmath>
#include <lapacke.h>
#include <cblas.h>

namespace material{

LinearElasticOrthotropic::LinearElasticOrthotropic(const std::string& name, const double density, std::vector<double> E, std::vector<double> nu, std::vector<double> G, std::vector<double> Smax, std::vector<double> Tmax):
    Material(name, std::move(Smax), std::move(Tmax)), density(density){
   
    this->S_2D.resize(9);
    S_2D[0] = 1/E[0];
    S_2D[1] = -nu[0]/E[1];
    S_2D[3] = -nu[0]/E[1];
    S_2D[4] = 1/E[1];
    S_2D[8] = 1/G[0];
    std::vector<int> ipiv(9);
    std::vector<double> S_2D_tmp = S_2D;
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 3, 3, S_2D_tmp.data(), 3, ipiv.data());
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, 3, S_2D_tmp.data(), 3, ipiv.data());
    this->D_2D = std::move(S_2D_tmp);

    this->S_3D.resize(36);
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
    std::vector<double> S_3D_tmp = S_3D;
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 6, 6, S_3D_tmp.data(), 6, ipiv.data());
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, 6, S_3D_tmp.data(), 6, ipiv.data());
    this->D_3D = std::move(S_3D_tmp);
}

double LinearElasticOrthotropic::beam_E_2D(const MeshElement* const e, const gp_Pnt& p, gp_Dir dir) const{
    // Anisotropic Elasticity: Theory and Applications
    // (Ting, 1996)
    (void)p;

    const double ca = dir.X(); 
    const double sa = dir.Y(); 
    const auto& d = this->D_2D;

    const double E = ca*ca*ca*ca*d[0] + 2*ca*ca*d[1]*sa*sa + ca*ca*d[8]*sa*sa + d[4]*sa*sa*sa*sa;

    return E;
}
double LinearElasticOrthotropic::beam_E_3D(const MeshElement* const e, const gp_Pnt& p, gp_Dir dir) const{
    // Anisotropic Elasticity: Theory and Applications
    // (Ting, 1996)
    (void)p;

    const auto& d = this->D_3D;

    gp_Dir z(0,0,1);
    gp_Dir x(1,0,0);
    if(dir.IsEqual(x, 0.001)){
        return this->D_3D[0];
    }
    double a = dir.AngleWithRef(x, z);
    gp_Dir cross(dir.Crossed(z));
    double b = M_PI/2 + dir.AngleWithRef(z, cross);
    double cb = std::cos(b);
    double sb = std::sin(b);
    double ca = std::cos(a);
    double sa = std::sin(a);

   const double E = ca*ca*ca*ca*d[7]*sb*sb*sb*sb + 2*ca*ca*cb*cb*d[1]*sb*sb - ca*ca*cb*cb*d[21]*sb*sb - ca*ca*d[35]*sa*sa*sb*sb*sb*sb + 2*ca*ca*d[8]*sa*sa*sb*sb*sb*sb + cb*cb*cb*cb*d[0] + 2*cb*cb*d[2]*sa*sa*sb*sb - cb*cb*d[28]*sa*sa*sb*sb + d[14]*sa*sa*sa*sa*sb*sb*sb*sb;
    
    return E;
}
std::array<double, 2> LinearElasticOrthotropic::beam_EG_2D(const MeshElement* const e, const gp_Pnt& p, gp_Dir dir) const{
    // Anisotropic Elasticity: Theory and Applications
    // (Ting, 1996)
    (void)p;

    const double ca = dir.X(); 
    const double sa = dir.Y(); 
    const auto& d = this->D_2D;

    const double E = ca*ca*ca*ca*d[0] + 2*ca*ca*d[1]*sa*sa + ca*ca*d[8]*sa*sa + d[4]*sa*sa*sa*sa;
    const double G = ca*ca*ca*ca*d[8] + ca*ca*d[0]*sa*sa - 2*ca*ca*d[1]*sa*sa + ca*ca*d[4]*sa*sa - 2*ca*ca*d[8]*sa*sa + d[8]*sa*sa*sa*sa;

    return {E, G};
}
std::array<double, 4> LinearElasticOrthotropic::beam_EG_3D(const MeshElement* const e, const gp_Pnt& p, gp_Dir dir) const{
    // Anisotropic Elasticity: Theory and Applications
    // (Ting, 1996)
    (void)p;

    const auto& d = this->D_3D;

    gp_Dir z(0,0,1);
    gp_Dir x(1,0,0);
    if(dir.IsEqual(x, 0.001)){
        return {D_3D[0], D_3D[21], D_3D[28], D_3D[35]};
    }
    double a = dir.AngleWithRef(x, z);
    gp_Dir cross(dir.Crossed(z));
    double b = M_PI/2 + dir.AngleWithRef(z, cross);
    double cb = std::cos(b);
    double sb = std::sin(b);
    double ca = std::cos(a);
    double sa = std::sin(a);

    const double E = ca*ca*ca*ca*d[7]*sb*sb*sb*sb + 2*ca*ca*cb*cb*d[1]*sb*sb + ca*ca*cb*cb*d[21]*sb*sb + ca*ca*d[35]*sa*sa*sb*sb*sb*sb + 2*ca*ca*d[8]*sa*sa*sb*sb*sb*sb + cb*cb*cb*cb*d[0] + 2*cb*cb*d[2]*sa*sa*sb*sb + cb*cb*d[28]*sa*sa*sb*sb + d[14]*sa*sa*sa*sa*sb*sb*sb*sb;
    const double G1 = ca*ca*ca*ca*cb*cb*d[7]*sb*sb + ca*ca*cb*cb*cb*cb*d[21] - 2*ca*ca*cb*cb*d[1]*sb*sb - 2*ca*ca*cb*cb*d[21]*sb*sb + ca*ca*cb*cb*d[35]*sa*sa*sb*sb + 2*ca*ca*cb*cb*d[8]*sa*sa*sb*sb + ca*ca*d[21]*sb*sb*sb*sb + cb*cb*cb*cb*d[28]*sa*sa + cb*cb*d[0]*sb*sb + cb*cb*d[14]*sa*sa*sa*sa*sb*sb - 2*cb*cb*d[2]*sa*sa*sb*sb - 2*cb*cb*d[28]*sa*sa*sb*sb + d[28]*sa*sa*sb*sb*sb*sb;
    const double G2 = ca*ca*ca*ca*d[35]*sb*sb + ca*ca*cb*cb*d[28] + ca*ca*d[14]*sa*sa*sb*sb - 2*ca*ca*d[35]*sa*sa*sb*sb + ca*ca*d[7]*sa*sa*sb*sb - 2*ca*ca*d[8]*sa*sa*sb*sb + cb*cb*d[21]*sa*sa + d[35]*sa*sa*sa*sa*sb*sb;
    const double G3 = ca*ca*ca*ca*cb*cb*d[35] + ca*ca*cb*cb*d[14]*sa*sa - 2*ca*ca*cb*cb*d[35]*sa*sa + ca*ca*cb*cb*d[7]*sa*sa - 2*ca*ca*cb*cb*d[8]*sa*sa + ca*ca*d[28]*sb*sb + cb*cb*d[35]*sa*sa*sa*sa + d[21]*sa*sa*sb*sb;

    return {E, G1, G2, G3};
}

std::vector<double> LinearElasticOrthotropic::get_max_stresses(gp_Dir d) const{
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
        double alpha = d.AngleWithRef(x, z);
        double s1 = std::pow(std::cos(alpha), 2);
        double s2 = std::pow(std::sin(alpha), 2);
        double s6 = -std::cos(alpha)*std::sin(alpha);

        double a = std::pow(s1, 2)/(this->Smax[0]*this->Smax[3])
                   + std::pow(s2, 2)/(this->Smax[1]*this->Smax[4])
                   + std::pow(s6, 2)/(this->Tmax[0]*this->Tmax[0])
                   + s1*s2/(-2*std::sqrt(this->Smax[0]*this->Smax[3]*this->Smax[1]*this->Smax[4]));
        double b = s1*(1/this->Smax[0] - 1/this->Smax[3])
                   + s2*(1/this->Smax[1] - 1/this->Smax[4]);
        double c = -1;

        double St = std::abs((-b + std::sqrt(b*b - 4*a*c))/(2*a));
        double Sc = std::abs((-b - std::sqrt(b*b - 4*a*c))/(2*a));

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

        double T1 = std::abs((-b + std::sqrt(b*b - 4*a*c))/(2*a));
        double T2 = std::abs((-b - std::sqrt(b*b - 4*a*c))/(2*a));
        double T = std::min(T1, T2);

        return {St, Sc, T};
    } else {
        // TODO 3D implementation
        // Need to find how to make that stress to principal stress transfor-
        // mation for 3D.
    }

    return std::vector<double>();
}


double LinearElasticOrthotropic::S12_2D(const MeshElement* const e, const gp_Pnt& p, gp_Dir dir) const{
    (void)p;
    const double ca = dir.X(); 
    const double sa = dir.Y(); 
    const auto& s = this->S_2D;

    const double S12 = ca*ca*ca*ca*s[1] + ca*ca*s[0]*sa*sa + ca*ca*s[4]*sa*sa - ca*ca*s[8]*sa*sa + s[1]*sa*sa*sa*sa;

    return S12;
}

std::array<double, 2> LinearElasticOrthotropic::S12_S13_3D(const MeshElement* const e, const gp_Pnt& p, gp_Dir dir) const{
    (void)p;
    const auto& s = this->S_3D;

    gp_Dir z(0,0,1);
    gp_Dir x(1,0,0);
    if(dir.IsEqual(x, 0.001)){
        return {S_3D[1], S_3D[2]};
    }
    double a = dir.AngleWithRef(x, z);
    gp_Dir cross(dir.Crossed(z));
    double b = M_PI/2 + dir.AngleWithRef(z, cross);
    double cb = std::cos(b);
    double sb = std::sin(b);
    double ca = std::cos(a);
    double sa = std::sin(a);

    const double S12 = ca*ca*ca*ca*cb*cb*s[7]*sb*sb + ca*ca*cb*cb*cb*cb*s[1] - ca*ca*cb*cb*s[21]*sb*sb + ca*ca*cb*cb*s[35]*sa*sa*sb*sb + 2*ca*ca*cb*cb*s[8]*sa*sa*sb*sb + ca*ca*s[1]*sb*sb*sb*sb + cb*cb*cb*cb*s[2]*sa*sa + cb*cb*s[0]*sb*sb + cb*cb*s[14]*sa*sa*sa*sa*sb*sb - cb*cb*s[28]*sa*sa*sb*sb + s[2]*sa*sa*sb*sb*sb*sb;
    const double S13 = ca*ca*ca*ca*s[8]*sb*sb + ca*ca*cb*cb*s[2] + ca*ca*s[14]*sa*sa*sb*sb - ca*ca*s[35]*sa*sa*sb*sb + ca*ca*s[7]*sa*sa*sb*sb + cb*cb*s[1]*sa*sa + s[8]*sa*sa*sa*sa*sb*sb;

    return {S12, S13};
}

}
