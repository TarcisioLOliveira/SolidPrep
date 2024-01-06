/*
 *   Copyright (C) 2023 Tarcísio Ladeia de Oliveira.
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

#include <lapacke.h>
#include "logger.hpp"
#include "utils/basis_tensor.hpp"
#include "utils/D_operations.hpp"
#include "material/linear_elastic_orthotropic_field.hpp"

namespace material{

LinearElasticOrthotropicField::LinearElasticOrthotropicField(const std::string& name, const double density, std::vector<double> E, std::vector<double> nu, std::vector<double> G, std::vector<double> Smax, std::vector<double> Tmax, const CoordinateField* field):
    Material(name, std::move(Smax), std::move(Tmax)), density(density), field(field){

    logger::log_assert(field->get_sub_type() == Field::SubType::DOMAIN || field->get_sub_type() == Field::SubType::PROJECTION, logger::ERROR, "field subtype for orthotropic field material must be DOMAIN or PROJECTION");
   
    this->S_2D.resize(9,0);
    S_2D[0] = 1/E[0]; S_2D[1] = -nu[0]/E[0];
    S_2D[3] = -nu[0]/E[0]; S_2D[4] = 1/E[1];
    S_2D[4] = 1/E[1];
    S_2D[8] = 1/G[0];
    this->D_2D = utils::D_op::invert_2D(S_2D);

    this->S_3D.resize(36,0);
    S_3D[ 0] = 1/E[0]; S_3D[ 1] = -nu[0]/E[0]; S_3D[ 2] = -nu[1]/E[0];
    S_3D[ 6] = -nu[0]/E[0]; S_3D[ 7] = 1/E[1]; S_3D[ 8] = -nu[2]/E[1];
    S_3D[12] = -nu[1]/E[0]; S_3D[13] = -nu[2]/E[1]; S_3D[14] = 1/E[2];
    S_3D[21] = 1/G[0];
    S_3D[28] = 1/G[1];
    S_3D[35] = 1/G[2];
    this->D_3D = utils::D_op::invert_3D(S_3D);
}

std::vector<double> LinearElasticOrthotropicField::stiffness_2D(const MeshElement* const e, const gp_Pnt& p) const{
    const auto M = this->field->get_matrix(e, p);
    const auto R = utils::basis_tensor_2D(M({0,1},{0,1}));
    auto D = this->D_2D;
    this->rotate_D_2D(D, R);

    return D;
}
std::vector<double> LinearElasticOrthotropicField::stiffness_3D(const MeshElement* const e, const gp_Pnt& p) const{
    const auto M = this->field->get_matrix(e, p);
    const auto R = utils::basis_tensor_3D(M);
    auto D = this->D_3D;
    this->rotate_D_3D(D, R);

    return D;
}
std::vector<double> LinearElasticOrthotropicField::stiffness_inverse_2D(const MeshElement* const e, const gp_Pnt& p) const{
    const auto M = this->field->get_matrix(e, p);
    const auto R = utils::basis_tensor_2D_inv_T(M({0,1},{0,1}));
    auto S = this->S_2D;
    this->rotate_S_2D(S, R);

    return S;
}
std::vector<double> LinearElasticOrthotropicField::stiffness_inverse_3D(const MeshElement* const e, const gp_Pnt& p) const{
    const auto M = this->field->get_matrix(e, p);
    const auto R = utils::basis_tensor_3D_inv_T(M);
    auto S = this->S_3D;
    this->rotate_S_3D(S, R);

    return S;
}

double LinearElasticOrthotropicField::beam_E_2D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 2, 2>& R) const{
    // Anisotropic Elasticity: Theory and Applications
    // (Ting, 1996)

    auto d = this->stiffness_2D(e, p);
    const auto Rt = utils::basis_tensor_2D(R);
    this->rotate_D_2D(d, Rt);

    return d[0];
}
double LinearElasticOrthotropicField::beam_E_3D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 3, 3>& R) const{
    // Anisotropic Elasticity: Theory and Applications
    // (Ting, 1996)

    auto d = this->stiffness_3D(e, p);

    const auto Rt = utils::basis_tensor_3D(R);
    this->rotate_D_3D(d, Rt);

    return d[0];
}
std::array<double, 2> LinearElasticOrthotropicField::beam_EG_2D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 2, 2>& R) const{
    // Anisotropic Elasticity: Theory and Applications
    // (Ting, 1996)

    auto d = this->stiffness_2D(e, p);

    const auto Rt = utils::basis_tensor_2D(R);
    this->rotate_D_2D(d, Rt);

    return {d[0], d[8]};
}
std::array<double, 4> LinearElasticOrthotropicField::beam_EG_3D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 3, 3>& R) const{
    // Anisotropic Elasticity: Theory and Applications
    // (Ting, 1996)

    auto d = this->stiffness_3D(e, p);

    const auto Rt = utils::basis_tensor_3D(R);
    this->rotate_D_3D(d, Rt);

    return {d[0], d[21], d[28], d[35]};
}

double LinearElasticOrthotropicField::S12_2D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 2, 2>& R) const{
    auto s = this->stiffness_inverse_2D(e, p);

    const auto Rt = utils::basis_tensor_2D_inv_T(R);
    this->rotate_S_2D(s, Rt);

    return s[1];
}

std::array<double, 2> LinearElasticOrthotropicField::S12_S13_3D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 3, 3>& R) const{

    auto s = this->stiffness_inverse_3D(e, p);

    const auto Rt = utils::basis_tensor_3D_inv_T(R);
    this->rotate_S_3D(s, Rt);

    return {s[1], s[2]};
}

std::vector<double> LinearElasticOrthotropicField::get_max_stresses(gp_Dir d) const{
    // Temporarily copied from LinearElasticOrthotropic
    //
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

}
