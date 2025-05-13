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
#include "material/linear_elastic_orthotropic_field.hpp"
#include "project_data.hpp"

namespace material{

LinearElasticOrthotropicField::LinearElasticOrthotropicField(const std::string& name, const double density, std::vector<double> E, std::vector<double> nu, std::vector<bool> nu_lower_half, std::vector<double> G, std::vector<double> Smax, std::vector<double> Tmax, CoordinateField* field):
    Material(name, std::move(Smax), std::move(Tmax)), density(density), field(field),
    D_2D(), D_3D(), S_2D(3,3), S_3D(6,6){

    logger::log_assert(field->get_sub_type() == Field::SubType::DOMAIN || field->get_sub_type() == Field::SubType::PROJECTION, logger::ERROR, "field subtype for orthotropic field material must be DOMAIN or PROJECTION");

    double Sxy = 0, Sxz = 0, Syz = 0;
    if(nu_lower_half[0]){
        Sxy = -nu[0]/E[0];
    } else {
        Sxy = -nu[0]/E[1];
    }
    if(nu_lower_half[1]){
        Sxz = -nu[1]/E[0];
    } else {
        Sxz = -nu[1]/E[2];
    }
    if(nu_lower_half[2]){
        Syz = -nu[2]/E[1];
    } else {
        Syz = -nu[2]/E[2];
    }

    S_2D.data()[0] = 1/E[0]; S_2D.data()[1] = Sxy;
    S_2D.data()[3] = Sxy; S_2D.data()[4] = 1/E[1];
    S_2D.data()[8] = 1/G[0];
    this->D_2D = S_2D.get_inverted_cholesky();

    S_3D.data()[ 0] = 1/E[0]; S_3D.data()[ 1] = Sxy; S_3D.data()[ 2] = Sxz;
    S_3D.data()[ 6] = Sxy; S_3D.data()[ 7] = 1/E[1]; S_3D.data()[ 8] = Syz;
    S_3D.data()[12] = Sxz; S_3D.data()[13] = Syz; S_3D.data()[14] = 1/E[2];
    S_3D.data()[21] = 1/G[0];
    S_3D.data()[28] = 1/G[1];
    S_3D.data()[35] = 1/G[2];
    this->D_3D = S_3D.get_inverted_cholesky();
}

LinearElasticOrthotropicField::LinearElasticOrthotropicField(const projspec::DataMap& data):
    Material(data.get_string("name"),
            {data.get_array("Smax")->get_double_array()},
            {data.get_array("Tmax")->get_double_array()}),
    density(data.get_double("density")),
    D_2D(), D_3D(), S_2D(3,3), S_3D(6,6){

    if(!data.get_bool("STUB", false)){
        Field* f = data.proj->fields[data.get_int("field")].get();
        logger::log_assert(f->get_type() == Field::Type::COORDINATE, logger::ERROR, "orthotropic_field material type requires a coordinate field");
        logger::log_assert(f->get_sub_type() == Field::SubType::DOMAIN || f->get_sub_type() == Field::SubType::PROJECTION, logger::ERROR, "field subtype for orthotropic field material must be DOMAIN or PROJECTION");

        this->field = static_cast<CoordinateField*>(f);

        const auto nu = data.get_array("nu")->get_double_array();
        const auto nu_lower_half = data.get_array("nu_lower_half")->get_bool_array();
        auto E = data.get_array("E")->get_double_array();
        auto G = data.get_array("G")->get_double_array();
        for(size_t i = 0; i < 3; ++i){
            E[i] *= 1e3;
            G[i] *= 1e3;
        }
       
        double Sxy = 0, Sxz = 0, Syz = 0;
        if(nu_lower_half[0]){
            Sxy = -nu[0]/E[0];
        } else {
            Sxy = -nu[0]/E[1];
        }
        if(nu_lower_half[1]){
            Sxz = -nu[1]/E[0];
        } else {
            Sxz = -nu[1]/E[2];
        }
        if(nu_lower_half[2]){
            Syz = -nu[2]/E[1];
        } else {
            Syz = -nu[2]/E[2];
        }

        S_2D.data()[0] = 1/E[0]; S_2D.data()[1] = Sxy;
        S_2D.data()[3] = Sxy; S_2D.data()[4] = 1/E[1];
        S_2D.data()[8] = 1/G[0];
        this->D_2D = S_2D.get_inverted_cholesky();

        S_3D.data()[ 0] = 1/E[0]; S_3D.data()[ 1] = Sxy; S_3D.data()[ 2] = Sxz;
        S_3D.data()[ 6] = Sxy; S_3D.data()[ 7] = 1/E[1]; S_3D.data()[ 8] = Syz;
        S_3D.data()[12] = Sxz; S_3D.data()[13] = Syz; S_3D.data()[14] = 1/E[2];
        S_3D.data()[21] = 1/G[0];
        S_3D.data()[28] = 1/G[1];
        S_3D.data()[35] = 1/G[2];
        this->D_3D = S_3D.get_inverted_cholesky();
    }
}
math::Matrix LinearElasticOrthotropicField::stiffness_2D(const MeshElement* const e, const gp_Pnt& p) const{
    const auto M = this->field->get_matrix(e, p);
    const auto R = utils::basis_tensor_2D(M);
    const auto& D = this->D_2D;

    return R*D*R.T();
}
math::Matrix LinearElasticOrthotropicField::stiffness_3D(const MeshElement* const e, const gp_Pnt& p) const{
    const auto M = this->field->get_matrix(e, p);
    const auto R = utils::basis_tensor_3D(M);
    const auto& D = this->D_3D;

    return R*D*R.T();
}
math::Matrix LinearElasticOrthotropicField::stiffness_inverse_2D(const MeshElement* const e, const gp_Pnt& p) const{
    const auto M = this->field->get_matrix(e, p);
    const auto R = utils::basis_tensor_2D_inv_T(M);
    const auto& S = this->S_2D;

    return R*S*R.T();
}
math::Matrix LinearElasticOrthotropicField::stiffness_inverse_3D(const MeshElement* const e, const gp_Pnt& p) const{
    const auto M = this->field->get_matrix(e, p);
    const auto R = utils::basis_tensor_3D_inv_T(M);
    const auto& S = this->S_3D;

    return R*S*R.T();
}

double LinearElasticOrthotropicField::beam_E_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const{
    // Anisotropic Elasticity: Theory and Applications
    // (Ting, 1996)

    const auto S = this->stiffness_inverse_2D(e, p);
    const auto Rt = utils::basis_tensor_2D_inv_T(R);

    return 1.0/(Rt*S*Rt.T())(0,0);
}
double LinearElasticOrthotropicField::beam_E_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const{
    // Anisotropic Elasticity: Theory and Applications
    // (Ting, 1996)

    auto S = this->stiffness_inverse_3D(e, p);
    const auto Rt = utils::basis_tensor_3D_inv_T(R);

    return 1.0/(Rt*S*Rt.T())(0,0);
}
std::array<double, 2> LinearElasticOrthotropicField::beam_EG_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const{
    // Anisotropic Elasticity: Theory and Applications
    // (Ting, 1996)

    auto S = this->stiffness_inverse_2D(e, p);
    const auto Rt = utils::basis_tensor_2D_inv_T(R);
    const auto Sr = Rt*S*Rt.T();

    return {1.0/Sr(0,0), 1.0/S(2,2)};
}
std::array<double, 4> LinearElasticOrthotropicField::beam_EG_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const{
    // Anisotropic Elasticity: Theory and Applications
    // (Ting, 1996)

    auto S = this->stiffness_inverse_3D(e, p);
    const auto Rt = utils::basis_tensor_3D_inv_T(R);
    const auto Sr = Rt*S*Rt.T();

    return {1.0/S(0,0), 1.0/S(3,3), 1.0/S(4,4), 1.0/S(5,5)};
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


using namespace projspec;
const bool LinearElasticOrthotropicField::reg = Factory<Material>::add(
    [](const DataMap& data){
        return std::make_unique<LinearElasticOrthotropicField>(data);
    },
    ObjectRequirements{
        "linear_elastic_orthotropic_field",
        {
            DataEntry{.name = "name", .type = TYPE_STRING, .required = true},
            DataEntry{.name = "field", .type = TYPE_INT, .required = true},
            DataEntry{.name = "Smax", .type = TYPE_ARRAY,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 3,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                         .required_if = 
                             [](const RequirementConditions& r){
                                 return r.generate_beams;
                             }
                     },
            DataEntry{.name = "Tmax", .type = TYPE_ARRAY,
                          .array_data = std::shared_ptr<ArrayRequirements>(
                              new ArrayRequirements{
                                  .size = 3,
                                  .type = TYPE_DOUBLE
                              }
                          ),
                          .required_if = 
                              [](const RequirementConditions& r){
                                  return r.generate_beams;
                              }
                      },
            DataEntry{.name = "E", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 3,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
            DataEntry{.name = "G", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 3,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
            DataEntry{.name = "nu", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 3,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
            DataEntry{.name = "nu_lower_half", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 3,
                                 .type = TYPE_BOOL
                             }
                         ),
                     },
            DataEntry{.name = "density", .type = TYPE_DOUBLE, .required = true}
        }
    }
);

}
