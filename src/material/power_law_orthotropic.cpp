/*
 *   Copyright (C) 2025 Tarcísio Ladeia de Oliveira.
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

#include "logger.hpp"
#include "project_data.hpp"
#include "utils/basis_tensor.hpp"
#include "material/power_law_orthotropic.hpp"

namespace material{

PowerLawOrthotropic::PowerLawOrthotropic(const projspec::DataMap& data):
    Material(data.get_string("name"),
            {0, 0, 0},
            {0, 0, 0}),
    direction_field(nullptr),
    plane_stress(data.get_bool("plane_stress")),
    E1(data.get_array("E1")->get_double_array_fixed<2>()),
    E2(data.get_array("E2")->get_double_array_fixed<2>()),
    E3(data.get_array("E3")->get_double_array_fixed<2>()),
    G12(data.get_array("G12")->get_double_array_fixed<2>()),
    G13(data.get_array("G13")->get_double_array_fixed<2>()),
    G23(data.get_array("G23")->get_double_array_fixed<2>()),
    nu(data.get_array("nu")->get_double_array_fixed<3>()){

    if(!data.get_bool("STUB", false)){
        E1[0] *= 1e3;
        E2[0] *= 1e3;
        E3[0] *= 1e3;
        G12[0] *= 1e3;
        G13[0] *= 1e3;
        G23[0] *= 1e3;

        Field* f = data.proj->fields[data.get_int("density_field")].get();
        logger::log_assert(f->get_type() == Field::Type::SCALAR, logger::ERROR, "bone material type requires a density field");
        logger::log_assert(f->get_sub_type() == Field::SubType::DOMAIN, logger::ERROR, "field subtype for density field material must be DOMAIN");

        this->density_field = static_cast<ScalarField*>(f);

        if(data.exists_int("direction_field")){
            f = data.proj->fields[data.get_int("direction_field")].get();
            logger::log_assert(f->get_type() == Field::Type::COORDINATE, logger::ERROR, "direction field must be a coordinate field");
            logger::log_assert(f->get_sub_type() == Field::SubType::DOMAIN || f->get_sub_type() == Field::SubType::PROJECTION, logger::ERROR, "field subtype for direction field material must be DOMAIN or PROJECTION");

            this->direction_field = static_cast<CoordinateField*>(f);
        }

    }
}
math::Matrix PowerLawOrthotropic::S_2D_base(const MeshElement* const e, const gp_Pnt& p) const{
    math::Matrix S_2D(3,3);
    const double rho = this->density_field->get(e, p);
    const std::array<double, 2> E = {
        E1[0]*std::pow(rho, E1[1]), 
        E2[0]*std::pow(rho, E2[1])
    };

    const double G = G12[0]*std::pow(rho, G12[1]);
    
    const double Sxy = -nu[0]/E[0];

    S_2D.data()[0] = 1/E[0]; S_2D.data()[1] = Sxy;
    S_2D.data()[3] = Sxy; S_2D.data()[4] = 1/E[1];
    S_2D.data()[8] = 1/G;

    return S_2D;
}
math::Matrix PowerLawOrthotropic::S_3D_base(const MeshElement* const e, const gp_Pnt& p) const{
    math::Matrix S_3D(6,6);
    const double rho = this->density_field->get(e, p);
    const std::array<double, 3> E = {
        E1[0]*std::pow(rho, E1[1]), 
        E2[0]*std::pow(rho, E2[1]),
        E3[0]*std::pow(rho, E3[1])
    };
    const std::array<double, 3> G = {
        G12[0]*std::pow(rho, G12[1]), 
        G13[0]*std::pow(rho, G13[1]),
        G23[0]*std::pow(rho, G23[1])
    };
    const double Sxy = -nu[0]/E[0];
    const double Sxz = -nu[1]/E[0];
    const double Syz = -nu[2]/E[1];

    S_3D.data()[ 0] = 1/E[0]; S_3D.data()[ 1] = Sxy; S_3D.data()[ 2] = Sxz;
    S_3D.data()[ 6] = Sxy; S_3D.data()[ 7] = 1/E[1]; S_3D.data()[ 8] = Syz;
    S_3D.data()[12] = Sxz; S_3D.data()[13] = Syz; S_3D.data()[14] = 1/E[2];
    S_3D.data()[21] = 1/G[0];
    S_3D.data()[28] = 1/G[1];
    S_3D.data()[35] = 1/G[2];

    return S_3D;
}
math::Matrix PowerLawOrthotropic::stiffness_2D(const MeshElement* const e, const gp_Pnt& p) const{
    auto D_2D = this->S_2D_base(e, p);
    D_2D.invert_cholesky();

    if(this->direction_field != nullptr){
        const auto M = this->direction_field->get_matrix(e, p);
        const auto R = utils::basis_tensor_2D(M.T());
        
        return R*D_2D*R.T();
    } else {
        return D_2D;
    }
}
math::Matrix PowerLawOrthotropic::stiffness_3D(const MeshElement* const e, const gp_Pnt& p) const{
    auto D_3D = this->S_3D_base(e, p);
    D_3D.invert_cholesky();

    if(this->direction_field != nullptr){
        const auto M = this->direction_field->get_matrix(e, p);
        const auto R = utils::basis_tensor_3D(M.T());
        
        return R*D_3D*R.T();
    } else {
        return D_3D;
    }
}
math::Matrix PowerLawOrthotropic::stiffness_inverse_2D(const MeshElement* const e, const gp_Pnt& p) const{
    const auto S_2D = this->S_2D_base(e, p);

    if(this->direction_field != nullptr){
        const auto M = this->direction_field->get_matrix(e, p);
        const auto R = utils::basis_tensor_2D_inv_T(M.T());
        
        return R*S_2D*R.T();
    } else {
        return S_2D;
    }
}
math::Matrix PowerLawOrthotropic::stiffness_inverse_3D(const MeshElement* const e, const gp_Pnt& p) const{
    const auto S_3D = this->S_3D_base(e, p);

    if(this->direction_field != nullptr){
        const auto M = this->direction_field->get_matrix(e, p);
        const auto R = utils::basis_tensor_3D_inv_T(M.T());
        
        return R*S_3D*R.T();
    } else {
        return S_3D;
    }
}

double PowerLawOrthotropic::beam_E_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const{
    // Anisotropic Elasticity: Theory and Applications
    // (Ting, 1996)

    const auto S = this->stiffness_inverse_2D(e, p);
    const auto Rt = utils::basis_tensor_2D_inv_T(R);

    return 1.0/(Rt*S*Rt.T())(0,0);
}
double PowerLawOrthotropic::beam_E_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const{
    // Anisotropic Elasticity: Theory and Applications
    // (Ting, 1996)

    auto S = this->stiffness_inverse_3D(e, p);
    const auto Rt = utils::basis_tensor_3D_inv_T(R);

    return 1.0/(Rt*S*Rt.T())(0,0);
}
std::array<double, 2> PowerLawOrthotropic::beam_EG_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const{
    // Anisotropic Elasticity: Theory and Applications
    // (Ting, 1996)

    auto S = this->stiffness_inverse_2D(e, p);
    const auto Rt = utils::basis_tensor_2D_inv_T(R);
    const auto Sr = Rt*S*Rt.T();

    return {1.0/Sr(0,0), 1.0/Sr(2,2)};
}
std::array<double, 4> PowerLawOrthotropic::beam_EG_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const{
    // Anisotropic Elasticity: Theory and Applications
    // (Ting, 1996)

    auto S = this->stiffness_inverse_3D(e, p);
    const auto Rt = utils::basis_tensor_3D_inv_T(R);
    const auto Sr = Rt*S*Rt.T();

    return {1.0/Sr(0,0), 1.0/Sr(3,3), 1.0/Sr(4,4), 1.0/Sr(5,5)};
}

using namespace projspec;
const bool PowerLawOrthotropic::reg = Factory<Material>::add(
    [](const DataMap& data){
        return std::make_unique<PowerLawOrthotropic>(data);
    },
    ObjectRequirements{
        "power_law_orthotropic",
        {
            DataEntry{.name = "name", .type = TYPE_STRING, .required = true},
            DataEntry{.name = "density_field", .type = TYPE_INT, .required = true},
            DataEntry{.name = "direction_field", .type = TYPE_INT, .required = false},
            DataEntry{.name = "plane_stress", .type = TYPE_BOOL, .required = true},
            DataEntry{.name = "E1", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 2,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
            DataEntry{.name = "E2", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 2,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
            DataEntry{.name = "E3", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 2,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
            DataEntry{.name = "G12", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 2,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
            DataEntry{.name = "G13", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 2,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
            DataEntry{.name = "G23", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 2,
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
        }
    }
);
}
