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
#include "material/power_law_isotropic.hpp"

namespace material{

PowerLawIsotropic::PowerLawIsotropic(const projspec::DataMap& data):
    Material(data.get_string("name"),
            {0, 0, 0},
            {0, 0, 0}),
    plane_stress(data.get_bool("plane_stress")),
    E(data.get_array("E")->get_double_array_fixed<2>()),
    nu(data.get_double("nu")){

    if(!data.get_bool("STUB", false)){
        E[0] *= 1e3;
        Field* f = data.proj->fields[data.get_int("density_field")].get();
        logger::log_assert(f->get_type() == Field::Type::SCALAR, logger::ERROR, "bone material type requires a density field");
        logger::log_assert(f->get_sub_type() == Field::SubType::DOMAIN, logger::ERROR, "field subtype for density field material must be DOMAIN");

        this->density_field = static_cast<ScalarField*>(f);
    }
}

math::Matrix PowerLawIsotropic::stiffness_2D(const MeshElement* const e, const gp_Pnt& p) const{
    const double rho = this->density_field->get(e, p);
    const double Ee = E[0]*std::pow(rho, E[1]);
    math::Matrix D_2D(3,3);
    if(plane_stress){
        const double Eii = Ee/(1-nu*nu);
        const double Eij = (Ee/(1-nu*nu))*nu;
        D_2D.data()[0] = Eii;
        D_2D.data()[1] = Eij;
        D_2D.data()[3] = Eij;
        D_2D.data()[4] = Eii;
        D_2D.data()[8] = (Ee/(1-nu*nu))*(1-nu)/2;
    } else {
        const double Eii = (Ee/((1+nu)*(1-2*nu)))*(1-nu);
        const double Eij = (Ee/((1+nu)*(1-2*nu)))*(nu);
        D_2D.data()[0] = Eii;
        D_2D.data()[1] = Eij;
        D_2D.data()[3] = Eij;
        D_2D.data()[4] = Eii;
        D_2D.data()[8] = Ee/(2*(1+nu));
    }

    return D_2D;
}
math::Matrix PowerLawIsotropic::stiffness_3D(const MeshElement* const e, const gp_Pnt& p) const{
    const double rho = this->density_field->get(e, p);
    const double Ee = E[0]*std::pow(rho, E[1]);
    math::Matrix D_3D(6,6);
    const double Eii = (Ee/((1+nu)*(1-2*nu)))*(1-nu);
    const double Eij = (Ee/((1+nu)*(1-2*nu)))*(nu);
    const double G   =  Ee/(2*(1+nu));
    D_3D.data()[ 0] = Eii;
    D_3D.data()[ 1] = Eij;
    D_3D.data()[ 2] = Eij;
    D_3D.data()[ 6] = Eij;
    D_3D.data()[ 7] = Eii;
    D_3D.data()[ 8] = Eij;
    D_3D.data()[12] = Eij;
    D_3D.data()[13] = Eij;
    D_3D.data()[14] = Eii;
    D_3D.data()[21] = G;
    D_3D.data()[28] = G;
    D_3D.data()[35] = G;

    return D_3D;
}
math::Matrix PowerLawIsotropic::stiffness_inverse_2D(const MeshElement* const e, const gp_Pnt& p) const{
    auto S_2D = this->stiffness_2D(e, p);
    S_2D.invert_cholesky();
    return S_2D;
}
math::Matrix PowerLawIsotropic::stiffness_inverse_3D(const MeshElement* const e, const gp_Pnt& p) const{
    auto S_3D = this->stiffness_3D(e, p);
    S_3D.invert_cholesky();
    return S_3D;
}

double PowerLawIsotropic::beam_E_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const{
    (void) R;
    const double rho = this->density_field->get(e, p);
    const double Ee = E[0]*std::pow(rho, E[1]);

    return Ee;
}
double PowerLawIsotropic::beam_E_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const{
    (void) R;
    const double rho = this->density_field->get(e, p);
    const double Ee = E[0]*std::pow(rho, E[1]);

    return Ee;
}
std::array<double, 2> PowerLawIsotropic::beam_EG_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const{
    (void) R;
    const double rho = this->density_field->get(e, p);
    const double Ee = E[0]*std::pow(rho, E[1]);
    double Ge = Ee/(2*(1 + nu));

    return {Ee, Ge};
}
std::array<double, 4> PowerLawIsotropic::beam_EG_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const{
    (void) R;
    const double rho = this->density_field->get(e, p);
    const double Ee = E[0]*std::pow(rho, E[1]);
    double Ge = Ee/(2*(1 + nu));

    return {Ee, Ge, Ge, Ge};
}

using namespace projspec;
const bool PowerLawIsotropic::reg = Factory<Material>::add(
    [](const DataMap& data){
        return std::make_unique<PowerLawIsotropic>(data);
    },
    ObjectRequirements{
        "power_law_isotropic",
        {
            DataEntry{.name = "name", .type = TYPE_STRING, .required = true},
            DataEntry{.name = "density_field", .type = TYPE_INT, .required = true},
            DataEntry{.name = "plane_stress", .type = TYPE_BOOL, .required = true},
            DataEntry{.name = "E", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 2,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
            DataEntry{.name = "nu", .type = TYPE_DOUBLE, .required = true},
        }
    }
);

}
