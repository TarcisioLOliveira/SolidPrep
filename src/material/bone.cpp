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

#include "material/bone.hpp"
#include "logger.hpp"
#include "project_data.hpp"

namespace material{

Bone::Bone(const projspec::DataMap& data):
    Material(data.get_string("name"),
            {0, 0, 0},
            {0, 0, 0}),
    outer(data.proj->get_material(data.get_string("material_outer"))),
    inner(data.proj->get_material(data.get_string("material_inner")))
{

    if(!data.get_bool("STUB", false)){
        Field* f = data.proj->fields[data.get_int("trabecular_map")].get();
        logger::log_assert(f->get_type() == Field::Type::SCALAR, logger::ERROR, "bone material type requires a scalar field");
        logger::log_assert(f->get_sub_type() == Field::SubType::DOMAIN, logger::ERROR, "field subtype for density field material must be DOMAIN");

        this->trabecular_field = static_cast<ScalarField*>(f);

    }
}
math::Matrix Bone::stiffness_2D(const MeshElement* const e, const gp_Pnt& p) const{
    const auto t = this->trabecular_field->get(e, p);
    if(t == 1){
        return this->inner->stiffness_2D(e, p);
    } else {
        return this->outer->stiffness_2D(e, p);
    }
}
math::Matrix Bone::stiffness_3D(const MeshElement* const e, const gp_Pnt& p) const{
    const auto t = this->trabecular_field->get(e, p);
    if(t == 1){
        return this->inner->stiffness_3D(e, p);
    } else {
        return this->outer->stiffness_3D(e, p);
    }
}
math::Matrix Bone::stiffness_inverse_2D(const MeshElement* const e, const gp_Pnt& p) const{
    const auto t = this->trabecular_field->get(e, p);
    if(t == 1){
        return this->inner->stiffness_inverse_2D(e, p);
    } else {
        return this->outer->stiffness_inverse_2D(e, p);
    }
}
math::Matrix Bone::stiffness_inverse_3D(const MeshElement* const e, const gp_Pnt& p) const{
    const auto t = this->trabecular_field->get(e, p);
    if(t == 1){
        return this->inner->stiffness_inverse_3D(e, p);
    } else {
        return this->outer->stiffness_inverse_3D(e, p);
    }
}
double Bone::get_density(const MeshElement* const e, const gp_Pnt& p) const{
    const auto t = this->trabecular_field->get(e, p);
    if(t == 1){
        return this->inner->get_density(e, p);
    } else {
        return this->outer->get_density(e, p);
    }
}

double Bone::beam_E_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const{
    const auto t = this->trabecular_field->get(e, p);
    if(t == 1){
        return this->inner->beam_E_2D(e, p, R);
    } else {
        return this->outer->beam_E_2D(e, p, R);
    }
}
double Bone::beam_E_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const{
    const auto t = this->trabecular_field->get(e, p);
    if(t == 1){
        return this->inner->beam_E_3D(e, p, R);
    } else {
        return this->outer->beam_E_3D(e, p, R);
    }
}
std::array<double, 2> Bone::beam_EG_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const{
    const auto t = this->trabecular_field->get(e, p);
    if(t == 1){
        return this->inner->beam_EG_2D(e, p, R);
    } else {
        return this->outer->beam_EG_2D(e, p, R);
    }
}
std::array<double, 4> Bone::beam_EG_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const{
    const auto t = this->trabecular_field->get(e, p);
    if(t == 1){
        return this->inner->beam_EG_3D(e, p, R);
    } else {
        return this->outer->beam_EG_3D(e, p, R);
    }
}

using namespace projspec;
const bool Bone::reg = Factory<Material>::add(
    [](const DataMap& data){
        return std::make_unique<Bone>(data);
    },
    ObjectRequirements{
        "bone",
        {
            DataEntry{.name = "name", .type = TYPE_STRING, .required = true},
            DataEntry{.name = "trabecular_map", .type = TYPE_INT, .required = true},
            DataEntry{.name = "material_outer", .type = TYPE_STRING, .required = true},
            DataEntry{.name = "material_inner", .type = TYPE_STRING, .required = true},
        }
    }
);

}
