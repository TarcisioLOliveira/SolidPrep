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

#ifndef BONE_HPP
#define BONE_HPP

#include "material.hpp"
#include "field.hpp"
#include "project_specification/data_map.hpp"

namespace material{

class Bone : public Material{
    public:
    Bone(const projspec::DataMap& data);

    virtual math::Matrix stiffness_2D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual math::Matrix stiffness_3D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual math::Matrix stiffness_inverse_2D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual math::Matrix stiffness_inverse_3D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual double get_density(const MeshElement* const e, const gp_Pnt& p) const override;

    virtual double beam_E_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override;
    virtual double beam_E_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override;
    virtual std::array<double, 2> beam_EG_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override;
    virtual std::array<double, 4> beam_EG_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override;

    virtual Class get_class() const override{
        // Temporary
        if(outer.get() != nullptr){
            return std::max(inner->get_class(), outer->get_class());
        } else {
            return Material::Class::ISOTROPIC;
        }
    }
    virtual Type get_type() const override{ return this->BONE; }
    virtual bool is_homogeneous() const override{ return false; }
    virtual Field* get_field() const override{
        return this->outer->get_field();
    }

    virtual std::vector<double> get_max_stresses(gp_Dir d) const override{
        (void)d;
        return {};
    }

    private:

    static const bool reg;
    ScalarField* trabecular_field;
    const utils::DelayedPointerView<Material> outer;
    const utils::DelayedPointerView<Material> inner;
};

}


#endif
