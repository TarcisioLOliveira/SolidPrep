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

#ifndef POWER_LAW_ISOTROPIC_HPP
#define POWER_LAW_ISOTROPIC_HPP

#include "material.hpp"
#include "field.hpp"
#include "project_specification/data_map.hpp"

namespace material{

class PowerLawIsotropic : public Material{
    public:
    PowerLawIsotropic(const projspec::DataMap& data);

    virtual math::Matrix stiffness_2D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual math::Matrix stiffness_3D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual math::Matrix stiffness_inverse_2D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual math::Matrix stiffness_inverse_3D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual double get_density(const MeshElement* const e, const gp_Pnt& p) const override{
        return this->density_field->get(e, p);
    }

    virtual double beam_E_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override;
    virtual double beam_E_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override;
    virtual std::array<double, 2> beam_EG_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override;
    virtual std::array<double, 4> beam_EG_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override;

    virtual Class get_class() const override{ return this->ISOTROPIC; }
    virtual Type get_type() const override{ return this->POWER_LAW_ISOTROPIC; }
    virtual bool is_homogeneous() const override{ return false; }
    virtual Field* get_field() const override{
        return this->density_field;
    }

    virtual std::vector<double> get_max_stresses(gp_Dir d) const override{
       (void) d; 
       return {};
    }

    private:
    static const bool reg;
    ScalarField* density_field;
    const bool plane_stress;
    std::array<double, 2> E;
    const double nu;

};

}

#endif
