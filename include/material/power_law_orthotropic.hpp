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

#ifndef POWER_LAW_ORTHOTROPIC_HPP
#define POWER_LAW_ORTHOTROPIC_HPP

#include "material.hpp"
#include "field.hpp"
#include "project_specification/data_map.hpp"

namespace material{

class PowerLawOrthotropic : public Material{
    public:
    PowerLawOrthotropic(const projspec::DataMap& data);

    virtual math::Matrix stiffness_2D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual math::Matrix stiffness_3D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual math::Matrix stiffness_inverse_2D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual math::Matrix stiffness_inverse_3D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual double get_density(const MeshElement* const e, const gp_Pnt& p) const override{
        return this->density_field->get(e, p);
    }

    virtual math::Matrix stiffness_inverse_2D_dsh(const MeshElement* const e, const gp_Pnt& p, const size_t n, const size_t dof) const override;
    virtual math::Matrix stiffness_inverse_3D_dsh(const MeshElement* const e, const gp_Pnt& p, const size_t n, const size_t dof) const override;
    virtual math::Matrix stiffness_2D_dsh(const MeshElement* const e, const gp_Pnt& p, const size_t n, const size_t dof) const override;
    virtual math::Matrix stiffness_3D_dsh(const MeshElement* const e, const gp_Pnt& p, const size_t n, const size_t dof) const override;
    virtual std::vector<math::Matrix> stiffness_2D_dsh_full(const MeshElement* const e, const gp_Pnt& p, const size_t n, const size_t dof) const override;
    virtual std::vector<math::Matrix> stiffness_3D_dsh_full(const MeshElement* const e, const gp_Pnt& p, const size_t n, const size_t dof) const override;

    virtual double beam_E_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override;
    virtual double beam_E_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override;
    virtual std::array<double, 2> beam_EG_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override;
    virtual std::array<double, 4> beam_EG_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override;

    virtual math::Matrix stiffness_2D_orig(const MeshElement* const e, const gp_Pnt& p) const override{
        auto D = this->S_2D_base(e, p);
        D.invert_cholesky();

        return D;
    }
    virtual math::Matrix stiffness_3D_orig(const MeshElement* const e, const gp_Pnt& p) const override{
        auto D = this->S_3D_base(e, p);
        D.invert_cholesky();

        return D;
    }
    virtual math::Matrix stiffness_inverse_2D_orig(const MeshElement* const e, const gp_Pnt& p) const override{
        return this->S_2D_base(e, p);
    }
    virtual math::Matrix stiffness_inverse_3D_orig(const MeshElement* const e, const gp_Pnt& p) const override{
        return this->S_3D_base(e, p);
    }

    virtual math::Matrix stiffness_2D_orig_dsh(const MeshElement* const e, const gp_Pnt& p, const size_t n, const size_t dof) const override{
        const auto S_2D = this->S_2D_base(e, p);
        const auto dS_2D = this->S_2D_base_dsh(e, p, n, dof);
        const auto D_2D = S_2D.get_inverted_cholesky();
        const auto dD_2D = -D_2D*dS_2D*D_2D;

        return dD_2D;
    }
    virtual math::Matrix stiffness_3D_orig_dsh(const MeshElement* const e, const gp_Pnt& p, const size_t n, const size_t dof) const override{
        const auto S_3D = this->S_3D_base(e, p);
        const auto dS_3D = this->S_3D_base_dsh(e, p, n, dof);
        const auto D_3D = S_3D.get_inverted_cholesky();
        const auto dD_3D = -D_3D*dS_3D*D_3D;

        return dD_3D;
    }
    virtual math::Matrix stiffness_inverse_2D_orig_dsh(const MeshElement* const e, const gp_Pnt& p, const size_t n, const size_t dof) const override{
        return this->S_2D_base_dsh(e, p, n, dof);
    }
    virtual math::Matrix stiffness_inverse_3D_orig_dsh(const MeshElement* const e, const gp_Pnt& p, const size_t n, const size_t dof) const override{
        return this->S_3D_base_dsh(e, p, n, dof);
    }

    virtual Class get_class() const override{ return this->ORTHOTROPIC; }
    virtual Type get_type() const override{ return this->POWER_LAW_ORTHOTROPIC; }
    virtual bool is_homogeneous() const override{ return false; }
    virtual Field* get_field() const override{
        if(this->direction_field != nullptr){
            return this->direction_field;
        } else {
            return this->density_field;
        }
    }

    virtual std::vector<double> get_max_stresses(gp_Dir d) const override{
       (void) d; 
       return {};
    }

    private:
    math::Matrix S_2D_base(const MeshElement* const e, const gp_Pnt& p) const;
    math::Matrix S_3D_base(const MeshElement* const e, const gp_Pnt& p) const;
    math::Matrix S_2D_base_dsh(const MeshElement* const e, const gp_Pnt& p, const size_t n, const size_t dof) const;
    math::Matrix S_3D_base_dsh(const MeshElement* const e, const gp_Pnt& p, const size_t n, const size_t dof) const;
    std::vector<math::Matrix> S_2D_base_dsh_full(const MeshElement* const e, const gp_Pnt& p, const size_t n, const size_t dof) const;
    std::vector<math::Matrix> S_3D_base_dsh_full(const MeshElement* const e, const gp_Pnt& p, const size_t n, const size_t dof) const;
    static const bool reg;
    ScalarField* density_field;
    CoordinateField* direction_field;
    const bool plane_stress;
    std::array<double, 2> E1;
    std::array<double, 2> E2;
    std::array<double, 2> E3;
    std::array<double, 2> G12;
    std::array<double, 2> G13;
    std::array<double, 2> G23;
    std::array<double, 3> nu;

};

}

#endif
