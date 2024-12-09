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

#ifndef LINEAR_ELASTIC_ISOTROPIC_HPP
#define LINEAR_ELASTIC_ISOTROPIC_HPP

#include "material.hpp"
#include "math/matrix.hpp"

namespace material{

class LinearElasticIsotropic : public Material{
    public:
    LinearElasticIsotropic(const std::string& name, const double density, double E, double nu, double Smax, double Tmax, bool plane_stress);

    inline virtual math::Matrix stiffness_2D(const MeshElement* const e, const gp_Pnt& p) const override{
        (void)p;
        (void)e;
        return this->D_2D;
    }
    inline virtual math::Matrix stiffness_3D(const MeshElement* const e, const gp_Pnt& p) const override{
        (void)p;
        (void)e;
        return this->D_3D;
    }
    inline virtual math::Matrix stiffness_inverse_2D(const MeshElement* const e, const gp_Pnt& p) const override{
        (void)p;
        (void)e;
        return this->S_2D;
    }
    inline virtual math::Matrix stiffness_inverse_3D(const MeshElement* const e, const gp_Pnt& p) const override{
        (void)p;
        (void)e;
        return this->S_3D;
    }
    inline virtual double get_density(const MeshElement* const e, const gp_Pnt& p) const override{
        (void)p;
        (void)e;
        return density;
    }

    inline virtual double beam_E_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override{
        (void)R;
        (void)p;
        (void)e;
        return this->E;
    }

    virtual double beam_E_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override{
        (void)R;
        (void)p;
        (void)e;
        return this->E;
    }
    virtual std::array<double, 2> beam_EG_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override{
        (void)R;
        (void)p;
        (void)e;
        return {E, G};
    }
    virtual std::array<double, 4> beam_EG_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override{
        (void)R;
        (void)p;
        (void)e;
        return {E, G, G, G};
    }

    virtual Type get_type() const override{ return this->LINEAR_ELASTIC_ISOTROPIC; }

    virtual std::vector<double> get_max_stresses(gp_Dir d) const override;

    private:
    const double E;
    const double G;
    const double nu;
    const double density;
    math::Matrix D_2D;
    math::Matrix D_3D;
    math::Matrix S_2D;
    math::Matrix S_3D;
};

}

#endif
