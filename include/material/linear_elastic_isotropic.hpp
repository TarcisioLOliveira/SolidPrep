/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   Foobar is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Foobar is distributed in the hope that it will be useful,
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

namespace material{

class LinearElasticIsotropic : public Material{
    public:
    LinearElasticIsotropic(float E, float nu, float Smax, float Tmax, bool plane_stress);

    virtual std::vector<float> stiffness_2D() const override;
    virtual std::vector<float> stiffness_3D() const override;

    virtual float beam_E_2D(gp_Dir d) const override;
    virtual float beam_E_3D(gp_Dir d) const override;

    virtual Type get_type() const override{ return this->LINEAR_ELASTIC_ISOTROPIC; }

    virtual std::vector<float> get_max_stresses(gp_Dir d) const override;

    private:
    const float E;
    const float nu;

    std::vector<float> D_2D;
    std::vector<float> D_3D;
};

}

#endif
