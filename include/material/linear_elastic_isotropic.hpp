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

namespace material{

class LinearElasticIsotropic : public Material{
    public:
    LinearElasticIsotropic(const std::string& name, const double density, double E, double nu, double Smax, double Tmax, bool plane_stress);

    virtual double beam_E_2D(gp_Dir d) const override;
    virtual double beam_E_3D(gp_Dir d) const override;

    virtual Type get_type() const override{ return this->LINEAR_ELASTIC_ISOTROPIC; }

    virtual std::vector<double> get_max_stresses(gp_Dir d) const override;

    private:
    const double E;
    const double nu;
};

}

#endif
