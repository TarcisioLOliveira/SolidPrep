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

#ifndef ORTHOTROPIC_ELASTIC_ISOTROPIC_HPP
#define ORTHOTROPIC_ELASTIC_ISOTROPIC_HPP

#include "material.hpp"

namespace material{

class LinearElasticOrthotropic : public Material{
    public:
    LinearElasticOrthotropic(std::vector<double> E, std::vector<double> nu, std::vector<double> G, std::vector<double> Smax, std::vector<double> Tmax);

    virtual std::vector<double> stiffness_2D() const override;
    virtual std::vector<double> stiffness_3D() const override;

    virtual double beam_E_2D(gp_Dir d) const override;
    virtual double beam_E_3D(gp_Dir d) const override;

    virtual Type get_type() const override{ return this->LINEAR_ELASTIC_ORTHOTROPIC; }

    virtual std::vector<double> get_max_stresses(gp_Dir d) const override;

    private:
    std::vector<double> D_2D;
    std::vector<double> D_3D;
};

}

#endif
