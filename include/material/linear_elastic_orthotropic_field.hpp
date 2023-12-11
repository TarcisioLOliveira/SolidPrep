/*
 *   Copyright (C) 2023 Tarcísio Ladeia de Oliveira.
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

#ifndef LINEAR_ELASTIC_ORTHOTROPIC_FIELD_HPP
#define LINEAR_ELASTIC_ORTHOTROPIC_FIELD_HPP

#include "material.hpp"
#include "field.hpp"

namespace material{

class LinearElasticOrthotropicField : public Material{
    public:
    LinearElasticOrthotropicField(const std::string& name, const double density, std::vector<double> E, std::vector<double> nu, std::vector<double> G, std::vector<double> Smax, std::vector<double> Tmax, const CoordinateField* field);

    virtual std::vector<double> stiffness_2D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual std::vector<double> stiffness_3D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual std::vector<double> stiffness_inverse_2D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual std::vector<double> stiffness_inverse_3D(const MeshElement* const e, const gp_Pnt& p) const override;
    inline virtual double get_density(const MeshElement* const e, const gp_Pnt& p) const override{
        (void)e;
        (void)p;
        return this->density;
    }

    virtual double beam_E_2D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 2, 2>& R) const override;
    virtual double beam_E_3D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 3, 3>& R) const override;
    virtual std::array<double, 2> beam_EG_2D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 2, 2>& R) const override;
    virtual std::array<double, 4> beam_EG_3D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 3, 3>& R) const override;
    virtual double S12_2D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 2, 2>& R) const override;
    virtual std::array<double, 2> S12_S13_3D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 3, 3>& R) const override;

    virtual Type get_type() const override{ return this->LINEAR_ELASTIC_ORTHOTROPIC_FIELD; }

    virtual std::vector<double> get_max_stresses(gp_Dir d) const override;

    private:
    const double density;
    const CoordinateField* field;
    std::vector<double> D_2D;
    std::vector<double> D_3D;
    std::vector<double> S_2D;
    std::vector<double> S_3D;
};

}

#endif
