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

#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <vector>
#include <array>
#include <cblas.h>
#include <gp_Dir.hxx>
#include "math/matrix.hpp"

class MeshElement;
class Field;

class Material{
    public:
    static std::string get_name(){
        return "material";
    }
    enum Type{
        NONE,
        LINEAR_ELASTIC_ISOTROPIC,
        LINEAR_ELASTIC_ORTHOTROPIC,
        MANDIBLE,
        LINEAR_ELASTIC_ORTHOTROPIC_FIELD,
        POWER_LAW_ISOTROPIC,
        POWER_LAW_ORTHOTROPIC,
    };
    enum Class{
        ISOTROPIC,
        TRANSVERSALLY_ISOTROPIC,
        ORTHOTROPIC,
        OTHER
    };

    virtual ~Material() = default;

    /**
     * Initializes a basic material function with only maximum design stress
     * values. Accounts for anisotropy (both for positive and negative stress
     * values) if so desired.
     *
     * @param name Name of the material (for material selection).
     * @param Smax Maximum normal stresses. Input either 1, 3 or 6 values.
     * @param Tmax Maximum shear stresses. Input either 1 or 3 values.
     */
    Material(const std::string& name, std::vector<double> Smax, std::vector<double> Tmax);

    virtual math::Matrix stiffness_2D(const MeshElement* const e, const gp_Pnt& p) const = 0;
    virtual math::Matrix stiffness_3D(const MeshElement* const e, const gp_Pnt& p) const = 0;
    virtual math::Matrix stiffness_inverse_2D(const MeshElement* const e, const gp_Pnt& p) const = 0;
    virtual math::Matrix stiffness_inverse_3D(const MeshElement* const e, const gp_Pnt& p) const = 0;
    virtual double get_density(const MeshElement* const e, const gp_Pnt& p) const = 0;

    virtual double beam_E_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const = 0;
    virtual double beam_E_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const = 0;
    virtual std::array<double, 2> beam_EG_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const = 0;
    virtual std::array<double, 4> beam_EG_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const = 0;

    virtual Class get_class() const = 0;
    virtual Type get_type() const{ return this->NONE; }
    virtual bool is_homogeneous() const{ return true; }
    virtual Field* get_field() const{ return nullptr; }

    /**
     * Considers each beam node as representing maximum stresses at its cross
     * section. Therefore, returns stresses acting on the plane with normal
     * d (maximum traction, maximum compression, maximum shear, in this order).
     *
     * @param d Cross section normal
     * @returns Maximum stresses at the cross section plane
     */
    virtual std::vector<double> get_max_stresses(gp_Dir d) const = 0;

    // Unused, may delete later
    virtual double get_max_Von_Mises_2D() const;
    virtual double get_max_Von_Mises_3D() const;
    
    inline void rotate_D(math::Matrix& D, const math::Matrix& R) const{
        D = R*D*R.T();
    }
    inline void rotate_S(math::Matrix& S, const math::Matrix& R) const{
        S = R*S*R.T();
    }

    const std::string material_name;
    protected:
    std::vector<double> Smax;
    std::vector<double> Tmax;
};

#endif
