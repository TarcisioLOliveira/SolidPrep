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

#include <Eigen/Core>
#include <vector>
#include <array>
#include <cblas.h>
#include <gp_Dir.hxx>

class MeshElement;

class Material{
    public:
    enum Type{
        NONE,
        LINEAR_ELASTIC_ISOTROPIC,
        LINEAR_ELASTIC_ORTHOTROPIC,
        MANDIBLE,
        LINEAR_ELASTIC_ORTHOTROPIC_FIELD,
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

    virtual std::vector<double> stiffness_2D(const MeshElement* const e, const gp_Pnt& p) const = 0;
    virtual std::vector<double> stiffness_3D(const MeshElement* const e, const gp_Pnt& p) const = 0;
    virtual std::vector<double> stiffness_inverse_2D(const MeshElement* const e, const gp_Pnt& p) const = 0;
    virtual std::vector<double> stiffness_inverse_3D(const MeshElement* const e, const gp_Pnt& p) const = 0;
    virtual double get_density(const MeshElement* const e, const gp_Pnt& p) const = 0;

    virtual double beam_E_2D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 2, 2>& R) const = 0;
    virtual double beam_E_3D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 3, 3>& R) const = 0;
    virtual std::array<double, 2> beam_EG_2D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 2, 2>& R) const = 0;
    virtual std::array<double, 4> beam_EG_3D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 3, 3>& R) const = 0;
    virtual double S12_2D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 2, 2>& R) const = 0;
    virtual std::array<double, 2> S12_S13_3D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 3, 3>& R) const = 0;

    virtual Type get_type() const{ return this->NONE; }
    virtual bool is_homogeneous() const{ return true; }

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
    
    inline void rotate_D_base(std::vector<double>& D, const std::vector<double>& R, const size_t N) const{
        std::vector<double> Dtmp(N*N, 0);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1, D.data(), N, R.data(), N, 0, Dtmp.data(), N);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, R.data(), N, Dtmp.data(), N, 0, D.data(), N);
    }
    inline void rotate_D_2D(std::vector<double>& D, const std::vector<double>& R) const{
        return rotate_D_base(D, R, 3);
    }
    inline void rotate_D_3D(std::vector<double>& D, const std::vector<double>& R) const{
        return rotate_D_base(D, R, 6);
    }

    inline void rotate_S_base(std::vector<double>& S, const std::vector<double>& R, const size_t N) const{
        std::vector<double> Stmp(N*N, 0);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1, S.data(), N, R.data(), N, 0, Stmp.data(), N);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, R.data(), N, Stmp.data(), N, 0, S.data(), N);
    }
    inline void rotate_S_2D(std::vector<double>& S, const std::vector<double>& R) const{
        return rotate_S_base(S, R, 3);
    }
    inline void rotate_S_3D(std::vector<double>& S, const std::vector<double>& R) const{
        return rotate_S_base(S, R, 6);
    }

    const std::string name;
    protected:
    std::vector<double> Smax;
    std::vector<double> Tmax;
};

#endif
