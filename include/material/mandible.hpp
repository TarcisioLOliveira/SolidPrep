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

#ifndef MATERIAL_MANDIBLE_HPP
#define MATERIAL_MANDIBLE_HPP

#include "material.hpp"

namespace material{

class Mandible : public Material{
    public:
    class ImplantRegion{
        public:
        gp_Pnt center_1;
        gp_Pnt center_2;
        double r1;
        double r2;

        void initialize(double decay_distance, double str_pnt, double str_pnt_dist);
        double get_implant_multiplier(const gp_Pnt& p) const;

        private:
        inline double f(const double x) const{
            return a*x*x*x + b*x*x + c*x + d;
        }

        gp_Dir normal;
        double decay_distance;
        double min_str;
        double a, b, c, d;
        double max_l;

    };

    Mandible(const std::string& name, Material* outer, Material* inner, const std::string& path_points1, const std::string& path_points2, double C, bool with_implant = false, ImplantRegion imp = ImplantRegion());

    virtual std::vector<double> stiffness_2D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual std::vector<double> stiffness_3D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual std::vector<double> stiffness_inverse_2D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual std::vector<double> stiffness_inverse_3D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual double get_density(const MeshElement* const e, const gp_Pnt& p) const override;

    virtual double beam_E_2D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 2, 2>& R) const override;
    virtual double beam_E_3D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 3, 3>& R) const override;
    virtual std::array<double, 2> beam_EG_2D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 2, 2>& R) const override;
    virtual std::array<double, 4> beam_EG_3D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 3, 3>& R) const override;
    virtual double S12_2D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 2, 2>& R) const override;
    virtual std::array<double, 2> S12_S13_3D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 3, 3>& R) const override;

    virtual Type get_type() const override{ return this->MANDIBLE; }
    virtual bool is_homogeneous() const override{ return false; }

    virtual std::vector<double> get_max_stresses(gp_Dir d) const override{
        (void)d;
        return std::vector<double>{};
    }

    private:
    struct RingPoint{
        double r;
        double theta;
        double lambda;

        inline bool operator<(const RingPoint& p){
            return this->theta < p.theta;
        }
    };
    struct RingPointCartesian{
        gp_Pnt p;
        double theta;
        double lambda;

        inline bool operator<(const RingPointCartesian& p){
            return this->theta < p.theta;
        }
    };
    struct Ring{
        Ring(Mandible* M):M(M){}
        const Mandible* const M;
        gp_Dir normal;
        gp_Pnt center;
        std::vector<RingPointCartesian> points;

        void initialize_center(const std::vector<gp_Pnt>& init);
        void initialize_points(const std::vector<gp_Pnt>& init);
        RingPoint get_r_max(double theta) const;
    };

    inline double heaviside(const double r, const double r_max) const{
        return std::tanh(-C*(r - r_max))/2.0 + 0.5;
    }

    std::vector<gp_Pnt> load_points(const std::string& path) const;
    RingPoint to_ring_point(const gp_Pnt& p) const;
    RingPointCartesian to_ring_point_cartesian(const gp_Pnt& p) const;
    double get_multiplier(const gp_Pnt& p) const;
    double angle_with_ref(const gp_Vec& dist) const;

    const Material* const outer;
    const Material* const inner;
    const double C;
    Ring ring1, ring2;
    gp_Dir center_normal;
    gp_Dir center_ref;
    const bool with_implant;
    ImplantRegion implant;
};

}

#endif
