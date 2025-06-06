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
#include "project_specification/data_map.hpp"
#include "utils/delayed_pointer.hpp"

namespace material{

class Mandible : public Material{
    public:
    class ImplantRegion{
        public:
        ImplantRegion() = default;
        ImplantRegion(const projspec::DataMap* const data);

        double get_implant_multiplier(const gp_Pnt& p) const;
        void set_maturation_alpha(double alpha);

        private:
        inline double f(const double x) const{
            double sum = 0;
            for(size_t i = 0; i < this->a_len; ++i){
                sum += a[i]*std::pow(x, i);
            }
            return sum;
        }

        gp_Pnt center_1;
        gp_Pnt center_2;
        double r1;
        double r2;
        gp_Dir normal;
        double decay_distance;
        std::vector<double> a;
        double min_str;
        std::vector<double> a_orig;
        size_t a_len;
        double max_l;
    };

    Mandible(const projspec::DataMap& data);

    virtual math::Matrix stiffness_2D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual math::Matrix stiffness_3D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual math::Matrix stiffness_inverse_2D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual math::Matrix stiffness_inverse_3D(const MeshElement* const e, const gp_Pnt& p) const override;
    virtual double get_density(const MeshElement* const e, const gp_Pnt& p) const override;

    virtual double beam_E_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override;
    virtual double beam_E_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override;
    virtual std::array<double, 2> beam_EG_2D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override;
    virtual std::array<double, 4> beam_EG_3D(const MeshElement* const e, const gp_Pnt& p, const math::Matrix& R) const override;

    virtual Type get_type() const override{ return this->MANDIBLE; }
    virtual bool is_homogeneous() const override{ return false; }
    virtual Field* get_field() const override{
        // Temporary, should return both fields, actually.
        return this->inner->get_field();
    }

    virtual std::vector<double> get_max_stresses(gp_Dir d) const override{
        (void)d;
        return std::vector<double>{};
    }

    inline void set_maturation_alpha(double alpha){
        this->implant.set_maturation_alpha(alpha);
    }

    private:
    static const bool reg;
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

    const utils::DelayedPointerView<Material> outer;
    const utils::DelayedPointerView<Material> inner;
    const double C;
    Ring ring1, ring2;
    gp_Dir center_normal;
    gp_Dir center_ref;
    const bool with_implant;
    ImplantRegion implant;
};

}

#endif
