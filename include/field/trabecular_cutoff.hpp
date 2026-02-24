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

#ifndef TRABECULAR_CUTOFF_HPP
#define TRABECULAR_CUTOFF_HPP

#include <map>
#include <unordered_map>
#include <vector>
#include "field.hpp"
#include "math/matrix.hpp"
#include "project_specification/data_map.hpp"

namespace field{

class TrabecularCutoff : public ScalarField {
    public:
    class BoundaryApproximation{
        public:
        BoundaryApproximation() = default;
        BoundaryApproximation(const projspec::DataMap* const data, const TopoDS_Shape& shell);
        bool is_inside(const gp_Pnt& p) const;

        private:
        bool defined = false;
        // Parameters
        math::Vector center1, center2;
        math::Vector long_dir, v_dir, h_dir;
        // Follows POV of long_dir
        double th_top, th_bot, th_left, th_right;
        size_t pts_long, pts_v;

        double spacing_l;

        // Data
        struct DistPair{
            // Longitudinal: top, bottom;
            // Vertical: left, right;
            double d1, d2; 
        };
        std::vector<std::vector<DistPair>> bounds;
        struct VerticalBounds{
            math::Vector p1;
            double l;
            double spacing_v;
        };
        std::vector<VerticalBounds> base_points;

        inline math::Vector to_vector(const gp_Pnt& p) const{
            return {p.X(), p.Y(), p.Z()};
        }

        inline gp_Pnt to_occt_point(const math::Vector& v) const{
            return gp_Pnt(v[0], v[1], v[2]);
        }
        inline gp_Vec to_occt_vec(const math::Vector& v) const{
            return gp_Vec(v[0], v[1], v[2]);
        }
        inline gp_Dir to_occt_dir(const math::Vector& v) const{
            return gp_Dir(v[0], v[1], v[2]);
        }
    };

    virtual ~TrabecularCutoff() = default;
    TrabecularCutoff(const projspec::DataMap& data);

    virtual void generate() override;
    virtual void initialize_views(Visualization* viz) override;
    virtual void display_views() const override;
    virtual double get(const MeshElement* e, const gp_Pnt& p) const override{
        (void) p;
        return this->trabecular_map.at(e->id);
    }

    inline virtual SubType get_sub_type() const override{
        return SubType::DOMAIN;
    }
    inline virtual Class get_class() const override{
        return Class::TRABECULAR_CUTOFF;
    }

    private:
    static const bool reg;

    std::vector<Geometry*> geoms;

    ScalarField* density;
    std::unordered_map<size_t, double> trabecular_map; // 1 == trabecular; 0 == cortical

    bool show;
    const double cutoff;

    TopoDS_Shape shell;

    BoundaryApproximation inner_bounds;

    ViewHandler* distr = nullptr;
};

}

#endif
