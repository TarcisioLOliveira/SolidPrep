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

#ifndef FIELD_MINIMUM_ENERGY_ORIENTATION_HPP
#define FIELD_MINIMUM_ENERGY_ORIENTATION_HPP

#include <map>
#include <vector>
#include "field.hpp"
#include "math/matrix.hpp"
#include "project_specification/data_map.hpp"
#include "utils/basis_tensor.hpp"

namespace field{

class MinimumEnergyOrientation : public CoordinateField {
    public:
    virtual ~MinimumEnergyOrientation() = default;
    MinimumEnergyOrientation(const projspec::DataMap& data);

    virtual void generate() override;
    virtual void initialize_views(Visualization* viz) override;
    virtual void display_views() const override;

    virtual std::array<gp_Dir, 3> get_array(const MeshElement* e, const gp_Pnt& p) const override;
    virtual math::Matrix get_matrix(const MeshElement* e, const gp_Pnt& p) const override;

    inline virtual Class get_class() const override{
        return Class::MINIMUM_ENERGY_ORIENTATION;
    }
    inline virtual SubType get_sub_type() const override{
        return SubType::DOMAIN;
    }
    inline virtual bool is_fea_dependent() const override{ return true; }

    private:

    inline math::Matrix Rtheta(const double theta) const{
        const double cos = std::cos(theta);
        const double sin = std::sin(theta);
        return math::Matrix({cos, -sin, 0,
                             sin,  cos, 0,
                               0,    0, 1}, 3, 3);
    };
    inline math::Matrix Romega(const double omega) const{
        const double cos = std::cos(omega);
        const double sin = std::sin(omega);
        return math::Matrix({ cos, 0,  sin,
                                0, 1,    0,
                             -sin, 0,  cos}, 3, 3);
    };
    inline math::Matrix Rphi(const double phi) const{
        const double cos = std::cos(phi);
        const double sin = std::sin(phi);
        return math::Matrix({1,   0,    0,
                             0, cos, -sin,
                             0, sin,  cos}, 3, 3);
    };
    inline math::Matrix Rtheta_d1(const double theta) const{
        const double cos = std::cos(theta);
        const double sin = std::sin(theta);
        return math::Matrix({-sin, -cos, 0,
                              cos, -sin, 0,
                                0,    0, 0}, 3, 3);
    };
    inline math::Matrix Romega_d1(const double omega) const{
        const double cos = std::cos(omega);
        const double sin = std::sin(omega);
        return math::Matrix({-sin, 0,  cos,
                                0, 0,    0,
                             -cos, 0, -sin}, 3, 3);
    };
    inline math::Matrix Rphi_d1(const double phi) const{
        const double cos = std::cos(phi);
        const double sin = std::sin(phi);
        return math::Matrix({0,    0,    0,
                             0, -sin, -cos,
                             0,  cos, -sin}, 3, 3);
    };

    static const bool reg;
    ProjectData* proj_data;
    const MeshElementFactory* elem_info;
    std::vector<Geometry*> geoms;
    std::vector<double> matrices;
    double thickness;
    bool show;
    size_t DIM;
    size_t NODES_PER_ELEM;
    bool first_run = true;
    bool lock = false;
    const size_t max_it;
    const bool use_strain;
    const math::Matrix R0;

    std::map<size_t, size_t> elem_pos_map;

    std::vector<double> theta, phi, omega;

    Visualization* viz;

    ViewHandler* L_view = nullptr;
    ViewHandler* R_view = nullptr;
    ViewHandler* T_view = nullptr;
};


}

#endif
