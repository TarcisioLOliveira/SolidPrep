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

#ifndef IMPLANT_REGION_HPP
#define IMPLANT_REGION_HPP

#include <vector>
#include "field.hpp"
#include "project_specification/data_map.hpp"

namespace field{

class ImplantRegion : public ScalarField {
    public:
    virtual ~ImplantRegion() = default;
    ImplantRegion(const projspec::DataMap& data);

    virtual void generate() override;
    virtual void initialize_views(Visualization* viz) override;
    virtual void display_views() const override;
    virtual double get(const MeshElement* e, const gp_Pnt& p) const override;

    inline virtual SubType get_sub_type() const override{
        return SubType::DOMAIN;
    }
    inline virtual Class get_class() const override{
        return Class::IMPLANT_REGION;
    }

    void freeze(std::vector<double>& current_values, std::vector<double>& maximum_values);

    inline void set_values(const std::vector<double>& new_vals){
        if(this->frozen){
            std::copy(new_vals.begin(), new_vals.end(), this->frozen_values.begin());
        }
    }

    private:
    static const bool reg;
    inline double f(const double x) const{
        double sum = a[0];
        for(size_t i = 1; i < this->a_len; ++i){
            sum += a[i]*std::pow(x, i);
        }
        return sum;
    }
    double get_implant_multiplier(const gp_Pnt& p) const;

    const MeshElementFactory* elem_info;
    std::vector<Geometry*> geoms;
    double thickness;
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
    ScalarField* density_field;

    // Value freezing
    bool frozen = false;
    std::vector<double> frozen_values;
    std::map<size_t, size_t> elem_id_to_value;
    //

    Visualization* viz;

    bool show;

    ViewHandler* density = nullptr;
};

}

#endif
