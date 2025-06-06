/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
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

#ifndef CONVOLUTION_HPP
#define CONVOLUTION_HPP

#include "density_filter.hpp"
#include "project_specification/data_map.hpp"

namespace density_filter{

class Convolution : public DensityFilter{
    public:
    Convolution(const projspec::DataMap& data);

    virtual ~Convolution() = default;

    virtual void initialize(const Meshing* const mesh, const size_t x_size) override;

    virtual void filter_densities(const std::vector<double>& x, std::vector<double>& new_x) override;

    virtual void filter_gradient(const std::vector<double>& df, std::vector<double>& new_df) override;

    virtual void filter_gradient_nodal(const std::vector<double>& df, std::vector<double>& new_df) override{
        (void) df;
        (void) new_df; 
        logger::log_assert(false, logger::ERROR, "chosen density filter does not support nodal densities");
    }

    virtual const std::vector<double>& get_nodal_densities() const override{
        logger::log_assert(false, logger::ERROR, "chosen density filter does not support nodal densities");
        return this->p; // Stop warning (this is unreachable)
    }

    virtual void get_gradient(std::vector<double>& gradx) const override{
        (void)gradx;
        logger::log_assert(false, logger::ERROR, "chosen density filter does not support nodal densities");
    }
    private:
    static const bool reg;
    const double radius;
    std::vector<std::vector<size_t>> neighbors;
    std::vector<double> p;
    std::vector<double> w;

    inline double get_distance(const size_t i, const size_t j) const{
        auto dx = this->p[3*i  ] - this->p[3*j  ];
        auto dy = this->p[3*i+1] - this->p[3*j+1];
        auto dz = this->p[3*i+2] - this->p[3*j+2];
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
};

}

#endif
