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

#ifndef DENSITY_FILTER_HPP
#define DENSITY_FILTER_HPP

#include <vector>
#include <map>
#include "logger.hpp"
#include "meshing.hpp"

class DensityFilter{
    public:
    virtual ~DensityFilter() = default;

    enum class FilterGradient{
        ELEMENTAL,
        NODAL
    };

    virtual void initialize(const Meshing* const mesh, const size_t x_size) = 0;

    virtual void filter_densities(const std::vector<double>& x, std::vector<double>& new_x) = 0;

    /**
     * Used when the function gradient can be calculated by differentiating
     * with regards to filtered density.
     */
    virtual void filter_gradient(const std::vector<double>& df, std::vector<double>& new_df){
        (void) df;
        (void) new_df; 
    }

    /**
     * Used when the function gradient must be calculated by differentiating
     * with regards to nodal densities, e.g. when the function involves using
     * density gradient or the problem uses nodal densities.
     */
    virtual void filter_gradient_nodal(const std::vector<double>& df, std::vector<double>& new_df){
        (void) df;
        (void) new_df; 
        logger::log_assert(false, logger::ERROR, "chosen density filter does not support nodal densities");
    }

    virtual const std::vector<double>& get_nodal_densities() const{
        logger::log_assert(false, logger::ERROR, "chosen density filter does not support nodal densities");
    }

    virtual void get_gradient(std::vector<double>& gradx) const{
        (void)gradx;
        logger::log_assert(false, logger::ERROR, "chosen density filter does not support nodal densities");
    }

    const inline std::map<std::pair<size_t, size_t>, size_t>& get_id_mapping() const{
        return this->id_mapping;
    }

    virtual size_t get_nodal_density_size() const{
        return 0;
    }

    protected:
    std::map<std::pair<size_t, size_t>, size_t> id_mapping;
};

#endif
