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

#ifndef AVERAGING_HPP
#define AVERAGING_HPP

#include "density_filter.hpp"

namespace density_filter{

class Averaging : public DensityFilter{
    public:
    virtual ~Averaging() = default;

    virtual void initialize(const Meshing* const mesh, const size_t x_size) override;

    virtual void filter_densities(const std::vector<double>& x, std::vector<double>& new_x) override;

    virtual void filter_gradient(const std::vector<double>& df, std::vector<double>& new_df) override;

    private:
    std::vector<double> D;
    std::vector<double> nodal_densities;
    std::vector<double> nodal_gradient;
    const Meshing* mesh;
};

}

#endif