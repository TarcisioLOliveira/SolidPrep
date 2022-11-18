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

#ifndef PROJECTION_HPP
#define PROJECTION_HPP

#include <cstddef>
#include <vector>

class Projection{
    public:
    struct Parameter{
        double value;
        double final_value;
        double value_step;
        size_t iteration_step;
    };

    virtual ~Projection() = default;
    
    virtual void update(const size_t iteration) = 0;

    virtual void project_densities(std::vector<double>& new_x) const = 0;

    virtual void project_gradient(std::vector<double>& new_df, const std::vector<double>& new_x) const = 0;

    protected:
    virtual void update_parameter(Parameter& info, const size_t iteration) const;
};

#endif
