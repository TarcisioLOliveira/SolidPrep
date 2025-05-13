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

#ifndef HEAVISIDE_HPP
#define HEAVISIDE_HPP

#include "project_specification/data_map.hpp"
#include "projection.hpp"

namespace projection{

class Heaviside : public Projection{
    public:
    Heaviside(Parameter beta, double eta);
    Heaviside(const projspec::DataMap& data);
    virtual ~Heaviside() = default;
    
    virtual void update(const size_t iteration) override;

    virtual void project_densities(std::vector<double>& new_x) const override;

    virtual void project_gradient(std::vector<double>& new_df, const std::vector<double>& new_x) const override;

    private:
    static const bool reg;
    Parameter beta;
    double eta;
};

}

#endif
