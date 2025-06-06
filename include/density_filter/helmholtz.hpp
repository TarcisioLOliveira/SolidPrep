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

#ifndef HELMHOLTZ_HPP
#define HELMHOLTZ_HPP

#include "density_filter.hpp"
#include "project_specification/data_map.hpp"

namespace density_filter{

class Helmholtz : public DensityFilter{
    public:
    Helmholtz(const projspec::DataMap& data);
    virtual ~Helmholtz() = default;

    virtual void initialize(const Meshing* const mesh, const size_t x_size) override;

    virtual void filter_densities(const std::vector<double>& x, std::vector<double>& new_x) override;

    virtual void filter_gradient(const std::vector<double>& df, std::vector<double>& new_df) override;

    virtual void filter_gradient_nodal(const std::vector<double>& df, std::vector<double>& new_df) override;

    virtual void get_gradient(std::vector<double>& gradx) const override;

    virtual const std::vector<double>& get_nodal_densities() const override{
        return this->nodal_densities;
    }

    virtual size_t get_nodal_density_size() const override{
        return this->nodal_densities.size();
    }

    private:
    static const bool reg;
    const double radius;
    size_t NN_n;
    size_t NN_kd;
    std::vector<double> NN;
    std::vector<double> nodal_densities;
    std::vector<double> nodal_gradient;
    std::vector<double> diff;
    std::vector<double> Hgrad;
    const Meshing* mesh;
};

}

#endif
