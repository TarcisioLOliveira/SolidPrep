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

#ifndef DENSITY_BASED_FUNCTION_COMPLIANCE_HPP
#define DENSITY_BASED_FUNCTION_COMPLIANCE_HPP

#include "function.hpp"
#include "meshing.hpp"
#include "project_specification/data_map.hpp"

namespace function::density_based{

class Compliance : public DensityBasedFunction{
    public:
    const double K_MIN = 1e-6;

    Compliance(const projspec::DataMap& data);

    virtual ~Compliance() = default;

    virtual double calculate(const DensityBasedOptimizer* const op, const std::vector<double>& u, const std::vector<double>& x) override;
    virtual double calculate_with_gradient(const DensityBasedOptimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad) override;

    private:
    static const bool reg;
    const double pc;
    const double psi;
    const Meshing* const mesh;
};

}

#endif
