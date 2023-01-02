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

#ifndef FUNCTION_GLOBAL_STRESS_PNORM_NORMALIZED_HPP
#define FUNCTION_GLOBAL_STRESS_PNORM_NORMALIZED_HPP

#include "function.hpp"
#include "meshing.hpp"
#include "finite_element.hpp"

namespace function{

class GlobalStressPnormNormalized : public DensityBasedFunction{
    public:
    GlobalStressPnormNormalized(const Meshing* const mesh, FiniteElement* fem, double pc, double P, double pt);

    virtual ~GlobalStressPnormNormalized() = default;

    virtual void initialize() override;
    virtual void update() override;
    virtual double calculate(const std::vector<double>& u, const std::vector<double>& x) override;
    virtual double calculate_with_gradient(const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad) override;

    private:
    const Meshing* const mesh;
    FiniteElement* fem;
    double pc;
    double P;
    double pt;
    std::vector<double> grad_V;
    size_t elem_number = 0;
    double c = 1;
    double alpha = 1;
    double Spn = 1;
    double Sm = 1;
};

}

#endif
