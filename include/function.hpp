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

#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <vector>
#include <cstdlib>
#include "optimizer.hpp"

class DensityBasedFunction{
    public:

    virtual void initialize(const Optimizer* const op){(void)op;}
    virtual void update(){}
    virtual double calculate(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x) = 0;
    virtual double calculate_with_gradient(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad) = 0;
    virtual size_t additional_steps() const = 0;
};

#endif
