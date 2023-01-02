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

#ifndef FUNCTION_VOLUME_HPP
#define FUNCTION_VOLUME_HPP

#include "function.hpp"
#include "meshing.hpp"

namespace function{

class Volume : public DensityBasedFunction{
    public:
    Volume(const Meshing* const mesh);

    virtual ~Volume() = default;

    virtual void initialize() override;
    virtual double calculate(const std::vector<double>& u, const std::vector<double>& x) override;
    virtual double calculate_with_gradient(const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad) override;

    private:
    const Meshing* const mesh;
    std::vector<double> grad_V;
    double max_V;
};

}

#endif
