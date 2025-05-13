/*
 *   Copyright (C) 2024 Tarcísio Ladeia de Oliveira.
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

#ifndef NODE_SHAPE_BASED_FUNCTION_VOLUME_HPP
#define NODE_SHAPE_BASED_FUNCTION_VOLUME_HPP

#include "function.hpp"
#include "meshing.hpp"
#include "project_specification/data_map.hpp"

namespace function::node_shape_based{

class Volume : public NodeShapeBasedFunction{
    public:
    Volume(const projspec::DataMap& data);

    virtual ~Volume() = default;

    virtual void initialize(const NodeShapeBasedOptimizer* const op) override;
    virtual double calculate(const NodeShapeBasedOptimizer* const op, const std::vector<double>& u) override;
    virtual double calculate_with_gradient(const NodeShapeBasedOptimizer* const op, const std::vector<double>& u, std::vector<double>& grad) override;

    private:
    static const bool reg;
    const Meshing* const mesh;
    double max_V;
};

}

#endif
