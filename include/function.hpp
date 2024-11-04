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
#include "density_filter.hpp"
#include "visualization.hpp"

class Optimizer;
class DensityBasedOptimizer;
class NodeShapeBasedOptimizer;


class Function {
    public:
    virtual ~Function() = default;

    virtual void initialize_views(Visualization* viz){(void)viz;}
    virtual void update(){}
};

class DensityBasedFunction : public Function {
    public:
    virtual ~DensityBasedFunction() = default;

    virtual void initialize(const DensityBasedOptimizer* const op){(void)op;}
    virtual double calculate(const DensityBasedOptimizer* const op, const std::vector<double>& u, const std::vector<double>& x){
        (void)op;
        (void)u;
        (void)x;
        return 0;
    }
    virtual double calculate_with_gradient(const DensityBasedOptimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
        (void)op;
        (void)u;
        (void)x;
        (void)grad;
        return 0;
    }
    virtual double calculate_with_gradient_nodal(const DensityBasedOptimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
        (void)op;
        (void)u;
        (void)x;
        (void)grad;
        return 0;
    }

    virtual DensityFilter::FilterGradient filter_gradient_type() const{
        return DensityFilter::FilterGradient::ELEMENTAL;
    }
};

class NodeShapeBasedFunction : public Function {
    public:
    virtual ~NodeShapeBasedFunction() = default;

    virtual void initialize(const NodeShapeBasedOptimizer* const op){(void)op;}
    virtual double calculate(const NodeShapeBasedOptimizer* const op, const std::vector<double>& u){
        (void)op;
        (void)u;
        return 0;
    }
    virtual double calculate_with_gradient(const NodeShapeBasedOptimizer* const op, const std::vector<double>& u, std::vector<double>& grad){
        (void)op;
        (void)u;
        (void)grad;
        return 0;
    }

};

#endif
