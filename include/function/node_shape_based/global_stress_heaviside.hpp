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

#ifndef NODE_SHAPE_BASED_FUNCTION_GLOBAL_STRESS_HEAVISIDE_HPP
#define NODE_SHAPE_BASED_FUNCTION_GLOBAL_STRESS_HEAVISIDE_HPP

#include "function.hpp"
#include "meshing.hpp"
#include "solver_manager.hpp"

namespace function::node_shape_based{

class GlobalStressHeaviside : public NodeShapeBasedFunction{
    public:
    GlobalStressHeaviside(const Meshing* const mesh, SolverManager* fem, double max_stress, double C);

    virtual ~GlobalStressHeaviside() = default;

    virtual double calculate(const NodeShapeBasedOptimizer* const op, const std::vector<double>& u) override;
    virtual double calculate_with_gradient(const NodeShapeBasedOptimizer* const op, const std::vector<double>& u, std::vector<double>& grad) override;

    private:
    const Meshing* const mesh;
    SolverManager* fem;
    const double max_stress;
    const double C;
    const double vm_eps = 1e-14;

    inline double heaviside(const double x) const{
        return std::atan(C*(x - max_stress))/M_PI + 0.5;
    }
    inline double heaviside_grad(const double x) const{
        const double X = C*(x - max_stress);
        const double datan = 1.0/(1.0 + X*X);
        return C*datan/M_PI;
    }
};

}

#endif
