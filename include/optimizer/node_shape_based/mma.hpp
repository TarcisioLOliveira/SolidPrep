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

#ifndef OPTIMIZER_NODE_SHAPE_BASED_MMA_HPP
#define OPTIMIZER_NODE_SHAPE_BASED_MMA_HPP

#include <memory>
#include "optimizer.hpp"
#include "function.hpp"

namespace optimizer::node_shape_based{

class MMA : public NodeShapeBasedOptimizer{
    public:
    MMA(ShapeHandler sh, ProjectData* data, std::vector<std::unique_ptr<NodeShapeBasedFunction>> objective, std::vector<double> objective_weights, std::vector<NodeShapeBasedConstraint> constraints, double asyminit, double asymdec, double asyminc, double minfac, double maxfac, double c, double xtol_abs, double ftol_rel, bool save);

    virtual void initialize_views(Visualization* viz) override;
    virtual TopoDS_Shape optimize(SolverManager* fem, Meshing* mesh) override;

    private:
    ProjectData* data;
    const double xtol_abs;
    const double ftol_rel;
    const double asyminit, asymdec, asyminc;
    const double minfac, maxfac;
    const double c;
    const bool save_result;
    std::vector<std::unique_ptr<NodeShapeBasedFunction>> objective;
    std::vector<double> objective_weights;
    std::vector<NodeShapeBasedConstraint> constraints;

    Visualization* viz;

    ViewHandler* stress_view = nullptr;
};

}

#endif
