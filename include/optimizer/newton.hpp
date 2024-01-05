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

#ifndef OPTIMIZER_NEWTON_HPP
#define OPTIMIZER_NEWTON_HPP

#include "density_filter.hpp"
#include "projection.hpp"
#include "optimizer.hpp"
#include "function.hpp"
#include <memory>

namespace optimizer{

class Newton : public Optimizer{
    public:
    Newton(DensityFilter* filter, Projection* projection, ProjectData* data, std::vector<Constraint> functions, double pc, double psi, double rho_init, double xtol_abs, double ftol_rel, double result_threshold, bool save);

    virtual void initialize_views(Visualization* viz) override;
    virtual TopoDS_Shape optimize(SolverManager* fem, Meshing* mesh) override;

    private:
    ProjectData* data;
    const double rho_init;
    const double xtol_abs;
    const double ftol_rel;
    const double pc;
    const double psi;
    const double result_threshold;
    const bool save_result;
    std::vector<Constraint> functions;

    DensityFilter* filter;
    Projection* projection;
    Visualization* viz;

    ViewHandler* stress_view = nullptr;
    ViewHandler* density_view = nullptr;
};

}

#endif
