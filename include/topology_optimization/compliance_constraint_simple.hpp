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

#ifndef COMPLIANCE_CONSTRAINT_SIMPLE_HPP
#define COMPLIANCE_CONSTRAINT_SIMPLE_HPP

#include "topology_optimization.hpp"
#include "density_filter.hpp"
#include "projection.hpp"

namespace topology_optimization{

class ComplianceConstraintSimple : public TopologyOptimization {
    public:

    ComplianceConstraintSimple(DensityFilter* filter, Projection* projection, double c_max, ProjectData* data, double rho_init, double xtol_abs, double result_threshold, bool save, int P, int pc);

    virtual void initialize_views(Visualization* viz) override;
    virtual TopoDS_Shape optimize(FiniteElement* fem, Meshing* mesh) override;

    private:
    double c_max;
    ProjectData* data;
    double rho_init;
    double xtol_abs;
    double result_threshold;
    bool save_result;
    int P;
    int pc;

    DensityFilter* filter;
    Projection* projection;
    Visualization* viz;
    FiniteElement* fem;
    Meshing* mesh;
    size_t elem_number;
    std::vector<double> grad_V;
    std::vector<double> u;

    ViewHandler* stress_view = nullptr;
    ViewHandler* density_view = nullptr;

    double fobj_grad(const std::vector<double>& x, std::vector<double>& grad);

    double fc_grad(const std::vector<double>& x, std::vector<double>& grad);
};

}

#endif
