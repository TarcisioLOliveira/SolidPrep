/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
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

#ifndef MINIMAL_COMPLIANCE_HPP
#define MINIMAL_COMPLIANCE_HPP

#include "topology_optimization.hpp"
#include "density_filter.hpp"
#include "projection.hpp"

namespace topology_optimization{

class MinimalCompliance : public TopologyOptimization{
    public:
    MinimalCompliance(DensityFilter* filter, Projection* projection, ProjectData* data, double Vfinal, double xtol_abs, double ftol_rel, double result_threshold, bool save, int pc);

    virtual void initialize_views(Visualization* viz) override;
    virtual TopoDS_Shape optimize(FiniteElement* fem, Meshing* mesh) override;

    private:
    ProjectData* data;
    double Vfinal;
    double xtol_abs;
    double ftol_rel;
    double result_threshold;
    bool save_result;
    int pc;
    size_t elem_number;

    DensityFilter* filter;
    Projection* projection;
    Visualization* viz;
    FiniteElement* fem;
    Meshing* mesh;
    double max_V;
    double cur_V;
    std::vector<double> grad_V;
    double alpha;

    ViewHandler* density_view = nullptr;
    ViewHandler* disp_view = nullptr;

    double fobj(const std::vector<double>& x);
    double fobj_grad(const std::vector<double>& x, std::vector<double>& grad);

    double fc_norm(const std::vector<double>& x);
    double fc_norm_grad(const std::vector<double>& x, std::vector<double>& grad);
};

}

#endif
