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

#ifndef MINIMAL_VOLUME_HPP
#define MINIMAL_VOLUME_HPP

#include "topology_optimization.hpp"

namespace topology_optimization{

class MinimalVolume : public TopologyOptimization{
    public:
    MinimalVolume(double r_o, double Smax, ProjectData* data, double rho_init, double xtol_abs, double Vfrac_abs, double result_threshold, bool save, int P, int pc);

    virtual TopoDS_Shape optimize(Visualization* viz, FiniteElement* fem, Meshing* mesh) override;

    private:
    double r_o;
    double Smax;
    ProjectData* data;
    double rho_init;
    double xtol_abs;
    double Vfrac_abs;
    double result_threshold;
    bool save_result;
    int P;
    int pc;

    Visualization* viz;
    FiniteElement* fem;
    Meshing* mesh;
    double c;
    std::vector<double> new_x;
    double max_V;
    double cur_V;
    std::vector<double> grad_V;
    double alpha;
    std::vector<std::vector<size_t>> neighbors;
    std::vector<double> p;
    std::vector<double> w;
    double Spn;
    double Sm;

    double fobj(const std::vector<double>& x);
    double fobj_grad(const std::vector<double>& x, std::vector<double>& grad);

    // Normalized stress norm (LE et al, 2009)
    void update_c();
    double fc_norm(const std::vector<double>& x);
    double fc_norm_grad(const std::vector<double>& x, std::vector<double>& grad);
};

}

#endif
