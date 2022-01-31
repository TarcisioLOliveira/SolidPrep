/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   Foobar is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Foobar is distributed in the hope that it will be useful,
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
};

}

#endif
