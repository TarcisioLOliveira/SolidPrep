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

namespace topology_optimization{

class MinimalCompliance : public TopologyOptimization{
    public:
    MinimalCompliance(double r_o, ProjectData* data, double Vfinal, double xtol_abs, double result_threshold, bool save, int pc);

    virtual TopoDS_Shape optimize(Visualization* viz, FiniteElement* fem, Meshing* mesh) override;

    private:
    double r_o;
    ProjectData* data;
    double Vfinal;
    double xtol_abs;
    double result_threshold;
    bool save_result;
    int pc;
};

}

#endif
