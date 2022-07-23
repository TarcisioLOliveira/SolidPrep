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
        MinimalCompliance(double r_o, ProjectData* data, double Vfinal, double xtol_abs, double ftol_rel, double result_threshold, bool save, int pc);

    virtual TopoDS_Shape optimize(Visualization* viz, FiniteElement* fem, Meshing* mesh) override;

    private:
    double r_o;
    ProjectData* data;
    double Vfinal;
    double xtol_abs;
    double ftol_rel;
    double result_threshold;
    bool save_result;
    int pc;
    size_t elem_number;

    Visualization* viz;
    FiniteElement* fem;
    Meshing* mesh;
    double max_V;
    double cur_V;
    std::vector<double> grad_V;
    double alpha;
    std::vector<std::vector<size_t>> neighbors;
    std::vector<double> p;
    std::vector<double> w;

    double fobj(const std::vector<double>& x);
    double fobj_grad(const std::vector<double>& x, std::vector<double>& grad);

    double fc_norm(const std::vector<double>& x);
    double fc_norm_grad(const std::vector<double>& x, std::vector<double>& grad);

    void init_convolution_filter(size_t x_size);
    std::vector<double> convolution_filter_density(const std::vector<double>& x);
    std::vector<double> convolution_grad_correction(const std::vector<double>& df);

    inline double get_distance(const size_t i, const size_t j) const{
        auto dx = this->p[3*i  ] - this->p[3*j+0];
        auto dy = this->p[3*i+1] - this->p[3*j+1];
        auto dz = this->p[3*i+2] - this->p[3*j+2];
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
};

}

#endif
