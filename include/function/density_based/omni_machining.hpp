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

#ifndef FUNCTION_RADIAL_MACHINING_HPP
#define FUNCTION_RADIAL_MACHINING_HPP

#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include "function.hpp"
#include "meshing.hpp"
#include "density_filter.hpp"
#include "project_specification/data_map.hpp"
#include "projection/threshold.hpp"

namespace function::density_based{

class OmniMachining : public DensityBasedFunction{
    public:
    OmniMachining(const projspec::DataMap& data);

    virtual ~OmniMachining() = default;

    virtual void initialize_views(Visualization* viz) override;
    virtual void initialize(const DensityBasedOptimizer* const op) override;
    virtual double calculate(const DensityBasedOptimizer* const op, const std::vector<double>& u, const std::vector<double>& x) override;
    virtual double calculate_with_gradient(const DensityBasedOptimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad) override;
    virtual DensityFilter::FilterGradient filter_gradient_type() const override{
        return DensityFilter::FilterGradient::ELEMENTAL;
    }

    private:
    static const bool reg;
    const Meshing* const mesh;
    const gp_Pnt center;
    const gp_Dir axis;
    const double v_norm;
    const double beta1;
    const double beta2;
    const double L;
    Eigen::VectorXd b;
    Eigen::VectorXd b_grad;
    std::vector<double> diff;
    std::vector<double> Hgrad;
    std::vector<long> id_mapping_linear;
    Eigen::SparseMatrix<double> Phi;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    projection::Threshold proj;
    bool first_time = true;
    ViewHandler* shadow_view = nullptr;
    ViewHandler* shadow_view_continuous = nullptr;
    ViewHandler* grad_view = nullptr;

    inline double heaviside(const double x, const double b, const double eta){
        return 1.0/(1.0+std::exp(-2*b*(x - eta)));
    }
    inline double heaviside_grad(const double x, const double b, const double eta){
        const double s = 1.0/(1.0+std::exp(-2*b*(x - eta)));
        const double ds = s*(1.0-s)*(2*b);
        return ds;
    }
};

}

#endif
