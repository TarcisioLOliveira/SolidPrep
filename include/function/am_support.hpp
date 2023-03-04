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

#ifndef FUNCTION_AM_SUPPORT_HPP
#define FUNCTION_AM_SUPPORT_HPP

#include <Eigen/SparseCore>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <Eigen/SparseLU>
#include "function.hpp"
#include "meshing.hpp"
#include "density_filter.hpp"
#include "projection/threshold.hpp"

namespace function{

class AMSupport : public DensityBasedFunction{
    public:
    AMSupport(const Meshing* const mesh, const DensityFilter* const filter, gp_Dir axis, double v_norm, double L, double beta, double angle);

    virtual ~AMSupport() = default;

    virtual void initialize_views(Visualization* viz) override;
    virtual void initialize(const Optimizer* const op) override;
    virtual double calculate(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x) override;
    virtual double calculate_with_gradient_nodal(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad) override;
    virtual size_t additional_steps() const override{
        return 0;
    }
    virtual DensityFilter::FilterGradient filter_gradient_type() const override{
        return DensityFilter::FilterGradient::NODAL;
    }

    private:
    const Meshing* const mesh;
    const DensityFilter* const filter;
    const gp_Dir axis;
    const double v_norm;
    const double L;
    const double beta;
    const double support_angle;
    Eigen::VectorXd b;
    Eigen::VectorXd b_grad;
    std::vector<double> diff;
    std::vector<double> Hgrad;
    std::vector<double> gradx;
    std::map<size_t, long> id_mapping;
    Eigen::SparseMatrix<double> Phi;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    projection::Threshold proj;
    bool first_time = true;
    ViewHandler* shadow_view = nullptr;

    inline double H(const double x, const double Hbeta, const double Heta, const double range, const double offset) const{
        return range*(std::tanh(Hbeta*Heta) + std::tanh(Hbeta*(x-Heta)))/(std::tanh(Hbeta*Heta) + std::tanh(Hbeta*(1.0-Heta))) - offset;
    }
    inline double dH(const double x, const double Hbeta, const double Heta, const double range) const{
        const double t = std::tanh(Hbeta*(x - Heta));
        return range*((1-t*t)*Hbeta)/(std::tanh(Hbeta*Heta) + std::tanh(Hbeta*(1.0-Heta)));
    }
    inline double H2(const double x, const double Hbeta, const double Heta) const{
        return std::tanh(Hbeta*(x - Heta));
    }
    inline double dH2(const double x, const double Hbeta, const double Heta) const{
        const double t = std::tanh(Hbeta*(x - Heta));
        return Hbeta*(1 - t*t);
    }
    inline double H3(const double x, const double Hbeta, const double Heta) const{
        return std::tanh(Hbeta*(x - Heta))/2 + 0.5;
    }
    inline double dH3(const double x, const double Hbeta, const double Heta) const{
        const double t = std::tanh(Hbeta*(x - Heta));
        return Hbeta*(1 - t*t)/2;
    }
};

}

#endif
