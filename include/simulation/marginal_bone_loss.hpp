/*
 *   Copyright (C) 2025 Tarcísio Ladeia de Oliveira.
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

#ifndef SIMULATION_MARGINAL_BONE_LOSS_HPP
#define SIMULATION_MARGINAL_BONE_LOSS_HPP

#include "math/matrix.hpp"
#include "meshing.hpp"
#include "project_specification/data_map.hpp"
#include "simulation.hpp"
#include "solver_manager.hpp"
#include "view_handler.hpp"

namespace simulation{

class MarginalBoneLoss : public Simulation{
    public:
    static const size_t RANGE_NUM = 4;
    typedef std::array<double, RANGE_NUM> Range;
    virtual ~MarginalBoneLoss() = default;

    MarginalBoneLoss(const projspec::DataMap& data);

    virtual void initialize_views(Visualization* viz) override;
    virtual void initialize() override;
    virtual void run() override;

    private:
    static const bool reg;
    Meshing* mesh;
    SolverManager* fem;
    ProjectData* proj;
    const Range base_range;
    const Range t, c, s;
    const double rho_eps = 0.2;
    const double lhs_eps = 1e-30;
    const Range K_e1, K_g, K_e2;
    const utils::ProblemType problem_type;
    const double MAX_RHO = 2.5;
    double rho_d = MAX_RHO/360.0;

    size_t elem_num = 0;
    size_t mand_geom_offset = 0;
    size_t density_num = 0;
    std::vector<double> v;
    double v_ref;

    const size_t geom_id;
    Geometry* mandible;
    const double time_step, time_limit;
    double dv0;
    double x0, x1, x2;
    math::Vector coeff;

    ViewHandler* stress_view = nullptr;
    ViewHandler* density_view = nullptr;
    ViewHandler* growth_view = nullptr;
    ViewHandler* strain_view = nullptr;
    ViewHandler* reactions_view = nullptr;

    Range get_K_e1(const Range& t, const Range& c) const;
    Range get_K_e2(const Range& t, const Range& c) const;
    Range get_K_g (const Range& s) const;

    inline Range interval_mult(double mult, const Range& s) const{
        Range r;
        for(size_t i = 0; i < RANGE_NUM; ++i){
            r[i] = mult*s[i];
        }
        return r;
    }

    typedef math::Vector StrainVector3D;

    double get_density_variation(const StrainVector3D& eps) const;
    double get_density_variation(const StrainVector3D& eps, const math::Matrix& D0, const math::Matrix& D) const;

    inline double LHS_3D(const size_t i, const StrainVector3D& e) const{
        const double dexy = e[0] - e[1];
        const double deyz = e[1] - e[2];
        const double dezx = e[2] - e[0];
        const double ex = e[0];
        const double ey = e[1];
        const double ez = e[2];
        const double gxy = 0.5*e[3];
        const double gyz = 0.5*e[4];
        const double gxz = 0.5*e[5];

        return std::sqrt(K_e1[i]*(dexy*dexy + deyz*deyz + dezx*dezx) + 2*K_g[i]*(gxy*gxy + gyz*gyz + gxz*gxz) + lhs_eps)
                + K_e2[i]*(ex + ey + ez);
    }

    inline double Ud(const size_t i, const StrainVector3D& eps, const math::Matrix& D0, const math::Matrix& D) const{
        const double A = D0(0,0) + 2*D0(1,2) + D0(1,1);
        const double B = 2*(D0(0,2) + D0(1,2));
        const double C = D0(2,2);

        const double e_min = t[i];
        const double e_max = c[i];

        const double K1 = A + B + C;
        const double K2 = -(B + 2*C);
        const double K3 = C;

        const double alpha = 0.5*(1.0/e_min - 1.0/e_max);
        const double se = (1/e_min + 1/e_max);
        const double beta = (se*se)/(4*(K1 + 3*K2 + 9*K3));

        const auto eps_tr = (eps[0] + eps[1] + eps[2]);
        const auto eps_h = eps_tr/3.0;
        const math::Vector eps_1({eps[0] - eps_h, eps[1] - eps_h, eps[2] - eps_h, eps[3], eps[4], eps[5]});
        const math::Vector eps_2({1.0/3.0 - alpha*eps_h, 1.0/3.0 - alpha*eps_h, 1.0/3.0 - alpha*eps_h, 0, 0, 0});

        return 3*std::sqrt(beta*(eps_1.T()*D*eps_1)) + alpha*eps_tr;
    }
};

}

#endif
