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
#include <vector>

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
    const double pc;
    const Range t, c, s;
    const double rho_eps = 0.2;
    const double lhs_eps = 1e-30;
    const Range K_e1, K_g, K_e2;
    const utils::ProblemType problem_type;
    const std::vector<double> a0, a1, a2, a3;
    const math::Vector eps_to_x;

    const size_t geom_id;
    Geometry* mandible;
    const double time_step, maximum_volume_variation, time_limit;
    const double maturation_rate;

    ViewHandler* stress_view = nullptr;
    ViewHandler* density_view = nullptr;
    ViewHandler* growth_view = nullptr;

    Range get_K_e1(const Range& t, const Range& c) const;
    Range get_K_e2(const Range& t, const Range& c) const;
    Range get_K_g (const Range& s) const;

    typedef math::Vector StrainVector3D;

    math::Vector make_eps_to_x() const;
    double get_density_variation(const StrainVector3D& eps) const;

    inline double LHS_3D(const size_t i, const StrainVector3D& e) const{
        const double dexy = e[0] - e[1];
        const double deyz = e[1] - e[2];
        const double dezx = e[2] - e[0];
        const double ex = e[0];
        const double ey = e[1];
        const double ez = e[2];
        const double gxy = e[3];
        const double gyz = e[4];
        const double gxz = e[5];

        return std::sqrt(K_e1[i]*(dexy*dexy + deyz*deyz + dezx*dezx) + 2*K_g[i]*(gxy*gxy + gyz*gyz + gxz*gxz) + lhs_eps)
                + K_e2[i]*(ex + ey + ez);
    }
};

}

#endif
