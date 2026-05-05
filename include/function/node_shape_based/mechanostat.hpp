/*
 *   Copyright (C) 2024 Tarcísio Ladeia de Oliveira.
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

#ifndef NODE_SHAPE_BASED_FUNCTION_MECHANOSTAT_HPP
#define NODE_SHAPE_BASED_FUNCTION_MECHANOSTAT_HPP

#include "function.hpp"
#include "math/matrix.hpp"
#include "meshing.hpp"
#include "project_specification/data_map.hpp"
#include "solver_manager.hpp"

namespace function::node_shape_based{

class Mechanostat: public NodeShapeBasedFunction{
    public:
    typedef std::array<double, 2> Range;
    Mechanostat(const projspec::DataMap& data);

    virtual void initialize_views(Visualization* viz) override;
    virtual void initialize(const NodeShapeBasedOptimizer* const op) override;
    virtual double calculate(const NodeShapeBasedOptimizer* const op, const std::vector<double>& u) override;
    virtual double calculate_with_gradient(const NodeShapeBasedOptimizer* const op, const std::vector<double>& u, std::vector<double>& grad) override;

    private:
    static const bool reg;
    const Meshing* const mesh;
    SolverManager* fem;
    const double beta;
    const Range t, c, s;
    const double rho_eps = 0.2;
    const double lhs_eps = 1e-30;
    const Range K_e1, K_g, K_e2;
    const utils::ProblemType problem_type;
    std::vector<double> He;
    std::vector<double> gradHe;
    ViewHandler* shadow_view = nullptr;
    ViewHandler* gradient_view = nullptr;

    typedef math::Vector StrainVector2D;
    typedef math::Vector StrainVector3D;

    inline double Hp(const double x, const double eta = 1.0) const{
        return std::atan(beta*(x - eta))/M_PI + 0.5;
    }
    inline double dHp(const double x, const double eta = 1.0) const{
        const double X = beta*(x - eta);
        const double datan = 1.0/(1.0 + X*X);
        return beta*datan/M_PI;
    }
    inline double Hm(const double x, const double eta = 1.0) const{
        (void) x;
        (void) eta;
        return 0;
        //return Hp(-x, -eta);
    }
    inline double dHm(const double x, const double eta = 1.0) const{
        (void) x;
        (void) eta;
        return 0;
        //return -dHp(-x, -eta);
    }

    inline double relaxed_rho(const double x) const{
        return x/(rho_eps*(1.0 - x) + x);
    }
    inline double relaxed_rho_grad(const double x) const{
        const double d = (rho_eps*(1.0 - x) + x);
        return rho_eps/(d*d);
    }

    inline double LHS_2D(const size_t i, const StrainVector2D& e) const{
        const double dexy = e[0] - e[1];
        const double deyz = e[1];
        const double dezx = e[0];
        const double ex = e[0];
        const double ey = e[1];
        const double gxy = e[2]/2.0;

        return std::sqrt(K_e1[i]*(dexy*dexy + deyz*deyz + dezx*dezx) + 2*K_g[i]*(gxy*gxy) + lhs_eps)
                + K_e2[i]*(ex + ey);
    }
    inline StrainVector2D dLHS_2D(const size_t i, const StrainVector2D& e) const{
        const double dexy = e[0] - e[1];
        const double deyz = e[1];
        const double dezx = e[0];
        const double ex = e[0];
        const double ey = e[1];
        const double gxy = e[2]/2.0;

        const double eVe = std::sqrt(K_e1[i]*(dexy*dexy + deyz*deyz + dezx*dezx) + 2*K_g[i]*(gxy*gxy) + lhs_eps);

        return {K_e1[i]*(2*ex - ey)/eVe + K_e2[i], 
                K_e1[i]*(-ex + 2*ey)/eVe + K_e2[i], 
                2*K_g[i]*gxy/eVe};
    }
    inline double LHS_3D(const size_t i, const StrainVector3D& e) const{
        const double dexy = e[0] - e[1];
        const double deyz = e[1] - e[2];
        const double dezx = e[2] - e[0];
        const double ex = e[0];
        const double ey = e[1];
        const double ez = e[2];
        const double gxy = e[3]/2.0;
        const double gyz = e[4]/2.0;
        const double gxz = e[5]/2.0;

        return std::sqrt(K_e1[i]*(dexy*dexy + deyz*deyz + dezx*dezx) + 2*K_g[i]*(gxy*gxy + gyz*gyz + gxz*gxz) + lhs_eps)
                + K_e2[i]*(ex + ey + ez);
    }
    inline StrainVector3D dLHS_3D(const size_t i, const StrainVector3D& e) const{
        const double dexy = e[0] - e[1];
        const double deyz = e[1] - e[2];
        const double dezx = e[2] - e[0];
        const double ex = e[0];
        const double ey = e[1];
        const double ez = e[2];
        const double gxy = e[3]/2.0;
        const double gyz = e[4]/2.0;
        const double gxz = e[5]/2.0;

        const double eVe = std::sqrt(K_e1[i]*(dexy*dexy + deyz*deyz + dezx*dezx) + 2*K_g[i]*(gxy*gxy + gyz*gyz + gxz*gxz) + lhs_eps);

        return {K_e1[i]*(2*ex - ey - ez)/eVe + K_e2[i], 
                K_e1[i]*(-ex + 2*ey - ez)/eVe + K_e2[i], 
                K_e1[i]*(-ex - ey + 2*ez)/eVe + K_e2[i], 
                2*K_g[i]*gxy/eVe,
                2*K_g[i]*gyz/eVe,
                2*K_g[i]*gxz/eVe};
    }

    inline double Ud_3D(const size_t i, const StrainVector3D& eps, const math::Matrix& D0, const math::Matrix& D) const{
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
        //const math::Vector eps_2({1.0/3.0 - alpha*eps_h, 1.0/3.0 - alpha*eps_h, 1.0/3.0 - alpha*eps_h, 0, 0, 0});

        return 3*std::sqrt(beta*(eps_1.T()*D*eps_1)) + alpha*eps_tr;
    }

    inline StrainVector3D dUd_3D(const size_t i, const StrainVector3D& eps, const math::Matrix& D0, const math::Matrix& D) const{
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
        //const math::Vector eps_2({1.0/3.0 - alpha*eps_h, 1.0/3.0 - alpha*eps_h, 1.0/3.0 - alpha*eps_h, 0, 0, 0});

        const math::Vector k({1, 1, 1, 0, 0, 0});
        const auto eDe = std::sqrt(beta*(eps_1.T()*D*eps_1));

        const math::Matrix I = math::Matrix::identity(6);
        const math::Matrix k2 = (1.0/3.0)*math::Matrix({
                1, 1, 1, 0, 0, 0,
                1, 1, 1, 0, 0, 0,
                1, 1, 1, 0, 0, 0,
                0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0
                }, 6, 6);

        return 3*beta*((I - k2)*D*eps_1)/eDe + alpha*k;
    }
};

}

#endif
