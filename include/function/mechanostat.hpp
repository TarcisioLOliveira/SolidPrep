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

#ifndef FUNCTION_MECHANOSTAT_HPP
#define FUNCTION_MECHANOSTAT_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include "function.hpp"
#include "meshing.hpp"
#include "finite_element.hpp"
#include "utils.hpp"

namespace function{

class Mechanostat : public DensityBasedFunction{
    public:
    typedef std::array<double, 2> Range;
    Mechanostat(const Meshing* const mesh, FiniteElement* fem, double pc, double psiK, double beta, Range traction, Range compression, Range shear, utils::ProblemType type);

    virtual ~Mechanostat() = default;

    virtual double calculate(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x) override;
    virtual double calculate_with_gradient(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad) override;
    virtual size_t additional_steps() const override{
        return 1;
    }

    private:
    const Meshing* const mesh;
    FiniteElement* fem;
    const double beta;
    const double pc, psiK;
    const Range t, c, s;
    const double rho_eps = 0.2;
    const double lhs_eps = 1e-14;
    const Range K_e1, K_g, K_e2;
    const utils::ProblemType problem_type;

    typedef Eigen::Vector<double, 3> StrainVector2D;
    typedef Eigen::Vector<double, 6> StrainVector3D;

    inline double Hp(const double x, const double eta = 1.0) const{
        return std::atan(beta*(x - eta))/M_PI + 0.5;
    }
    inline double dHp(const double x, const double eta = 1.0) const{
        const double X = beta*(x - eta);
        const double datan = 1.0/(1.0 + X*X);
        return beta*datan/M_PI;
    }
    inline double Hm(const double x, const double eta = 1.0) const{
        return Hp(-x, -eta);
    }
    inline double dHm(const double x, const double eta = 1.0) const{
        return -dHp(-x, -eta);
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
        const double gxy = e[2];

        return std::sqrt(K_e1[i]*(dexy*dexy + deyz*deyz + dezx*dezx) + 2*K_g[i]*(gxy*gxy) + lhs_eps)
                + K_e2[i]*(ex + ey);
    }
    inline StrainVector2D dLHS_2D(const size_t i, const StrainVector2D& e) const{
        const double ex = e[0];
        const double ey = e[1];
        const double gxy = e[2];

        const double eVe = this->LHS_2D(i, e);

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
        const double gxy = e[3];
        const double gyz = e[4];
        const double gxz = e[5];

        return std::sqrt(K_e1[i]*(dexy*dexy + deyz*deyz + dezx*dezx) + 2*K_g[i]*(gxy*gxy + gyz*gyz + gxz*gxz) + lhs_eps)
                + K_e2[i]*(ex + ey + ez);
    }
    inline StrainVector3D dLHS_3D(const size_t i, const StrainVector3D& e) const{
        const double ex = e[0];
        const double ey = e[1];
        const double ez = e[2];
        const double gxy = e[3];
        const double gyz = e[4];
        const double gxz = e[5];

        const double eVe = this->LHS_3D(i, e);

        return {K_e1[i]*(2*ex - ey - ez)/eVe + K_e2[i], 
                K_e1[i]*(-ex + 2*ey - ez)/eVe + K_e2[i], 
                K_e1[i]*(-ex - ey + 2*ez)/eVe + K_e2[i], 
                2*K_g[i]*gxy/eVe,
                2*K_g[i]*gyz/eVe,
                2*K_g[i]*gxz/eVe};
    }

    // inline double I1_2D(const StrainVector2D& e) const{
    //     return e[0] + e[1];
    // }
    // inline double J2_2D(const StrainVector2D& e) const{
    //     Eigen::Matrix<double, 3, 3> V
    //         {{1, -0.5, 0},
    //          {-0.5, 1, 0},
    //          {   0, 0, 3}};

    //     return (e.transpose()*V*e).value()/3.0;
    // }
    // inline double J3_2D(const double I1, const double J2) const{
    //     return (I1*(J2 - I1*I1/9.0)/3.0);
    // }

    // inline double I1_3D(const StrainVector3D& e) const{
    //     return e[0] + e[1] + e[2];
    // }
    // inline double J2_3D(const StrainVector3D& e) const{
    //     Eigen::Matrix<double, 6, 6> V
    //         {{1, -0.5, -0.5, 0, 0, 0},
    //          {-0.5, 1, -0.5, 0, 0, 0},
    //          {-0.5, -0.5, 1, 0, 0, 0},
    //          {   0,    0, 0, 3, 0, 0},
    //          {   0,    0, 0, 0, 3, 0},
    //          {   0,    0, 0, 0, 0, 3}};

    //     return (e.transpose()*V*e).value()/3.0;
    // }
    // inline double J3_3D(const StrainVector3D& e) const{
    //     const double p = (e[0] + e[1] + e[2])/3.0;
    //     const Eigen::Matrix<double, 3, 3> s
    //         {{e[0] - p, e[3], e[4]},
    //          {e[3], e[1] - p, e[5]},
    //          {e[4], e[5], e[2] - p}};

    //     return s.determinant();
    //     //return s[0]*s[1]*s[2] + 2*e[3]*e[4]*e[5]
    //     //        - (s[0]*e[5]*e[5] + s[1]*e[4]*e[4] + s[2]*e[3]*e[3]);
    // }

    // inline double theta(const double J2, const double J3) const{
    //     return std::asin(-3*std::sqrt(3)*J3/(2*std::sqrt(J2*J2*J2)))/3;
    // }

    // inline std::array<double, 3> principal_strains(const double th, const double I1, const double J2) const{
    //     const double e1 = 2.0*std::sqrt(J2/3)*std::sin(th + 2.0*M_PI/3.0) + I1/3.0;
    //     const double e3 = 2.0*std::sqrt(J2/3)*std::sin(th - 2.0*M_PI/3.0) + I1/3.0;
    //     const double g = (e1 - e3)/2.0;

    //     return {e1, e3, g};
    // }
    // inline StrainVector2D dI1_2D() const{
    //     return {1, 1, 0};
    // }
    // inline double dJ2_2D(const StrainVector2D& e) const{
    //     Eigen::Matrix<double, 3, 3> V
    //         {{1, -0.5, 0},
    //          {-0.5, 1, 0},
    //          {   0, 0, 3}};

    //     return 2.0*(V*e).value()/3.0;
    // }
    // inline StrainVector2D dJ3_2D(const double I1, const double J2, const StrainVector2D& dI1, const StrainVector2D& dJ2) const{
    //     return (dI1*(J2 - I1*I1/9.0) + I1*(dJ2 - 2.0*I1*dI1/9.0))/3.0;
    // }

    // inline StrainVector3D dI1_3D() const{
    //     return {1, 1, 1, 0, 0, 0};
    // }
    // inline double dJ2_3D(const StrainVector3D& e) const{
    //     Eigen::Matrix<double, 6, 6> V
    //         {{1, -0.5, -0.5, 0, 0, 0},
    //          {-0.5, 1, -0.5, 0, 0, 0},
    //          {-0.5, -0.5, 1, 0, 0, 0},
    //          {   0,    0, 0, 3, 0, 0},
    //          {   0,    0, 0, 0, 3, 0},
    //          {   0,    0, 0, 0, 0, 3}};

    //     return (e.transpose()*V*e).value()/3.0;
    // }
};

}

#endif
