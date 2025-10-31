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

#ifndef CTRI3_HPP
#define CTRI3_HPP

#include "element.hpp"
#include "logger.hpp"
#include "material.hpp"
#include <vector>
#include "math/matrix.hpp"
#include "utils.hpp"

namespace contact_element{

class CTRI3 : public ContactMeshElement{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 2;
    static const size_t NODE_DOF       = 3;
    static const size_t NODES_PER_ELEM = 3;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;
    static const size_t S_SIZE         = 6; // Size of the stress and strain vectors
    static const size_t DIM            = 3; // Number of dimensions
    static const Element::Shape SHAPE_TYPE = Element::Shape::TRI;
    static const utils::ProblemType PROBLEM_TYPE = utils::PROBLEM_TYPE_3D;

    CTRI3(ElementShape s, const MeshElement* const e1, const MeshElement* const e2, bool e1_base);

    virtual math::Matrix get_frictionless_Ge() const override;

    virtual math::Matrix fl3_uL(const math::Matrix& D, const math::Vector& ln_e) const override;
    virtual math::Matrix fl3_LL(const math::Matrix& D, const math::Vector& ln_e, const math::Vector& lp1_e, const math::Vector& lp2_e, const std::vector<double>& u) const override;
    virtual void fl3_Ku(const math::Matrix& D, const std::vector<long> u_pos, const std::vector<long>& lu_pos, const std::vector<double>& u, std::vector<double>& Ku) const override;
    virtual void fl3_dKu(const math::Matrix& D, const double eta, const std::vector<long> u_pos, const std::vector<long>& lu_pos, const std::vector<double>& u, const std::vector<double>& du, std::vector<double>& Ku) const override;
    virtual math::Vector fl3_eq(const math::Vector& ln_e, const math::Vector& lp1_e, const math::Vector& lp2_e, const math::Vector& u_e) const override;
    virtual math::Vector fl3_eq(const math::Vector& ln_e, const math::Vector& lp1_e, const math::Vector& lp2_e, const math::Vector& u_e, const size_t dof) const override;

    virtual math::Matrix fl2_uu(const math::Vector& l_e, const math::Vector& u1, const math::Vector& u2) const override;
    virtual math::Matrix fl2_uL(const math::Vector& l_e, const math::Vector& u1, const math::Vector& u2) const override;
    virtual math::Matrix fl2_LL(const math::Vector& l_e, const math::Vector& u1, const math::Vector& u2) const override;

    virtual void fl2_Ku_lambda(const double EPS, const std::vector<long> u1_pos, const std::vector<long> u2_pos, const std::vector<long>& lu_pos, const std::vector<double>& u, std::vector<double>& Ku) const override;
    virtual void fl2_dKu_lambda(const double EPS, const std::vector<long> u1_pos, const std::vector<long> u2_pos, const std::vector<long>& lu_pos, const std::vector<double>& u, const std::vector<double>& du, std::vector<double>& Ku) const override;

    virtual double get_area() const override{
        return this->delta;
    }

    virtual gp_Pnt get_centroid() const override{
        const size_t N = CTRI3::NODES_PER_ELEM;

        double x = 0;
        double y = 0;
        double z = 0;
        for(size_t i = 0; i < N; ++i){
            const auto& n = this->nodes[i];
            x += n->point.X();
            y += n->point.Y();
            z += n->point.Z();
        }

        return gp_Pnt(x/N, y/N, z/N);
    }
    virtual gp_Dir get_normal() const override{
        return gp_Dir(R(0,2), R(1,2), R(2,2));
    }
    virtual const math::Matrix& get_R() const override{
        return this->R;
    }

    protected:
    double a[3], b[3], c[3], delta;
    size_t E_KW;
    math::Matrix R;
    const double FL2_EPS = 1e-14;//1e-14;
    const double FL2_OFF = 1e-5;
    const double FL2_K = 1e5;

    double A(const double x) const{
        return (1 + x/std::sqrt(x*x + FL2_EPS))/2;
    }
    double dA(const double x) const{
        return FL2_EPS/(2*std::pow(x*x + FL2_EPS, 1.5));
    }
    double ddA(const double x) const{
        return -3*FL2_EPS*x/(2*std::pow(x*x + FL2_EPS, 2.5));
    }

    double h(const double x) const{
        //return A(x)*(x + FL2_OFF);
        const double dx = x - FL2_OFF;
        if(std::abs(FL2_K*dx) < 345){
            return std::log(1.0 + std::exp(FL2_K*dx))/FL2_K;
        } else if(x > 0){
            return x;
        } else {
            return 0;
        }
    }

    double dh(const double x) const{
        //return dA(x)*(x + FL2_OFF) + A(x);
        const double dx = x - FL2_OFF;
        if(std::abs(FL2_K*dx) < 345){
            return 1.0/(1.0 + std::exp(-FL2_K*dx));
        } else if(x > 0){
            return 1;
        } else {
            return 0;
        }
    }

    double ddh(const double x) const{
        //return ddA(x)*(x + FL2_OFF) + 2*dA(x);
        const double dx = x - FL2_OFF;
        if(std::abs(FL2_K*dx) < 345){
            const double ekx = std::exp(-FL2_K*dx);
            const double ekx1 = 1.0 + ekx;
            return FL2_K*ekx/(ekx1*ekx1);
        } else {
            return 0;
        }
    }

    virtual gp_Dir get_d1() const{
        return gp_Dir(R(0,1), R(1,1), R(2,1));
    }
    virtual gp_Dir get_d2() const{
        return gp_Dir(R(0,0), R(1,0), R(2,0));
    }

    inline gp_Pnt R_GS_point(double c1, double c2, double c3) const{
        const size_t N = CTRI3::NODES_PER_ELEM;
        std::array<double, N> x, y;//, z;
        std::fill(x.begin(), x.end(), 0);
        std::fill(y.begin(), y.end(), 0);
        for(size_t i = 0; i < N; ++i){
            for(size_t j = 0; j < NODE_DOF; ++j){
                x[i] += R(j,0)*this->nodes[i]->point.Coord(1+j);
                y[i] += R(j,1)*this->nodes[i]->point.Coord(1+j);
            }
        }
        return gp_Pnt(
            c1*x[0] + c2*x[1] + c3*x[2],
            c1*y[0] + c2*y[1] + c3*y[2],
            0
        );
    }
    inline gp_Pnt GS_point(double c1, double c2, double c3) const{
        return gp_Pnt(
            c1*this->nodes[0]->point.X() + c2*this->nodes[1]->point.X() + c3*this->nodes[2]->point.X(),
            c1*this->nodes[0]->point.Y() + c2*this->nodes[1]->point.Y() + c3*this->nodes[2]->point.Y(),
            c1*this->nodes[0]->point.Z() + c2*this->nodes[1]->point.Z() + c3*this->nodes[2]->point.Z()
        );
    }
    inline double N(const gp_Pnt& p, size_t i) const{
        return a[i] + b[i]*p.X() + c[i]*p.Y();
    }
    inline double dNdx(size_t i) const{
        return b[i];
    }
    inline double dNdy(size_t i) const{
        return c[i];
    }
    inline math::Vector N_mat_1dof(const gp_Pnt& p) const{
        return math::Vector{N(p, 0), N(p, 1), N(p, 2)};
    }
    inline math::Vector dNdx_mat_1dof() const{
        return math::Vector{dNdx(0), dNdx(1), dNdx(2)};
    }
    inline math::Vector dNdy_mat_1dof() const{
        return math::Vector{dNdy(0), dNdy(1), dNdy(2)};
    }
    inline math::Matrix gradN_1dof() const{
        return math::Matrix(
            {dNdx(0), dNdx(1), dNdx(2),
             dNdy(0), dNdy(1), dNdy(2)}, 2, 3);
    }
    inline math::Vector gradN_1dof_vec(const math::Vector& dXI) const{
        return math::Vector
            {dNdx(0)*dXI[0] + dNdy(0)*dXI[1],
             dNdx(1)*dXI[0] + dNdy(1)*dXI[1],
             dNdx(2)*dXI[0] + dNdy(2)*dXI[1]};
    }

    inline math::Matrix eps_mat_lambda(const gp_Dir& n) const{
        math::Matrix eps(S_SIZE, NODES_PER_ELEM);
        std::vector<math::Matrix> didj(NODES_PER_ELEM, math::Matrix(NODE_DOF, NODE_DOF));

        math::Matrix Rl(
                {R(0,0), R(1,0), R(2,0),
                 R(0,1), R(1,1), R(2,1)}, 2, 3);

        // Calculate the gradient of lambda times n, without the lambda vector,
        // then transform the rank-3 tensor into the strain matrix, which
        // multiplied to the lambda vector should result in a strain vector.
        math::Matrix gradN(Rl.T()*this->gradN_1dof());

        // for(size_t i = 0; i < NODE_DOF; ++i){
        //     const math::Matrix gradNni(n.Coord(1+i)*gradN);
        //     for(size_t j = 0; j < NODE_DOF; ++j){
        //         for(size_t k = 0; k < NODES_PER_ELEM; ++k){
        //             didj[k](i, j) = gradNni(j, k);
        //         }
        //     }
        // }

        // for(size_t i = 0; i < NODE_DOF; ++i){
        //     for(size_t j = 0; j < NODES_PER_ELEM; ++j){
        //         eps(i, j) = didj[j](i,i);
        //     }
        // }
        // for(size_t j = 0; j < NODES_PER_ELEM; ++j){
        //     eps(3, j) = didj[j](0, 1) + didj[j](1, 0);
        // }
        // for(size_t j = 0; j < NODES_PER_ELEM; ++j){
        //     eps(4, j) = didj[j](0, 2) + didj[j](2, 0);
        // }
        // for(size_t j = 0; j < NODES_PER_ELEM; ++j){
        //     eps(5, j) = didj[j](1, 2) + didj[j](2, 1);
        // }

        for(size_t i = 0; i < NODE_DOF; ++i){
            for(size_t j = 0; j < NODES_PER_ELEM; ++j){
                eps(i, j) = gradN(i, j)*n.Coord(1+i);
            }
        }
        for(size_t j = 0; j < NODES_PER_ELEM; ++j){
            eps(3, j) = gradN(0, j)*n.Y() + gradN(1, j)*n.X();
        }
        for(size_t j = 0; j < NODES_PER_ELEM; ++j){
            eps(4, j) = gradN(0, j)*n.Z() + gradN(2, j)*n.X();
        }
        for(size_t j = 0; j < NODES_PER_ELEM; ++j){
            eps(5, j) = gradN(1, j)*n.Z() + gradN(2, j)*n.Y();
        }

        return eps;
    }

    inline std::array<math::Matrix, 3> square_lambda_deriv(const gp_Pnt& p) const{
        std::array<math::Matrix, 3> result;

        const auto N = this->N_mat_1dof(p);
        math::Matrix Rl(
                {R(0,0), R(1,0), R(2,0),
                 R(0,1), R(1,1), R(2,1)}, 2, 3);

        math::Matrix gradN(Rl.T()*this->gradN_1dof());
        for(size_t i = 0; i < 3; ++i){
            math::Vector dN({gradN(i, 0), gradN(i, 1), gradN(i, 2)});
            result[i] = dN*N.T() + N*dN.T();
        }

        return result;
    }

    inline std::array<math::Matrix, S_SIZE> eps_mat_lambda_2(const gp_Dir& n, const gp_Pnt& p) const{
        std::array<math::Matrix, S_SIZE> eps;

        const auto grad = this->square_lambda_deriv(p);

        eps[0] = n.X()*grad[0];
        eps[1] = n.Y()*grad[1];
        eps[2] = n.Z()*grad[2];

        eps[3] = n.X()*grad[1] + n.Y()*grad[0];
        eps[4] = n.X()*grad[2] + n.Z()*grad[0];
        eps[5] = n.Y()*grad[2] + n.Z()*grad[1];

        return eps;
    }

    inline math::Vector eps_vec(const gp_Pnt& p, const math::Vector& le_n, const math::Vector& le_p1, const math::Vector& le_p2, const gp_Dir& n, const gp_Dir& p1, const gp_Dir& p2) const{
        const math::Vector N(this->N_mat_1dof(p));

        const auto eps_n = this->eps_mat_lambda_2(n, p);
        const auto eps_p1 = this->eps_mat_lambda(p1);
        const auto eps_p2 = this->eps_mat_lambda(p2);

        math::Vector eps_n_full(S_SIZE);
        for(size_t i = 0; i < S_SIZE; ++i){
            eps_n_full[i] = le_n.T()*eps_n[i]*le_n;
        }

        return eps_n_full + eps_p1*le_p1 + eps_p2*le_p2;
    }

    inline math::Matrix delta_eps_1(const gp_Pnt& p, const math::Vector& le_n, const gp_Dir& n, const gp_Dir& p1, const gp_Dir& p2) const{
        const math::Vector N(this->N_mat_1dof(p));

        const auto eps_n = this->eps_mat_lambda_2(n, p);
        const auto eps_p1 = this->eps_mat_lambda(p1);
        const auto eps_p2 = this->eps_mat_lambda(p2);

        math::Matrix eps_n_full(S_SIZE, NODES_PER_ELEM);
        for(size_t i = 0; i < S_SIZE; ++i){
            const auto mult1(eps_n[i]*le_n);
            for(size_t j = 0; j < NODES_PER_ELEM; ++j){
                eps_n_full(i, j) = mult1[j];
            }
        }

        math::Matrix eps_f(S_SIZE, 3*NODES_PER_ELEM);

        for(size_t i = 0; i < S_SIZE; ++i){
            for(size_t j = 0; j < NODES_PER_ELEM; ++j){
                eps_f(i,j + 0*NODES_PER_ELEM) = 2*eps_n_full(i,j);
                eps_f(i,j + 1*NODES_PER_ELEM) = eps_p1(i,j);
                eps_f(i,j + 2*NODES_PER_ELEM) = eps_p2(i,j);
            }
        }

        return eps_f;
    }

    inline math::Matrix delta_eps_2(const gp_Pnt& p, const gp_Dir& n, const math::Vector& eps, const math::Matrix& D) const{
        const math::Vector N(this->N_mat_1dof(p));

        const auto eps_n = this->eps_mat_lambda_2(n, p);

        const auto eps_D(eps.T()*D);

        math::Matrix result(3*NODES_PER_ELEM, 3*NODES_PER_ELEM);
        math::Matrix sum(NODES_PER_ELEM, NODES_PER_ELEM);
        for(size_t k = 0; k < S_SIZE; ++k){
            sum += eps_D[k]*eps_n[k];
        }            
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            for(size_t j = 0; j < NODES_PER_ELEM; ++j){
                result(i,j) = 2*sum(i,j);
            } 
        }

        return result;
    }

    inline math::Vector eps_vec_eta(const gp_Pnt& p, const double eta, const math::Vector& le_n, const math::Vector& dle_n, const math::Vector& dle_p1, const math::Vector& dle_p2, const gp_Dir& n, const gp_Dir& p1, const gp_Dir& p2) const{
        (void) eta;
        const math::Vector N(this->N_mat_1dof(p));

        const auto eps_n = this->eps_mat_lambda_2(n, p);
        const auto eps_p1 = this->eps_mat_lambda(p1);
        const auto eps_p2 = this->eps_mat_lambda(p2);

        math::Vector eps_n_full(S_SIZE);
        for(size_t i = 0; i < S_SIZE; ++i){
            eps_n_full[i] = 2*((le_n.T())*eps_n[i]*dle_n);
        }

        return eps_n_full + eps_p1*dle_p1 + eps_p2*dle_p2;
    }

    inline math::Vector eps_D_eta(const gp_Pnt& p, const math::Vector& dle_n, const gp_Dir& n, const math::Vector& eps, const math::Matrix& D) const{
        const math::Vector N(this->N_mat_1dof(p));

        const auto eps_n = this->eps_mat_lambda_2(n, p);

        const auto eps_D(eps.T()*D);

        math::Vector result(3*NODES_PER_ELEM);
        math::Matrix sum(NODES_PER_ELEM, NODES_PER_ELEM);
        for(size_t k = 0; k < S_SIZE; ++k){
            sum += eps_D[k]*eps_n[k];
        }    
        math::Vector sum_mult(sum*dle_n);        
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            result[i] = 2*sum_mult[i];
        }

        return result;
    }

    inline math::Vector fl4_eps_vec(const gp_Pnt& p, const math::Vector& le_n, const gp_Dir& n) const{
        const math::Vector N(this->N_mat_1dof(p));

        const auto eps_n = this->eps_mat_lambda_2(n, p);

        math::Vector eps_n_full(S_SIZE);
        for(size_t i = 0; i < S_SIZE; ++i){
            eps_n_full[i] = le_n.T()*eps_n[i]*le_n;
        }

        return eps_n_full;
    }

    inline math::Matrix fl4_delta_eps_1(const gp_Pnt& p, const math::Vector& le_n, const gp_Dir& n) const{
        const math::Vector N(this->N_mat_1dof(p));

        const auto eps_n = this->eps_mat_lambda_2(n, p);

        math::Matrix eps_n_full(S_SIZE, NODES_PER_ELEM);
        for(size_t i = 0; i < S_SIZE; ++i){
            const auto mult1(eps_n[i]*le_n);
            for(size_t j = 0; j < NODES_PER_ELEM; ++j){
                eps_n_full(i, j) = mult1[j];
            }
        }

        return eps_n_full;
    }

    inline math::Matrix fl4_delta_eps_2(const gp_Pnt& p, const gp_Dir& n, const math::Vector& eps, const math::Matrix& D) const{
        const math::Vector N(this->N_mat_1dof(p));

        const auto eps_n = this->eps_mat_lambda_2(n, p);

        const auto eps_D(eps.T()*D);

        math::Matrix sum(NODES_PER_ELEM, NODES_PER_ELEM);
        for(size_t k = 0; k < S_SIZE; ++k){
            sum += eps_D[k]*eps_n[k];
        }            

        return sum;
    }
};

}

#endif
