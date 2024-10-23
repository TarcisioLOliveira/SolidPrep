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

#include <lapacke.h>
#include <Eigen/Dense>
#include "boundary_element/BQ4.hpp"
#include "element/H8.hpp"
#include "logger.hpp"
#include "utils/gauss_legendre.hpp"
#include "element_factory.hpp"

namespace boundary_element{

BQ4::BQ4(ElementShape s, const MeshElement* const parent):
    BoundaryMeshElement(s.nodes, parent){
    const size_t N = BQ4::NODES_PER_ELEM;
   
    std::array<double, N> x, y, z;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
        z[i] = this->nodes[i]->point.Z();
    }

    std::array<double, N*N> M = 
        {1, x[0], y[0], x[0]*y[0],
         1, x[1], y[1], x[1]*y[1],
         1, x[2], y[2], x[2]*y[2],
         1, x[3], y[3], x[3]*y[3]};

    std::array<int, N> ipiv;

    // M*C = I -> C=M^-1
    // C = {a[0], a[1], a[2], a[3],
    //      b[0], b[1], b[2], b[3],
    //      c[0], c[1], c[2], c[3],
    //      d[0], d[1], d[2], d[3]}
    int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating LU in BQ4.", info);
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating computing inverse from LU in BQ4.", info);

    for(size_t i = 0; i < N; ++i){
        this->a[i] = M[i];
        this->b[i] = M[i+N];
        this->c[i] = M[i+2*N];
        this->d[i] = M[i+3*N];
    }

    gp_Vec AC(this->nodes[0]->point, this->nodes[2]->point);
    gp_Vec BD(this->nodes[1]->point, this->nodes[3]->point);

    this->A = 0.5*(AC.Crossed(BD)).Magnitude();
}

std::vector<double> BQ4::get_K_ext(const Eigen::MatrixXd& D, const gp_Pnt& center) const{
    Eigen::Matrix<double, K_DIM, K_DIM + 6> Km;
    Km.fill(0);
    const auto& gsi = utils::GaussLegendre<4>::get();
    for(auto xi = gsi.begin(); xi < gsi.end(); ++xi){
        for(auto eta = gsi.begin(); eta < gsi.end(); ++eta){
            const gp_Pnt pi = this->norm_to_nat(xi->x, eta->x);
            const auto eps = this->get_eps(pi, center);
            const auto deps = this->get_deps(pi);
            const double detJ = std::abs(this->J(xi->x, eta->x).determinant());
            Km += (xi->w*eta->w*detJ)*(deps.transpose()*D*eps);
        }
    }
    std::vector<double> K(K_DIM*(K_DIM+6));
    for(size_t i = 0; i < K_DIM; ++i){
        for(size_t j = 0; j < K_DIM + 6; ++j){
            K[i*(K_DIM + 6) + j] = Km(i,j);
        }
    }

    return K;
}
std::vector<double> BQ4::get_normal_stresses(const Eigen::MatrixXd& D, const std::vector<double>& u, const gp_Pnt& p, const gp_Pnt& center) const{
    const auto eps = this->get_eps(p, center);
    std::vector<double> E(S_SIZE, 0);
    std::vector<double> S(3, 0);
    std::vector<double> vals(K_DIM + 6, 0);
    const size_t offset = u.size() - 6;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        size_t pos = this->nodes[i]->id;
        for(size_t j = 0; j < NODE_DOF; ++j){
            vals[i*NODE_DOF + j] = u[NODE_DOF*pos + j];
        }
    }
    for(size_t i = 0; i < 6; ++i){
        vals[K_DIM + i] = u[offset + i];
    }

    for(size_t i = 0; i < S_SIZE; ++i){
        for(size_t j = 0; j < K_DIM + 6; ++j){
            E[i] += eps(i, j)*vals[j];
        }
    }
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < S_SIZE; ++j){
            S[i] += D((i + 2), j)*E[j];
        }
    }

    return S;
}
std::vector<double> BQ4::get_stress_integrals(const Eigen::MatrixXd& D, const gp_Pnt& center) const{
    // M_y M_x V_z M_z V_y V_x
    std::vector<double> M(6*(K_DIM + 6), 0);
    Eigen::Matrix<double, S_SIZE, K_DIM + 6> E;
    const auto& gsi = utils::GaussLegendre<4>::get();
    for(auto xi = gsi.begin(); xi < gsi.end(); ++xi){
        for(auto eta = gsi.begin(); eta < gsi.end(); ++eta){
            const gp_Pnt pi = this->norm_to_nat(xi->x, eta->x);
            const double dx = pi.X() - center.X();
            const double dy = pi.Y() - center.Y();
            const auto eps = this->get_eps(pi, center);
            const double detJ = std::abs(this->J(xi->x, eta->x).determinant());
            E = D*eps;
            for(size_t i = 0; i < K_DIM + 6; ++i){
                M[0*(K_DIM + 6) + i] += (xi->w*eta->w*detJ)*E(2, i)*dx;
                M[1*(K_DIM + 6) + i] += (xi->w*eta->w*detJ)*E(2, i)*dy;
                M[2*(K_DIM + 6) + i] += (xi->w*eta->w*detJ)*E(2, i);
                M[3*(K_DIM + 6) + i] += (xi->w*eta->w*detJ)*(E(3, i)*dx - E(4, i)*dy);
                M[4*(K_DIM + 6) + i] += (xi->w*eta->w*detJ)*E(3, i);
                M[5*(K_DIM + 6) + i] += (xi->w*eta->w*detJ)*E(4, i);
            }
        }
    }

    return M;
}
std::vector<double> BQ4::get_equilibrium_partial(const Eigen::MatrixXd& D, const gp_Pnt& center, const std::vector<size_t>& stresses) const{
    const size_t KW = this->K_DIM; // workaround that's necessary for some reason
    Eigen::MatrixXd Km(KW, K_DIM + 6);
    Km.fill(0);
    const auto& gsi = utils::GaussLegendre<4>::get();
    for(auto xi = gsi.begin(); xi < gsi.end(); ++xi){
        for(auto eta = gsi.begin(); eta < gsi.end(); ++eta){
            const gp_Pnt pi = this->norm_to_nat(xi->x, eta->x);
            const auto eps = this->get_eps(pi, center);
            const auto deps = this->get_deps(pi);
            const Eigen::MatrixXd Ktmp = D*eps;
            const double detJ = std::abs(this->J(xi->x, eta->x).determinant());
            for(size_t i = 0; i < K_DIM; ++i){
                for(size_t j = 0; j < K_DIM + 6; ++j){
                    for(size_t k = 0; k < stresses.size(); ++k){
                        Km(i,j) += (xi->w*eta->w*detJ)*(deps(stresses[k], i)*Ktmp(stresses[k], j));
                    }
                }
            }
        }
    }
    std::vector<double> K(K_DIM*(K_DIM+6));
    for(size_t i = 0; i < K_DIM; ++i){
        for(size_t j = 0; j < K_DIM + 6; ++j){
            K[i*(K_DIM + 6) + j] = Km(i,j);
        }
    }

    return K;
}
std::vector<double> BQ4::get_dz_vector(const Eigen::MatrixXd& S, const Eigen::MatrixXd& D, const double Az, const double Bz, const gp_Pnt& center) const{
    const size_t KW = this->K_DIM; // workaround that's necessary for some reason
    std::vector<double> vec(KW, 0);
    const auto& gsi = utils::GaussLegendre<4>::get();
    for(auto xi = gsi.begin(); xi < gsi.end(); ++xi){
        for(auto eta = gsi.begin(); eta < gsi.end(); ++eta){
            const gp_Pnt pi = this->norm_to_nat(xi->x, eta->x);
            const double dx = pi.X() - center.X();
            const double dy = pi.Y() - center.Y();
            const double detJ = std::abs(this->J(xi->x, eta->x).determinant());
            for(size_t i = 0; i < NODES_PER_ELEM; ++i){
                const double N = this->N(pi.X(), pi.Y(), i);
                vec[i*NODE_DOF + 0] += -(xi->w*eta->w*detJ)*N*D(4,2)*(Az*dx + Bz*dy)/(S(2,2)*D(2,2));
                vec[i*NODE_DOF + 1] += -(xi->w*eta->w*detJ)*N*D(3,2)*(Az*dx + Bz*dy)/(S(2,2)*D(2,2));
                vec[i*NODE_DOF + 2] += -(xi->w*eta->w*detJ)*N*(Az*dx + Bz*dy)/S(2,2);
            }
        }
    }

    return vec;
}
std::vector<double> BQ4::get_force_vector(const Eigen::MatrixXd& D, const std::vector<double>& u, const gp_Pnt& center, const Eigen::MatrixXd& rot) const{
    const size_t KW_P = this->parent->get_element_info()->get_k_dimension();
    const bool is_H8 = this->parent->get_element_info()->get_gmsh_element_type() == 5;

    logger::log_assert(is_H8, logger::ERROR, "BQ4 currently assumes parent element is H8, modifications will be necessary to make it more general");

    const element::H8* const parent = static_cast<const element::H8*>(this->parent);

    std::vector<gp_Pnt> points(NODES_PER_ELEM);
    Eigen::Matrix<double, 3, NODES_PER_ELEM> ptransf;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        for(size_t j = 0; j < 3; ++j){
            ptransf(j, i) = this->nodes[i]->point.Coord(1+j);
        }
    }
    ptransf = rot*ptransf;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        for(size_t j = 0; j < 3; ++j){
            points[i].SetCoord(1+j, ptransf(j, i));
        }
    }

    const element::H8::CubeSide cube_side = parent->get_cube_side(points);

    std::vector<double> F(KW_P, 0);

    const auto& gsi = utils::GaussLegendre<5>::get();
    for(auto xi = gsi.begin(); xi < gsi.end(); ++xi){
        for(auto eta = gsi.begin(); eta < gsi.end(); ++eta){
            const gp_Pnt pi = this->norm_to_nat(xi->x, eta->x);
            const Eigen::Vector<double, 3> pv(rot*Eigen::Vector<double, 3>(pi.X(), pi.Y(), pi.Z()));
            const gp_Pnt pr(pv[0], pv[1], pv[2]);
            auto S = this->get_normal_stresses(D, u, pi, center);
            const Eigen::Vector<double, 3> Sv{S[2], S[1], S[0]};
            const Eigen::Vector<double, 3> Sr(rot*Sv);
            const gp_Pnt p_H8 = parent->to_surface_point(xi->x, eta->x, cube_side);
            const auto N = this->parent->get_Ni(p_H8);
            const double detJ = std::abs(this->J(xi->x, eta->x).determinant());
            for(size_t i = 0; i < 3; ++i){
                for(size_t j = 0; j < KW_P; ++j){
                    F[j] += (xi->w*eta->w*detJ)*N[i*KW_P + j]*Sr[i];
                }
            }
        }
    }

    return F;
}

Eigen::MatrixXd BQ4::diffusion_1dof(const Eigen::MatrixXd& M) const{
    Eigen::Matrix<double, NODES_PER_ELEM, NODES_PER_ELEM> result;
    result.fill(0);
    const auto& gsi = utils::GaussLegendre<5>::get();
    for(auto xi = gsi.begin(); xi < gsi.end(); ++xi){
        for(auto eta = gsi.begin(); eta < gsi.end(); ++eta){
            const gp_Pnt p = this->norm_to_nat(xi->x, eta->x);
            const auto dNN = this->dN_mat_1dof(p);
            const double detJ = std::abs(this->J(xi->x, eta->x).determinant());
            result += (xi->w*eta->w*detJ)*dNN.transpose()*M*dNN;
        }
    }
    return result;
}
Eigen::MatrixXd BQ4::advection_1dof(const Eigen::VectorXd& v) const{
    Eigen::Matrix<double, NODES_PER_ELEM, NODES_PER_ELEM> result;
    result.fill(0);
    const auto& gsi = utils::GaussLegendre<6>::get();
    for(auto xi = gsi.begin(); xi < gsi.end(); ++xi){
        for(auto eta = gsi.begin(); eta < gsi.end(); ++eta){
            const gp_Pnt p = this->norm_to_nat(xi->x, eta->x);
            const auto dNN = this->dN_mat_1dof(p);
            const auto  NN = this->N_mat_1dof(p);
            const double detJ = std::abs(this->J(xi->x, eta->x).determinant());
            result += (xi->w*eta->w*detJ)*NN*(v.transpose()*dNN);
        }
    }
    return result;
}
Eigen::MatrixXd BQ4::absorption_1dof() const{
    Eigen::Matrix<double, NODES_PER_ELEM, NODES_PER_ELEM> result;
    result.fill(0);
    const auto& gsi = utils::GaussLegendre<7>::get();
    for(auto xi = gsi.begin(); xi < gsi.end(); ++xi){
        for(auto eta = gsi.begin(); eta < gsi.end(); ++eta){
            const gp_Pnt p = this->norm_to_nat(xi->x, eta->x);
            const auto NN = this->N_mat_1dof(p);
            const double detJ = std::abs(this->J(xi->x, eta->x).determinant());
            result += (xi->w*eta->w*detJ)*NN*NN.transpose();
        }
    }
    return result;
}
Eigen::VectorXd BQ4::source_1dof() const{
    Eigen::Vector<double, NODES_PER_ELEM> result;
    result.fill(0);
    const auto& gsi = utils::GaussLegendre<5>::get();
    for(auto xi = gsi.begin(); xi < gsi.end(); ++xi){
        for(auto eta = gsi.begin(); eta < gsi.end(); ++eta){
            const gp_Pnt p = this->norm_to_nat(xi->x, eta->x);
            const auto NN = this->N_mat_1dof(p);
            const double detJ = std::abs(this->J(xi->x, eta->x).determinant());
            result += (xi->w*eta->w*detJ)*NN;
        }
    }
    return result;
}
Eigen::VectorXd BQ4::flow_1dof(const std::array<const Node*, 2>& nodes) const{
    std::array<double, 2> x{nodes[0]->point.X(), nodes[1]->point.X()};
    std::array<double, 2> y{nodes[0]->point.Y(), nodes[1]->point.Y()};

    const double rnorm = 0.5*nodes[0]->point.Distance(nodes[1]->point);

    Eigen::Vector<double, NODES_PER_ELEM> M;
    M.fill(0);
    const auto& GL = utils::GaussLegendre<2>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        const double s = xi->x;
        const double X = 0.5*(x[0]*(1-s) + x[1]*(1+s));
        const double Y = 0.5*(y[0]*(1-s) + y[1]*(1+s));
        const gp_Pnt p(X, Y, 0);
        M += xi->w*N_mat_1dof(p);
    }
    M *= rnorm;

    return M;
}

}
