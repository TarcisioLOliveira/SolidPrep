/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
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

#include <vector>
#include <cblas.h>
#include <lapacke.h>
#include "utils/gauss_legendre.hpp"
#include "element/Q4.hpp"
#include "logger.hpp"

namespace element{

Q4::Q4(ElementShape s):
    MeshElementCommon2DQuad<Q4>(s.nodes){
    
    constexpr size_t N = Q4::NODES_PER_ELEM;
    std::array<double, N> x, y;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
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
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating LU in Q4.", info);
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating computing inverse from LU in Q4.", info);

    for(size_t i = 0; i < N; ++i){
        this->a[i] = M[i];
        this->b[i] = M[i+4];
        this->c[i] = M[i+8];
        this->d[i] = M[i+12];
    }

    this->A = this->get_volume(1.0);
}

std::vector<double> Q4::get_k(const std::vector<double>& D, const double t) const{
    constexpr size_t N = Q4::NODES_PER_ELEM;
    std::array<double, N> x, y;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
    }

    const Eigen::Matrix<double, S_SIZE, S_SIZE> Dm = Eigen::Map<const Eigen::Matrix<double, S_SIZE, S_SIZE>>(D.data(), S_SIZE, S_SIZE);

    Eigen::Matrix<double, K_DIM, K_DIM> k;
    k.fill(0);
    constexpr size_t GN = 3;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto X = this->norm_to_nat(xi->x, eta->x, x, y);
            const auto B = this->B_mat_nat(X[0], X[1]);
            const double detJ = this->J(xi->x, eta->x, x, y).determinant();

            k += (xi->w*eta->w*detJ)*B.transpose()*Dm*B;
        }
    }
    k *= t;

    std::vector<double> K(K_DIM*K_DIM);
    std::copy(k.data(), k.data()+K_DIM*K_DIM, K.begin());

    return K;
}

std::vector<double> Q4::get_Nf(const double t, const std::vector<gp_Pnt>& points) const{
    const double x[]{points[0].X(), points[1].X()};
    const double y[]{points[0].Y(), points[1].Y()};
    const double rnorm = 0.5*points[0].Distance(points[1]);

    Eigen::Matrix<double, NODE_DOF, K_DIM> Nf;
    Nf.fill(0);
    constexpr size_t GN = 1;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        const double s = xi->x;
        const double X = 0.5*(x[0]*(1-s) + x[1]*(1+s));
        const double Y = 0.5*(y[0]*(1-s) + y[1]*(1+s));
        Nf += xi->w*N_mat(X, Y);
    }
    Nf *= t*rnorm;
    std::vector<double> Nf_vec(K_DIM*NODE_DOF);
    std::copy(Nf.data(), Nf.data()+K_DIM*NODE_DOF, Nf_vec.begin());

    return Nf_vec;
}

std::vector<double> Q4::get_R(const std::vector<double>& K, const double t, const std::vector<gp_Pnt>& points) const{
    const double x[]{points[0].X(), points[1].X()};
    const double y[]{points[0].Y(), points[1].Y()};
    const double rnorm = 0.5*points[0].Distance(points[1]);

    Eigen::Matrix<double, K_DIM, K_DIM> R;
    Eigen::Matrix<double, DIM, DIM> Km = Eigen::Map<const Eigen::Matrix<double, DIM, DIM>>(K.data(), DIM, DIM);
    R.fill(0);
    constexpr size_t GN = 3;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        const double s = xi->x;
        const double X = 0.5*(x[0]*(1-s) + x[1]*(1+s));
        const double Y = 0.5*(y[0]*(1-s) + y[1]*(1+s));
        const auto NN = N_mat(X, Y);
        R += xi->w*NN.transpose()*Km*NN;
    }
    R *= t*rnorm;
    std::vector<double> R_vec(K_DIM*K_DIM);
    std::copy(R.data(), R.data()+K_DIM*K_DIM, R_vec.begin());

    return R_vec;
}

Eigen::MatrixXd Q4::diffusion_1dof(const double t, const std::vector<double>& A) const{
    constexpr size_t N = Q4::NODES_PER_ELEM;
    std::array<double, N> x, y;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
    }

    const Eigen::Matrix<double, NODE_DOF, NODE_DOF> Am{{A[0], A[1]},
                                                       {A[3], A[4]}};

    Eigen::Matrix<double, NODES_PER_ELEM, NODES_PER_ELEM> M;
    M.fill(0);
    constexpr size_t GN = 3;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto X = this->norm_to_nat(xi->x, eta->x, x, y);
            const auto B = this->dN_mat_1dof(X[0], X[1]);
            const double detJ = this->J(xi->x, eta->x, x, y).determinant();

            M += (xi->w*eta->w*detJ)*B.transpose()*Am*B;
        }
    }
    M *= t;

    return M;
}
Eigen::MatrixXd Q4::advection_1dof(const double t, const std::vector<double>& v) const{
    constexpr size_t N = Q4::NODES_PER_ELEM;
    std::array<double, N> x, y;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
    }

    const Eigen::Vector<double, NODE_DOF> vv{v[0], v[1]};

    Eigen::Matrix<double, NODES_PER_ELEM, NODES_PER_ELEM> M;
    M.fill(0);
    constexpr size_t GN = 4;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto X = this->norm_to_nat(xi->x, eta->x, x, y);
            const auto B = this->dN_mat_1dof(X[0], X[1]);
            const auto Nv = this->N_mat_1dof(X[0], X[1]);
            const double detJ = this->J(xi->x, eta->x, x, y).determinant();

            M += (xi->w*eta->w*detJ)*B.transpose()*vv*Nv.transpose();
        }
    }
    M *= t;

    return M;
}
Eigen::MatrixXd Q4::absorption_1dof(const double t) const{
    constexpr size_t N = Q4::NODES_PER_ELEM;
    std::array<double, N> x, y;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
    }

    Eigen::Matrix<double, NODES_PER_ELEM, NODES_PER_ELEM> M;
    M.fill(0);
    constexpr size_t GN = 5;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto X = this->norm_to_nat(xi->x, eta->x, x, y);
            const auto Nv = this->N_mat_1dof(X[0], X[1]);
            const double detJ = this->J(xi->x, eta->x, x, y).determinant();

            M += (xi->w*eta->w*detJ)*Nv*Nv.transpose();
        }
    }
    M *= t;

    return M;
}
Eigen::VectorXd Q4::source_1dof(const double t) const{
    constexpr size_t N = Q4::NODES_PER_ELEM;
    std::array<double, N> x, y;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
    }

    Eigen::Vector<double, NODES_PER_ELEM> M;
    M.fill(0);
    constexpr size_t GN = 3;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto X = this->norm_to_nat(xi->x, eta->x, x, y);
            const auto Nv = this->N_mat_1dof(X[0], X[1]);
            const double detJ = this->J(xi->x, eta->x, x, y).determinant();

            M += (xi->w*eta->w*detJ)*Nv;
        }
    }
    M *= t;

    return M;
}

std::vector<double> Q4::get_B(const gp_Pnt& point) const{
    const double x = point.X();
    const double y = point.Y();

    std::vector<double> B{
    b[0] + d[0]*y
    ,
    0
    ,
    b[1] + d[1]*y
    ,
    0
    ,
    b[2] + d[2]*y
    ,
    0
    ,
    b[3] + d[3]*y
    ,
    0
    ,
    0
    ,
    c[0] + d[0]*x
    ,
    0
    ,
    c[1] + d[1]*x
    ,
    0
    ,
    c[2] + d[2]*x
    ,
    0
    ,
    c[3] + d[3]*x
    ,
    c[0] + d[0]*x
    ,
    b[0] + d[0]*y
    ,
    c[1] + d[1]*x
    ,
    b[1] + d[1]*y
    ,
    c[2] + d[2]*x
    ,
    b[2] + d[2]*y
    ,
    c[3] + d[3]*x
    ,
    b[3] + d[3]*y
    };
    return B;
}

std::vector<double> Q4::get_DB(const std::vector<double>& D, const gp_Pnt& point) const{
    const double x = point.X();
    const double y = point.Y();

    std::vector<double> DB{
    D[0]*b[0] + D[0]*d[0]*y + D[2]*c[0] + D[2]*d[0]*x
    ,
    D[1]*c[0] + D[1]*d[0]*x + D[2]*b[0] + D[2]*d[0]*y
    ,
    D[0]*b[1] + D[0]*d[1]*y + D[2]*c[1] + D[2]*d[1]*x
    ,
    D[1]*c[1] + D[1]*d[1]*x + D[2]*b[1] + D[2]*d[1]*y
    ,
    D[0]*b[2] + D[0]*d[2]*y + D[2]*c[2] + D[2]*d[2]*x
    ,
    D[1]*c[2] + D[1]*d[2]*x + D[2]*b[2] + D[2]*d[2]*y
    ,
    D[0]*b[3] + D[0]*d[3]*y + D[2]*c[3] + D[2]*d[3]*x
    ,
    D[1]*c[3] + D[1]*d[3]*x + D[2]*b[3] + D[2]*d[3]*y
    ,
    D[3]*b[0] + D[3]*d[0]*y + D[5]*c[0] + D[5]*d[0]*x
    ,
    D[4]*c[0] + D[4]*d[0]*x + D[5]*b[0] + D[5]*d[0]*y
    ,
    D[3]*b[1] + D[3]*d[1]*y + D[5]*c[1] + D[5]*d[1]*x
    ,
    D[4]*c[1] + D[4]*d[1]*x + D[5]*b[1] + D[5]*d[1]*y
    ,
    D[3]*b[2] + D[3]*d[2]*y + D[5]*c[2] + D[5]*d[2]*x
    ,
    D[4]*c[2] + D[4]*d[2]*x + D[5]*b[2] + D[5]*d[2]*y
    ,
    D[3]*b[3] + D[3]*d[3]*y + D[5]*c[3] + D[5]*d[3]*x
    ,
    D[4]*c[3] + D[4]*d[3]*x + D[5]*b[3] + D[5]*d[3]*y
    ,
    D[6]*b[0] + D[6]*d[0]*y + D[8]*c[0] + D[8]*d[0]*x
    ,
    D[7]*c[0] + D[7]*d[0]*x + D[8]*b[0] + D[8]*d[0]*y
    ,
    D[6]*b[1] + D[6]*d[1]*y + D[8]*c[1] + D[8]*d[1]*x
    ,
    D[7]*c[1] + D[7]*d[1]*x + D[8]*b[1] + D[8]*d[1]*y
    ,
    D[6]*b[2] + D[6]*d[2]*y + D[8]*c[2] + D[8]*d[2]*x
    ,
    D[7]*c[2] + D[7]*d[2]*x + D[8]*b[2] + D[8]*d[2]*y
    ,
    D[6]*b[3] + D[6]*d[3]*y + D[8]*c[3] + D[8]*d[3]*x
    ,
    D[7]*c[3] + D[7]*d[3]*x + D[8]*b[3] + D[8]*d[3]*y
    };

    return DB;
}

std::vector<double> Q4::get_nodal_density_gradient(gp_Pnt p) const{
    (void)p;

    const double x = p.X();
    const double y = p.Y();

    return std::vector<double>{b[0] + d[0]*y, b[1] + d[1]*y, b[2] + d[2]*y, b[3] + d[3]*y,
                               c[0] + d[0]*x, c[1] + d[1]*x, c[2] + d[2]*x, c[3] + d[3]*x};
}

}


