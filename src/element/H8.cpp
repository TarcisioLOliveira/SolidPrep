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

#include <lapacke.h>
#include "element/H8.hpp"
#include "logger.hpp"
#include "utils/gauss_legendre.hpp"

namespace element{

H8::H8(ElementShape s):
    MeshElementCommon3DHex<H8>(s.nodes){
}

std::vector<double> H8::get_k(const std::vector<double>& D, const double t) const{
    (void)t;
    constexpr size_t N = H8::NODES_PER_ELEM;
    std::array<double, N> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X() - c.X();
        y[i] = this->nodes[i]->point.Y() - c.Y();
        z[i] = this->nodes[i]->point.Z() - c.Z();
    }

    const Eigen::Matrix<double, S_SIZE, S_SIZE> Dm = Eigen::Map<const Eigen::Matrix<double, S_SIZE, S_SIZE>>(D.data(), S_SIZE, S_SIZE);

    Eigen::Matrix<double, K_DIM, K_DIM> k;
    k.fill(0);
    constexpr size_t GN = 3;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            for(auto zeta = GL.begin(); zeta < GL.end(); ++zeta){
                const auto B = this->B_mat_norm(xi->x, eta->x, zeta->x, x, y, z);
                const double detJ = this->J(xi->x, eta->x, zeta->x, x, y, z).determinant();

                k += B.transpose()*((xi->w*eta->w*zeta->w*detJ)*Dm)*B;
            }
        }
    }

    std::vector<double> K(K_DIM*K_DIM);
    std::copy(k.data(), k.data()+K_DIM*K_DIM, K.begin());

    return K;
}

std::vector<double> H8::get_nodal_density_gradient(gp_Pnt p) const{
    constexpr size_t N = H8::NODES_PER_ELEM;
    std::array<double, N> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X() - c.X();
        y[i] = this->nodes[i]->point.Y() - c.Y();
        z[i] = this->nodes[i]->point.Z() - c.Z();
    }

    auto M = this->dN_mat_1dof(p.X(), p.Y(), p.Z(), x, y, z);
    std::vector<double> Mv(NODE_DOF*NODES_PER_ELEM);
    for(size_t i = 0; i < NODE_DOF; ++i){
        for(size_t j = 0; j < NODES_PER_ELEM; ++j){
            Mv[i*NODES_PER_ELEM + j] = M(i, j);
        }
    }

    return Mv;
}

std::vector<double> H8::get_DB(const std::vector<double>& D, const gp_Pnt& point) const{
    constexpr size_t N = H8::NODES_PER_ELEM;
    std::array<double, N> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X() - c.X();
        y[i] = this->nodes[i]->point.Y() - c.Y();
        z[i] = this->nodes[i]->point.Z() - c.Z();
    }
    const double px = point.X() - c.X();
    const double py = point.Y() - c.Y();
    const double pz = point.Z() - c.Z();
    const auto B = this->B_mat_norm(px, py, pz, x, y, z);
    const Eigen::Matrix<double, S_SIZE, S_SIZE> Dm = Eigen::Map<const Eigen::Matrix<double, S_SIZE, S_SIZE>>(D.data(), S_SIZE, S_SIZE);
    const Eigen::Matrix<double, S_SIZE, K_DIM> DBv = Dm*B;
    std::vector<double> DB(DBv.size(),0);
    for(size_t i = 0; i < S_SIZE; ++i){
        for(size_t j = 0; j < K_DIM; ++j){
            DB[i*K_DIM + j] = DBv(i, j);
        }
    }

    return DB;
}

std::vector<double> H8::get_B(const gp_Pnt& point) const{
    constexpr size_t N = H8::NODES_PER_ELEM;
    std::array<double, N> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X() - c.X();
        y[i] = this->nodes[i]->point.Y() - c.Y();
        z[i] = this->nodes[i]->point.Z() - c.Z();
    }
    const double px = point.X() - c.X();
    const double py = point.Y() - c.Y();
    const double pz = point.Z() - c.Z();
    const auto Bv = this->B_mat_norm(px, py, pz, x, y, z);
    std::vector<double> B(Bv.size(),0);
    for(size_t i = 0; i < S_SIZE; ++i){
        for(size_t j = 0; j < K_DIM; ++j){
            B[i*K_DIM + j] = Bv(i, j);
        }
    }

    return B;
}

std::vector<double> H8::get_Nf(const double t, const std::vector<gp_Pnt>& points) const{
    (void)t;

    // As I was unable to make this work using natural coordinates, this method
    // has additional code to find which face of the hex the loads are being
    // applied onto.
    //
    // This is probably (hopefully) not the best solution, but it works for now.

    int xp[8] = {-1, 1, 1, -1, -1, 1, 1, -1};
    int yp[8] = {1, 1, -1, -1, 1, 1, -1, -1};
    int zp[8] = {-1, -1, -1, -1, 1, 1, 1, 1};

    std::array<double, BOUNDARY_NODES_PER_ELEM> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
        x[i] = points[i].X() - c.X();
        y[i] = points[i].Y() - c.Y();
        z[i] = points[i].Z() - c.Z();
    }
    std::array<double, NODES_PER_ELEM> x2, y2, z2;
    int idx[4];
    int idx_i = 0;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        x2[i] = this->nodes[i]->point.X() - c.X();
        y2[i] = this->nodes[i]->point.Y() - c.Y();
        z2[i] = this->nodes[i]->point.Z() - c.Z();

        for(size_t j = 0; j < BOUNDARY_NODES_PER_ELEM; ++j){
            if(this->nodes[i]->point.IsEqual(points[j], Precision::Confusion())){
                idx[idx_i] = i;
                ++idx_i;
            }
        }
    }
    int xt = xp[idx[0]];
    int yt = yp[idx[0]]; 
    int zt = zp[idx[0]]; 
    for(size_t i = 1; i < 4; ++i){
        if(xp[idx[i]] != xt){
            xt = 0;
        }
        if(yp[idx[i]] != yt){
            yt = 0;
        }
        if(zp[idx[i]] != zt){
            zt = 0;
        }
    }

    Eigen::Matrix<double, NODE_DOF, K_DIM> Nf;
    Nf.fill(0);
    constexpr size_t GN = 2;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto drnorm = this->surface_drnorm(xi->x, eta->x, x, y, z);
            if(xt != 0){
                Nf += (xi->w*eta->w*drnorm)*N_mat_norm(xt, xi->x, eta->x);
            } else if(yt != 0){
                Nf += (xi->w*eta->w*drnorm)*N_mat_norm(xi->x, yt, eta->x);
            } else if(zt != 0){
                Nf += (xi->w*eta->w*drnorm)*N_mat_norm(xi->x, eta->x, zt);
            }
        }
    }
    std::vector<double> Nf_vec(K_DIM*NODE_DOF);
    for(size_t i = 0; i < NODE_DOF; ++i){
        for(size_t j = 0; j < K_DIM; ++j){
            Nf_vec[j*NODE_DOF + i] = Nf(i, j);
        }
    }

    return Nf_vec;
}

std::vector<double> H8::get_R(const std::vector<double>& K, const double t, const std::vector<gp_Pnt>& points) const{
    (void)t;

    // As I was unable to make this work using natural coordinates, this method
    // has additional code to find which face of the hex the loads are being
    // applied onto.
    //
    // This is probably (hopefully) not the best solution, but it works for now.

    int xp[8] = {-1, 1, 1, -1, -1, 1, 1, -1};
    int yp[8] = {1, 1, -1, -1, 1, 1, -1, -1};
    int zp[8] = {-1, -1, -1, -1, 1, 1, 1, 1};

    std::array<double, BOUNDARY_NODES_PER_ELEM> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
        x[i] = points[i].X() - c.X();
        y[i] = points[i].Y() - c.Y();
        z[i] = points[i].Z() - c.Z();
    }
    std::array<double, NODES_PER_ELEM> x2, y2, z2;
    int idx[4];
    int idx_i = 0;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        x2[i] = this->nodes[i]->point.X() - c.X();
        y2[i] = this->nodes[i]->point.Y() - c.Y();
        z2[i] = this->nodes[i]->point.Z() - c.Z();

        for(size_t j = 0; j < BOUNDARY_NODES_PER_ELEM; ++j){
            if(this->nodes[i]->point.IsEqual(points[j], Precision::Confusion())){
                idx[idx_i] = i;
                ++idx_i;
            }
        }
    }
    int xt = xp[idx[0]];
    int yt = yp[idx[0]]; 
    int zt = zp[idx[0]]; 
    for(size_t i = 1; i < 4; ++i){
        if(xp[idx[i]] != xt){
            xt = 0;
        }
        if(yp[idx[i]] != yt){
            yt = 0;
        }
        if(zp[idx[i]] != zt){
            zt = 0;
        }
    }

    Eigen::Matrix<double, K_DIM, K_DIM> R;
    Eigen::Matrix<double, DIM, DIM> Km = Eigen::Map<const Eigen::Matrix<double, DIM, DIM>>(K.data(), DIM, DIM);
    R.fill(0);
    constexpr size_t GN = 5;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto drnorm = this->surface_drnorm(xi->x, eta->x, x, y, z);
            if(xt != 0){
                const auto NN = N_mat_norm(xt, xi->x, eta->x);
                R += (xi->w*eta->w*drnorm)*NN.transpose()*Km*NN;
            } else if(yt != 0){
                const auto NN = N_mat_norm(xi->x, yt, eta->x);
                R += (xi->w*eta->w*drnorm)*NN.transpose()*Km*NN;
            } else if(zt != 0){
                const auto NN = N_mat_norm(xi->x, eta->x, zt);
                R += (xi->w*eta->w*drnorm)*NN.transpose()*Km*NN;
            }
        }
    }
    std::vector<double> R_vec(K_DIM*K_DIM);
    std::copy(R.data(), R.data()+K_DIM*K_DIM, R_vec.begin());

    return R_vec;
}

std::vector<double> H8::get_Rf(const std::vector<double>& S, const std::vector<double>& F, const gp_Pnt& C, const double t, const std::vector<gp_Pnt>& points) const{
    (void)t;

    // As I was unable to make this work using natural coordinates, this method
    // has additional code to find which face of the hex the loads are being
    // applied onto.
    //
    // This is probably (hopefully) not the best solution, but it works for now.

    int xp[8] = {-1, 1, 1, -1, -1, 1, 1, -1};
    int yp[8] = {1, 1, -1, -1, 1, 1, -1, -1};
    int zp[8] = {-1, -1, -1, -1, 1, 1, 1, 1};

    std::array<double, BOUNDARY_NODES_PER_ELEM> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
        x[i] = points[i].X() - c.X();
        y[i] = points[i].Y() - c.Y();
        z[i] = points[i].Z() - c.Z();
    }
    std::array<double, NODES_PER_ELEM> x2, y2, z2;
    int idx[4];
    int idx_i = 0;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        x2[i] = this->nodes[i]->point.X() - c.X();
        y2[i] = this->nodes[i]->point.Y() - c.Y();
        z2[i] = this->nodes[i]->point.Z() - c.Z();

        for(size_t j = 0; j < BOUNDARY_NODES_PER_ELEM; ++j){
            if(this->nodes[i]->point.IsEqual(points[j], Precision::Confusion())){
                idx[idx_i] = i;
                ++idx_i;
            }
        }
    }
    int xt = xp[idx[0]];
    int yt = yp[idx[0]]; 
    int zt = zp[idx[0]]; 
    for(size_t i = 1; i < 4; ++i){
        if(xp[idx[i]] != xt){
            xt = 0;
        }
        if(yp[idx[i]] != yt){
            yt = 0;
        }
        if(zp[idx[i]] != zt){
            zt = 0;
        }
    }

    Eigen::Vector<double, DIM> x_vec;
    Eigen::Vector<double, DIM> Fv{F[0], F[1], F[2]};
    Eigen::Vector<double, K_DIM> Rf;
    Eigen::Matrix<double, DIM, DIM> Sm = Eigen::Map<const Eigen::Matrix<double, DIM, DIM>>(S.data(), DIM, DIM);
    Rf.fill(0);
    constexpr size_t GN = 4;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto drnorm = this->surface_drnorm(xi->x, eta->x, x, y, z);
            if(xt != 0){
                x_vec[0] = xt - C.X();
                x_vec[1] = xi->x - C.Y();
                x_vec[2] = eta->x - C.Z();
                const auto NN = N_mat_norm(xt, xi->x, eta->x);
                Rf += (xi->w*eta->w*drnorm)*NN.transpose()*(Sm*x_vec + Fv);
            } else if(yt != 0){
                x_vec[0] = xi->x - C.X();
                x_vec[1] = yt - C.Y();
                x_vec[2] = eta->x - C.Z();
                const auto NN = N_mat_norm(xi->x, yt, eta->x);
                Rf += (xi->w*eta->w*drnorm)*NN.transpose()*(Sm*x_vec + Fv);
            } else if(zt != 0){
                x_vec[0] = xi->x - C.X();
                x_vec[1] = eta->x - C.Y();
                x_vec[2] = xt - C.Z();
                const auto NN = N_mat_norm(xi->x, eta->x, zt);
                Rf += (xi->w*eta->w*drnorm)*NN.transpose()*(Sm*x_vec + Fv);
            }
        }
    }
    std::vector<double> Rf_vec(K_DIM);
    std::copy(Rf.data(), Rf.data()+K_DIM, Rf_vec.begin());

    return Rf_vec;
}

Eigen::MatrixXd H8::diffusion_1dof(const double t, const std::vector<double>& A) const{
    (void)t;
    constexpr size_t N = H8::NODES_PER_ELEM;
    std::array<double, N> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X() - c.X();
        y[i] = this->nodes[i]->point.Y() - c.Y();
        z[i] = this->nodes[i]->point.Z() - c.Z();
    }

    const Eigen::Matrix<double, NODE_DOF, NODE_DOF> Am{{A[0], A[1], A[2]},
                                                       {A[3], A[4], A[5]},
                                                       {A[6], A[7], A[8]}};
    Eigen::Matrix<double, NODES_PER_ELEM, NODES_PER_ELEM> M;
    M.fill(0);
    constexpr size_t GN = 4;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            for(auto zeta = GL.begin(); zeta < GL.end(); ++zeta){
                const auto B = this->dN_mat_1dof(xi->x, eta->x, zeta->x, x, y, z);
                const double detJ = this->J(xi->x, eta->x, zeta->x, x, y, z).determinant();

                M += B.transpose()*((xi->w*eta->w*zeta->w*detJ)*Am)*B;
            }
        }
    }

    return M;
}
Eigen::MatrixXd H8::advection_1dof(const double t, const std::vector<double>& v) const{
    (void)t;
    constexpr size_t N = H8::NODES_PER_ELEM;
    std::array<double, N> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X() - c.X();
        y[i] = this->nodes[i]->point.Y() - c.Y();
        z[i] = this->nodes[i]->point.Z() - c.Z();
    }

    const Eigen::Vector<double, NODE_DOF> vv{v[0], v[1], v[2]};

    Eigen::Matrix<double, NODES_PER_ELEM, NODES_PER_ELEM> M;
    M.fill(0);
    constexpr size_t GN = 5;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            for(auto zeta = GL.begin(); zeta < GL.end(); ++zeta){
                const auto B = this->dN_mat_1dof(xi->x, eta->x, zeta->x, x, y, z);
                const auto Nv = this->N_mat_1dof(xi->x, eta->x, zeta->x);
                const double detJ = this->J(xi->x, eta->x, zeta->x, x, y, z).determinant();

                M += B.transpose()*((xi->w*eta->w*zeta->w*detJ)*vv)*Nv.transpose();
            }
        }
    }

    return M;

}
Eigen::MatrixXd H8::absorption_1dof(const double t) const{
    (void)t;
    constexpr size_t N = H8::NODES_PER_ELEM;
    std::array<double, N> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X() - c.X();
        y[i] = this->nodes[i]->point.Y() - c.Y();
        z[i] = this->nodes[i]->point.Z() - c.Z();
    }

    Eigen::Matrix<double, NODES_PER_ELEM, NODES_PER_ELEM> M;
    M.fill(0);
    constexpr size_t GN = 6;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            for(auto zeta = GL.begin(); zeta < GL.end(); ++zeta){
                const auto Nv = this->N_mat_1dof(xi->x, eta->x, zeta->x);
                const double detJ = this->J(xi->x, eta->x, zeta->x, x, y, z).determinant();

                M += Nv*((xi->w*eta->w*zeta->w*detJ))*Nv.transpose();
            }
        }
    }

    return M;
}
Eigen::VectorXd H8::source_1dof(const double t) const{
    (void)t;
    constexpr size_t N = H8::NODES_PER_ELEM;
    std::array<double, N> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X() - c.X();
        y[i] = this->nodes[i]->point.Y() - c.Y();
        z[i] = this->nodes[i]->point.Z() - c.Z();
    }

    Eigen::Vector<double, NODES_PER_ELEM> M;
    M.fill(0);
    constexpr size_t GN = 2;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            for(auto zeta = GL.begin(); zeta < GL.end(); ++zeta){
                const auto Nv = this->N_mat_1dof(xi->x, eta->x, zeta->x);
                const double detJ = this->J(xi->x, eta->x, zeta->x, x, y, z).determinant();

                M += (xi->w*eta->w*zeta->w*detJ)*Nv;
            }
        }
    }

    return M;
}

}
