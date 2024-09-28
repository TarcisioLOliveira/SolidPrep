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
#include "element/Q4.hpp"
#include "logger.hpp"
#include "utils/gauss_legendre.hpp"

namespace element{

H8::H8(ElementShape s):
    MeshElementCommon3DHex<H8>(s.nodes){
}

std::unique_ptr<BoundaryMeshElementFactory> H8::get_boundary_element_info() {
    logger::log_assert(false, logger::ERROR, "BQ4 ELEMENT TYPE NOT IMPLEMENTED");
    return std::unique_ptr<BoundaryMeshElementFactory>();
    //return std::unique_ptr<BoundaryMeshElementFactory>(new BoundaryMeshElementFactoryImpl<boundary_element::BQ4>());
}
std::unique_ptr<ContactMeshElementFactory> H8::get_contact_element_info() {
    logger::log_assert(false, logger::ERROR, "CQ4 ELEMENT TYPE NOT IMPLEMENTED");
    return std::unique_ptr<ContactMeshElementFactory>();
    //return std::unique_ptr<ContactMeshElementFactory>(new ContactMeshElementFactoryImpl<contact_element::CQ4>());
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

    std::array<double, BOUNDARY_NODES_PER_ELEM> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
        x[i] = points[i].X() - c.X();
        y[i] = points[i].Y() - c.Y();
        z[i] = points[i].Z() - c.Z();
    }
    std::array<double, NODES_PER_ELEM> x2, y2, z2;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        x2[i] = this->nodes[i]->point.X() - c.X();
        y2[i] = this->nodes[i]->point.Y() - c.Y();
        z2[i] = this->nodes[i]->point.Z() - c.Z();
    }

    const CubeSide cs = this->get_cube_side(points);

    Eigen::Matrix<double, NODE_DOF, K_DIM> Nf;
    Nf.fill(0);
    constexpr size_t GN = 2;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto drnorm = this->surface_drnorm(xi->x, eta->x, x, y, z);
            gp_Pnt pi = this->to_surface_point(xi->x, eta->x, cs);
            Nf += (xi->w*eta->w*drnorm)*N_mat_norm(pi.X(), pi.Y(), pi.Z());
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

    std::array<double, BOUNDARY_NODES_PER_ELEM> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
        x[i] = points[i].X() - c.X();
        y[i] = points[i].Y() - c.Y();
        z[i] = points[i].Z() - c.Z();
    }
    std::array<double, NODES_PER_ELEM> x2, y2, z2;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        x2[i] = this->nodes[i]->point.X() - c.X();
        y2[i] = this->nodes[i]->point.Y() - c.Y();
        z2[i] = this->nodes[i]->point.Z() - c.Z();
    }

    const CubeSide cs = this->get_cube_side(points);

    Eigen::Matrix<double, K_DIM, K_DIM> R;
    Eigen::Matrix<double, DIM, DIM> Km = Eigen::Map<const Eigen::Matrix<double, DIM, DIM>>(K.data(), DIM, DIM);
    R.fill(0);
    constexpr size_t GN = 5;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto drnorm = this->surface_drnorm(xi->x, eta->x, x, y, z);
            gp_Pnt pi = this->to_surface_point(xi->x, eta->x, cs);
            const auto NN = N_mat_norm(pi.X(), pi.Y(), pi.Z());
            R += (xi->w*eta->w*drnorm)*NN.transpose()*Km*NN;
        }
    }
    std::vector<double> R_vec(K_DIM*K_DIM);
    std::copy(R.data(), R.data()+K_DIM*K_DIM, R_vec.begin());

    return R_vec;
}

std::vector<double> H8::get_Rf(const std::vector<double>& S, const std::vector<double>& F, const gp_Pnt& C, const double t, const std::vector<gp_Pnt>& points) const{
    (void)t;

    std::array<double, BOUNDARY_NODES_PER_ELEM> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
        x[i] = points[i].X() - c.X();
        y[i] = points[i].Y() - c.Y();
        z[i] = points[i].Z() - c.Z();
    }
    std::array<double, NODES_PER_ELEM> x2, y2, z2;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        x2[i] = this->nodes[i]->point.X() - c.X();
        y2[i] = this->nodes[i]->point.Y() - c.Y();
        z2[i] = this->nodes[i]->point.Z() - c.Z();
    }

    const CubeSide cs = this->get_cube_side(points);

    Eigen::Vector<double, DIM> x_vec;
    Eigen::Vector<double, DIM> Fv{F[0], F[1], F[2]};
    Eigen::Vector<double, K_DIM> Rf;
    Eigen::Matrix<double, DIM, DIM> Sm = Eigen::Map<const Eigen::Matrix<double, DIM, DIM>>(S.data(), DIM, DIM).transpose();
    Rf.fill(0);
    constexpr size_t GN = 4;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto drnorm = this->surface_drnorm(xi->x, eta->x, x, y, z);
            auto Xs = this->surface_to_nat(xi->x, eta->x, x, y, z);
            gp_Pnt pi = this->to_surface_point(xi->x, eta->x, cs);
            x_vec[0] = Xs[0] + c.X() - C.X();
            x_vec[1] = Xs[1] + c.Y() - C.Y();
            x_vec[2] = Xs[2] + c.Z() - C.Z();
            const auto NN = N_mat_norm(pi.X(), pi.Y(), pi.Z());
            Rf += (xi->w*eta->w*drnorm)*NN.transpose()*(Sm*x_vec + Fv);
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

Eigen::VectorXd H8::flow_1dof(const double t, const MeshNode** nodes) const{
    (void)t;

    std::vector<gp_Pnt> points(BOUNDARY_NODES_PER_ELEM);
    std::array<double, BOUNDARY_NODES_PER_ELEM> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
        x[i] = nodes[i]->point.X() - c.X();
        y[i] = nodes[i]->point.Y() - c.Y();
        z[i] = nodes[i]->point.Z() - c.Z();
        points[i] = nodes[i]->point;
    }
    std::array<double, NODES_PER_ELEM> x2, y2, z2;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        x2[i] = this->nodes[i]->point.X() - c.X();
        y2[i] = this->nodes[i]->point.Y() - c.Y();
        z2[i] = this->nodes[i]->point.Z() - c.Z();
    }

    const CubeSide cs = this->get_cube_side(points);

    Eigen::Vector<double, NODES_PER_ELEM> M;
    M.fill(0);
    constexpr size_t GN = 2;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto drnorm = this->surface_drnorm(xi->x, eta->x, x, y, z);
            gp_Pnt pi = this->to_surface_point(xi->x, eta->x, cs);
            M += (xi->w*eta->w*drnorm)*N_mat_1dof(pi.X(), pi.Y(), pi.Z());
        }
    }

    return M;
}

std::vector<double> H8::get_MnMn(const MeshElement* const e2, const std::vector<double>& u_ext, const std::vector<gp_Pnt>& bounds, const gp_Dir n) const{
    const size_t DOF = NODE_DOF;
    const size_t KW = K_DIM;

    std::array<double, BOUNDARY_NODES_PER_ELEM> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
        x[i] = bounds[i].X() - c.X();
        y[i] = bounds[i].Y() - c.Y();
        z[i] = bounds[i].Z() - c.Z();
    }

    std::vector<double> uv1(KW), uv2(KW);
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        for(size_t j = 0; j < DOF; ++j){
            uv1[i*DOF + j] = u_ext[this->nodes[i]->u_pos[j]];
            uv2[i*DOF + j] = u_ext[e2->nodes[i]->u_pos[j]];
        }
    }
    std::vector<double> NN(2*KW, 0);
    std::vector<double> MnMn(4*KW*KW, 0);

    const CubeSide cs = this->get_cube_side(bounds);

    std::vector<double> up1(DIM), up2(DIM);
    constexpr size_t GN = 2*ORDER + 2;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            std::fill(up1.begin(), up1.end(), 0);
            std::fill(up2.begin(), up2.end(), 0);

            const auto drnorm = this->surface_drnorm(xi->x, eta->x, x, y, z);
            const gp_Pnt pi = this->to_surface_point(xi->x, eta->x, cs);
            const auto N1 = this->get_Ni(pi);
            const auto N2 = e2->get_Ni(pi);
            double gp = 0;
            for(size_t j = 0; j < DIM; ++j){
                for(size_t i = 0; i < KW; ++i){
                    up1[j] += N1[j*KW + i]*uv1[i];
                    up2[j] += N2[j*KW + i]*uv2[i];
                }
                gp += (up2[j] - up1[j])*n.Coord(1+j);
            }
            if(gp < 1e-7){
                std::fill(NN.begin(), NN.end(), 0);
                for(size_t i = 0; i < KW; ++i){
                    for(size_t j = 0; j < DIM; ++j){
                        NN[i] += N1[i + j*KW]*n.Coord(1+j);
                        NN[i + KW] += N2[i + j*KW]*n.Coord(1+j);
                    }
                }
                for(size_t i = 0; i < 2*KW; ++i){
                    for(size_t j = 0; j < 2*KW; ++j){
                        MnMn[2*KW*i + j] += (xi->w*eta->w*drnorm)*NN[i]*NN[j];
                    }
                }
            }
        }
    }
    for(size_t i = 0; i < KW; ++i){
        for(size_t j = KW; j < 2*KW; ++j){
            MnMn[2*KW*i + j] *= -1.0;
            MnMn[2*KW*j + i] *= -1.0;
        }
    }
    const gp_Vec v1(bounds[0], bounds[1]);
    const gp_Vec v2(bounds[0], bounds[2]);
    const double A = 0.5*v1.Crossed(v2).Magnitude();
    cblas_dscal(MnMn.size(), A, MnMn.data(), 1);

    return MnMn;
}

H8::CubeSide H8::get_cube_side(const std::vector<gp_Pnt>& points) const{
    // As I was unable to make this work using natural coordinates, this method
    // has additional code to find which face of the hex the loads are being
    // applied onto.
    //
    // This is probably (hopefully) not the best solution, but it works for now.

    int xp[8] = {-1, 1, 1, -1, -1, 1, 1, -1};
    int yp[8] = {1, 1, -1, -1, 1, 1, -1, -1};
    int zp[8] = {-1, -1, -1, -1, 1, 1, 1, 1};

    int idx[4];
    int idx_i = 0;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
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
    logger::log_assert(xt*yt*zt == 0 && xt*yt == 0 && xt*zt == 0 && yt*zt == 0,
                       logger::ERROR,
                       "wrong assumptions getting cube side: {} {} {}", xt, yt, zt);

    if(xt == 1){
        return CubeSide::X_MAX;
    } else if(xt == -1){
        return CubeSide::X_MIN;
    }
    if(yt == 1){
        return CubeSide::Y_MAX;
    } else if(yt == -1){
        return CubeSide::Y_MIN;
    }
    if(zt == 1){
        return CubeSide::Z_MAX;
    } else if(zt == -1){
        return CubeSide::Z_MIN;
    }

    return CubeSide::UNKNOWN;
}
gp_Pnt H8::to_surface_point(const double xi, const double eta, const H8::CubeSide side) const{
    switch(side){
        case CubeSide::X_MAX:
            return {1, xi, eta};
        case CubeSide::X_MIN:
            return {-1, xi, eta};
        case CubeSide::Y_MAX:
            return {xi, 1, eta};
        case CubeSide::Y_MIN:
            return {xi, -1, eta};
        case CubeSide::Z_MAX:
            return {xi, eta, 1};
        case CubeSide::Z_MIN:
            return {xi, eta, -1};
        case CubeSide::UNKNOWN:
            logger::log_assert(false, logger::ERROR,  "unknown cube side for H8");
    }
    return {0,0,0};
}

std::vector<double> H8::get_Ni(const gp_Pnt& p) const{
    const auto Nm = this->N_mat_norm(p.X(), p.Y(), p.Z());
    std::vector<double> Nv(NODE_DOF*K_DIM, 0);
    for(size_t i = 0; i < NODE_DOF; ++i){
        for(size_t j = 0; j < K_DIM; ++j){
            Nv[i*K_DIM + j] = Nm(i, j);
        }
    }
    return Nv;
}

}
