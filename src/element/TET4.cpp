/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
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

#include <Eigen/Core>
#include <Eigen/Dense>
#include "element/TET4.hpp"
#include "boundary_element/BTRI3.hpp"
#include "contact_element/CTRI3.hpp"
#include "cblas.h"
#include "logger.hpp"
#include "utils/gauss_legendre.hpp"
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <Eigen/src/Core/Matrix.h>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS.hxx>
#include <IntTools_EdgeEdge.hxx>
#include <lapacke.h>

namespace element{

TET4::TET4(ElementShape s):
    MeshElementCommon3DTet<TET4>(s.nodes){

    this->calculate_coefficients();    
}

std::unique_ptr<BoundaryMeshElementFactory> TET4::get_boundary_element_info() {
    return std::unique_ptr<BoundaryMeshElementFactory>(new BoundaryMeshElementFactoryImpl<boundary_element::BTRI3>());
}

std::unique_ptr<ContactMeshElementFactory> TET4::get_contact_element_info() {
    return std::unique_ptr<ContactMeshElementFactory>(new ContactMeshElementFactoryImpl<contact_element::CTRI3>());
}

void TET4::calculate_coefficients(){
    constexpr size_t N = TET4::NODES_PER_ELEM;

    for(size_t i = 0; i < N; ++i){
        this->C[i*N + 0] = 1;
        this->C[i*N + 1] = this->nodes[i]->point.X();
        this->C[i*N + 2] = this->nodes[i]->point.Y();
        this->C[i*N + 3] = this->nodes[i]->point.Z();
    }

    std::array<int, N> ipiv;

    // M*C = I -> C=M^-1
    // C = {a[0], a[1], a[2], a[3],
    //      b[0], b[1], b[2], b[3],
    //      c[0], c[1], c[2], c[3],
    //      d[0], d[1], d[2], d[3]}
    int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, this->C.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating LU in TET4.", info);
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR, N, this->C.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating computing inverse from LU in TET4.", info);

    this->V = this->get_volume(1.0);
}

std::vector<double> TET4::get_k(const std::vector<double>& D, const double t) const{
    (void)t;
    Eigen::Matrix<double, S_SIZE, S_SIZE> Dm = Eigen::Map<const Eigen::Matrix<double, S_SIZE, S_SIZE>>(D.data(), S_SIZE, S_SIZE);

    const auto B = this->B(this->C);

    Eigen::Matrix<double, K_DIM, K_DIM> Km = this->V*(B.transpose()*Dm*B);

    std::vector<double> K(K_DIM*K_DIM);
    for(size_t i = 0; i < K_DIM; ++i){
        for(size_t j = 0; j < K_DIM; ++j){
            K[K_DIM*i + j] = Km(i,j);
        }
    }

    return K;
}
std::vector<double> TET4::get_B(const gp_Pnt& point) const{
    (void)point;
    const double* const a = this->C.data();
    const double* const b = a + NODES_PER_ELEM;
    const double* const c = b + NODES_PER_ELEM;
    const double* const d = c + NODES_PER_ELEM;
    std::vector<double> B{
        b[0], 0, 0, b[1], 0, 0, b[2], 0, 0, b[3], 0, 0,
        0, c[0], 0, 0, c[1], 0, 0, c[2], 0, 0, c[3], 0,
        0, 0, d[0], 0, 0, d[1], 0, 0, d[2], 0, 0, d[3],
        c[0], b[0], 0, c[1], b[1], 0, c[2], b[2], 0, c[3], b[3], 0,
        d[0], 0, b[0], d[1], 0, b[1], d[2], 0, b[2], d[3], 0, b[3],
        0, d[0], c[0], 0, d[1], c[1], 0, d[2], c[2], 0, d[3], c[3]
    }; 
    return B;
}

std::vector<double> TET4::get_DB(const std::vector<double>& D, const gp_Pnt& point) const{
    (void)point;

    const auto B = this->get_B(point);
    
    std::vector<double> DB(S_SIZE*K_DIM);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, S_SIZE, K_DIM, S_SIZE, 1, D.data(), S_SIZE, B.data(), K_DIM, 0, DB.data(), K_DIM);

    return DB;
}

std::vector<double> TET4::get_Nf(const double t, const std::vector<gp_Pnt>& points) const{
    (void)t;

    const gp_Vec v1(points[0], points[1]);
    const gp_Vec v2(points[0], points[2]);

    const double AA = v1.Crossed(v2).Magnitude()/6;
    double A[4] = {0,0,0,0};
    for(size_t i = 0; i < 4; ++i){
        for(size_t j = 0; j < 3; ++j){
            if(points[j].IsEqual(this->nodes[i]->point, Precision::Confusion())){
                A[i] = AA;
                break;
            }
        }
    }

    std::vector<double> Nf{
        A[0],  0,  0,
         0, A[0],  0,
         0,  0, A[0],
        A[1],  0,  0,
         0, A[1],  0,
         0,  0, A[1],
        A[2],  0,  0,
         0, A[2],  0,
         0,  0, A[2],
        A[3],  0,  0,
         0, A[3],  0,
         0,  0, A[3]
    };

    return Nf;
}

std::vector<double> TET4::get_nodal_density_gradient(gp_Pnt p) const{
    (void)p;
    const double* const a = C.data();
    const double* const b = a + NODES_PER_ELEM;
    const double* const c = b + NODES_PER_ELEM;
    const double* const d = c + NODES_PER_ELEM;
    
    return std::vector<double>{b[0], b[1], b[2], b[3],
                               c[0], c[1], c[2], c[3],
                               d[0], d[1], d[2], d[3]};
}

std::vector<double> TET4::get_R(const std::vector<double>& K, const double t, const std::vector<gp_Pnt>& points) const{
    (void)t;
    Eigen::Matrix<double, K_DIM, K_DIM> R;
    Eigen::Matrix<double, DIM, DIM> Km = Eigen::Map<const Eigen::Matrix<double, DIM, DIM>>(K.data(), DIM, DIM);
    R.fill(0);

    const double GL[3][3] = {{0.5, 0.5, 0},
                             {0, 0.5, 0.5},
                             {0.5, 0, 0.5}};

    const auto& p = points;
    gp_Vec v1(p[1], p[0]);
    gp_Vec v2(p[2], p[0]);
    const double drnorm = (v1.Crossed(v2)).Magnitude();
    for(size_t i = 0; i < 3; ++i){
        const double xi = GL[i][0]*p[0].X() + GL[i][1]*p[1].X() + GL[i][2]*p[2].X();
        const double eta = GL[i][0]*p[0].Y() + GL[i][1]*p[1].Y() + GL[i][2]*p[2].Y();
        const double zeta = GL[i][0]*p[0].Z() + GL[i][1]*p[1].Z() + GL[i][2]*p[2].Z();
        const auto NN = N_mat(xi, eta, zeta, this->C);
        R += (drnorm*NN.transpose()*Km*NN)/3.0;
    }
    std::vector<double> R_vec(K_DIM*K_DIM);
    std::copy(R.data(), R.data()+K_DIM*K_DIM, R_vec.begin());

    return R_vec;
}

std::vector<double> TET4::get_Rf(const std::vector<double>& S, const std::vector<double>& F, const gp_Pnt& C, const double t, const std::vector<gp_Pnt>& points) const{
    (void)t;
    Eigen::Vector<double, DIM> x_vec;
    Eigen::Vector<double, DIM> Fv{F[0], F[1], F[2]};
    Eigen::Vector<double, K_DIM> Rf;
    Eigen::Matrix<double, DIM, DIM> Sm = Eigen::Map<const Eigen::Matrix<double, DIM, DIM>>(S.data(), DIM, DIM).transpose();
    Rf.fill(0);

    const double GL[3][3] = {{0.5, 0.5, 0},
                             {0, 0.5, 0.5},
                             {0.5, 0, 0.5}};

    const auto& p = points;
    gp_Vec v1(p[1], p[0]);
    gp_Vec v2(p[2], p[0]);
    const double drnorm = (v1.Crossed(v2)).Magnitude()/2;
    for(size_t i = 0; i < 3; ++i){
        const double xi = GL[i][0]*p[0].X() + GL[i][1]*p[1].X() + GL[i][2]*p[2].X();
        const double eta = GL[i][0]*p[0].Y() + GL[i][1]*p[1].Y() + GL[i][2]*p[2].Y();
        const double zeta = GL[i][0]*p[0].Z() + GL[i][1]*p[1].Z() + GL[i][2]*p[2].Z();
        x_vec[0] = xi - C.X();
        x_vec[1] = eta - C.Y();
        x_vec[2] = zeta - C.Z();
        const auto NN = N_mat(xi, eta, zeta, this->C);
        Rf += (drnorm*NN.transpose()*(Sm*x_vec + Fv))/3.0;
    }
    std::vector<double> Rf_vec(K_DIM);
    std::copy(Rf.data(), Rf.data()+K_DIM, Rf_vec.begin());

    return Rf_vec;
}

Eigen::MatrixXd TET4::diffusion_1dof(const double t, const std::vector<double>& A) const{
    (void)t;
    const auto B = this->dN_mat_1dof(this->C);

    Eigen::Matrix<double, DIM, DIM> Am = Eigen::Map<const Eigen::Matrix<double, DIM, DIM>>(A.data(), DIM, DIM).transpose();

    Eigen::Matrix<double, NODES_PER_ELEM, NODES_PER_ELEM> M = this->V*(B.transpose()*Am*B);

    return M;
}
Eigen::MatrixXd TET4::advection_1dof(const double t, const std::vector<double>& v) const{
    (void)t;
    Eigen::Vector<double, DIM> vv = Eigen::Map<const Eigen::Vector<double, DIM>>(v.data(), DIM).transpose();
    Eigen::Matrix<double, NODES_PER_ELEM, NODES_PER_ELEM> M;
    M.fill(0);

    auto& gl = utils::GaussLegendreTet<1>::get();

    for(auto it = gl.begin(); it < gl.end(); ++it){
        const gp_Pnt p = this->GL_point(it->a, it->b, it->c, it->d);
        const auto B = this->dN_mat_1dof(this->C);
        const auto N = this->N_mat_1dof(p.X(), p.Y(), p.Z(), this->C);
        M += it->w*(B.transpose()*vv*N.transpose());
    }

    return this->V*M;
}
Eigen::MatrixXd TET4::absorption_1dof(const double t) const{
    (void)t;
    Eigen::Matrix<double, NODES_PER_ELEM, NODES_PER_ELEM> M;
    M.fill(0);

    auto& gl = utils::GaussLegendreTet<1>::get();

    for(auto it = gl.begin(); it < gl.end(); ++it){
        const gp_Pnt p = this->GL_point(it->a, it->b, it->c, it->d);
        const auto N = this->N_mat_1dof(p.X(), p.Y(), p.Z(), this->C);
        M += it->w*(N*N.transpose());
    }

    return this->V*M;
}
Eigen::VectorXd TET4::source_1dof(const double t) const{
    (void)t;

    Eigen::Vector<double, 4> M{
        V/4
        ,
        V/4
        ,
        V/4
        ,
        V/4
    };
    return M;
}

Eigen::VectorXd TET4::flow_1dof(const double t, const MeshNode** nodes) const{
    (void)t;

    const gp_Vec v1(nodes[0]->point, nodes[1]->point);
    const gp_Vec v2(nodes[0]->point, nodes[2]->point);

    const double AA = v1.Crossed(v2).Magnitude()/6;
    double A[4] = {0,0,0,0};
    for(size_t i = 0; i < 4; ++i){
        for(size_t j = 0; j < 3; ++j){
            if(nodes[j]->point.IsEqual(this->nodes[i]->point, Precision::Confusion())){
                A[i] = AA;
                break;
            }
        }
    }

    Eigen::Vector<double, 4> M{
        A[0],
        A[1],
        A[2],
        A[3]
    };

    return M;
}

std::vector<double> TET4::get_Ni(const gp_Pnt& p) const{
    const auto Nm = this->N_mat(p.X(), p.Y(), p.Z(), this->C);

    std::vector<double> Ni(DIM*K_DIM);
    for(size_t i = 0; i < DIM; ++i){
        for(size_t j = 0; j < K_DIM; ++j){
            Ni[i*K_DIM + j] = Nm(i,j);
        }
    }

    return Ni;
}

std::vector<double> TET4::get_dk_sh(const std::vector<double>& D, const double t, const size_t n, const size_t dof) const{
    (void)t;
    Eigen::Matrix<double, S_SIZE, S_SIZE> Dm = Eigen::Map<const Eigen::Matrix<double, S_SIZE, S_SIZE>>(D.data(), S_SIZE, S_SIZE);

    const auto C2 = this->get_C_derivative(n, dof);

    const auto B1 = this->B(this->C);
    const auto B2 = this->B(C2);

    Eigen::Matrix<double, K_DIM, K_DIM> Km = this->V*(B2.transpose()*Dm*B1 + B1.transpose()*Dm*B2);

    std::vector<double> K(K_DIM*K_DIM);
    for(size_t i = 0; i < K_DIM; ++i){
        for(size_t j = 0; j < K_DIM; ++j){
            K[K_DIM*i + j] = Km(i,j);
        }
    }

    return K;
}

TET4::CoeffMat TET4::get_C_derivative(const size_t n, const size_t dof) const{
    constexpr size_t N = TET4::NODES_PER_ELEM;
    CoeffMat C2;
    C2.fill(0);
    C2[N*n + (dof + 1)] = 1;

    std::vector<double> Ctmp(N*N, 0);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, C2.data(), N, C.data(), N, 0, Ctmp.data(), N);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, C.data(), N, Ctmp.data(), N, 0, C2.data(), N);
    cblas_dscal(N*N, -1, C2.data(), 1);

    return C2;
}

std::vector<double> TET4::get_dB_sh(const gp_Pnt& p, const size_t n, const size_t dof) const{
    (void) p;
    const auto C2 = this->get_C_derivative(n, dof);

    const double* const a = C2.data();
    const double* const b = a + NODES_PER_ELEM;
    const double* const c = b + NODES_PER_ELEM;
    const double* const d = c + NODES_PER_ELEM;
    std::vector<double> B{
        b[0], 0, 0, b[1], 0, 0, b[2], 0, 0, b[3], 0, 0,
        0, c[0], 0, 0, c[1], 0, 0, c[2], 0, 0, c[3], 0,
        0, 0, d[0], 0, 0, d[1], 0, 0, d[2], 0, 0, d[3],
        c[0], b[0], 0, c[1], b[1], 0, c[2], b[2], 0, c[3], b[3], 0,
        d[0], 0, b[0], d[1], 0, b[1], d[2], 0, b[2], d[3], 0, b[3],
        0, d[0], c[0], 0, d[1], c[1], 0, d[2], c[2], 0, d[3], c[3]
    }; 

    return B;
}

}
