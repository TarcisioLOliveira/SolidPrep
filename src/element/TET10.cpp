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

#include <Eigen/Core>
#include <Eigen/Dense>
#include "element/TET10.hpp"
#include "boundary_element/BTRI6.hpp"
#include "cblas.h"
#include "logger.hpp"
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
#include <vector>
#include "utils/gauss_legendre.hpp"

namespace element{

TET10::TET10(ElementShape s):
    MeshElementCommon3DTet<TET10>(s.nodes){

    this->get_coeffs();    
}

std::unique_ptr<BoundaryMeshElementFactory> TET10::get_boundary_element_info() {
    return std::unique_ptr<BoundaryMeshElementFactory>(new BoundaryMeshElementFactoryImpl<boundary_element::BTRI6>());
}

void TET10::get_coeffs(){
    constexpr size_t N = 4;
    std::array<double, N> x, y, z;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
        z[i] = this->nodes[i]->point.Z();
    }

    std::array<double, N*N> M = 
        {1, x[0], y[0], z[0],
         1, x[1], y[1], z[1],
         1, x[2], y[2], z[2],
         1, x[3], y[3], z[3]};

    std::array<int, N> ipiv;

    // M*C = I -> C=M^-1
    // C = {a[0], a[1], a[2], a[3],
    //      b[0], b[1], b[2], b[3],
    //      c[0], c[1], c[2], c[3],
    //      d[0], d[1], d[2], d[3]}
    int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating LU in TET10.", info);
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating computing inverse from LU in TET10.", info);

    this->V = this->get_volume(1.0);
    for(auto& m:M){
        m *= 6*V;
    }

    std::copy(M.begin(), M.begin()+4, a);
    std::copy(M.begin()+4, M.begin()+8, b);
    std::copy(M.begin()+8, M.begin()+12, c);
    std::copy(M.begin()+12, M.begin()+16, d);
}

std::vector<double> TET10::get_k(const std::vector<double>& D, const double t) const{
    (void)t;
    
    Eigen::Matrix<double, K_DIM, K_DIM> K;
    K.fill(0);
    Eigen::Matrix<double, S_SIZE, S_SIZE> Dm = Eigen::Map<const Eigen::Matrix<double, S_SIZE, S_SIZE>>(D.data(), S_SIZE, S_SIZE);

    const auto& gli = utils::GaussLegendreTet<2>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GL_point(it->a, it->b, it->c, it->d);
        const auto B(this->B_mat(p));
        K += it->w*B.transpose()*Dm*B;
    }
    K *= this->V;

    std::vector<double> k(K_DIM*K_DIM);
    std::copy(K.data(), K.data()+k.size(), k.begin());

    return k;
}
std::vector<double> TET10::get_B(const gp_Pnt& point) const{

    const auto b = this->B_mat(point);
    std::vector<double> B(S_SIZE*K_DIM);
    for(size_t i = 0; i < S_SIZE; ++i){
        for(size_t j = 0; j < K_DIM; ++j){
            B[i*K_DIM + j] = b(i,j);
        }
    }

    return B;
}

std::vector<double> TET10::get_DB(const std::vector<double>& D, const gp_Pnt& point) const{
    Eigen::Matrix<double, S_SIZE, S_SIZE> Dm = Eigen::Map<const Eigen::Matrix<double, 6, 6>>(D.data(), 6, 6);

    const Eigen::Matrix<double, S_SIZE, K_DIM> db = Dm*this->B_mat(point);
    std::vector<double> DB(S_SIZE*K_DIM);
    for(size_t i = 0; i < S_SIZE; ++i){
        for(size_t j = 0; j < K_DIM; ++j){
            DB[i*K_DIM + j] = db(i,j);
        }
    }

    return DB;
}

std::vector<double> TET10::get_Nf(const double t, const std::vector<gp_Pnt>& points) const{
    (void)t;

    const gp_Vec v1(points[0], points[1]);
    const gp_Vec v2(points[0], points[2]);

    const double delta = v1.Crossed(v2).Magnitude()/2;
    Eigen::Matrix<double, DIM, K_DIM> NNN;
    NNN.fill(0);
    const auto& gli = utils::GaussLegendreTri<3>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GL_point_2D(it->a, it->b, it->c, points);
        const Eigen::Matrix<double, DIM, K_DIM> NN(this->N_mat(p));
        NNN += it->w*NN;
    }
    NNN *= delta;

    std::vector<double> Nf(DIM*K_DIM);
    for(size_t i = 0; i < DIM; ++i){
        for(size_t j = 0; j < K_DIM; ++j){
            Nf[j*DIM + i] = NNN(i,j);
        }
    }

    return Nf;
}

std::vector<double> TET10::get_nodal_density_gradient(gp_Pnt p) const{
    const Eigen::Matrix<double, NODE_DOF, NODES_PER_ELEM> dn(this->dN_mat_1dof(p));

    std::vector<double> dN(NODE_DOF*NODES_PER_ELEM);
    for(size_t i = 0; i < NODE_DOF; ++i){
        for(size_t j = 0; j < NODES_PER_ELEM; ++j){
            dN[i*NODES_PER_ELEM + j] = dn(i,j);
        }
    }

    return dN;

}

std::vector<double> TET10::get_R(const std::vector<double>& K, const double t, const std::vector<gp_Pnt>& points) const{
    (void)t;
    Eigen::Matrix<double, K_DIM, K_DIM> R;
    Eigen::Matrix<double, DIM, DIM> Km = Eigen::Map<const Eigen::Matrix<double, DIM, DIM>>(K.data(), DIM, DIM);
    R.fill(0);

    const auto& p = points;
    gp_Vec v1(p[1], p[0]);
    gp_Vec v2(p[2], p[0]);
    const double delta = v1.Crossed(v2).Magnitude()/2;
    const auto& gli = utils::GaussLegendreTri<4>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GL_point_2D(it->a, it->b, it->c, points);
        const auto NN(this->N_mat(p));
        R += it->w*NN.transpose()*Km*NN;
    }
    R *= delta;

    std::vector<double> R_vec(K_DIM*K_DIM);
    std::copy(R.data(), R.data()+K_DIM*K_DIM, R_vec.begin());

    return R_vec;
}

std::vector<double> TET10::get_Rf(const std::vector<double>& S, const std::vector<double>& F, const gp_Pnt& C, const double t, const std::vector<gp_Pnt>& points) const{
    (void)t;

    (void)S;
    (void)F;
    (void)C;
    (void)points;

    logger::log_assert(false, logger::ERROR, "get_Rf() should be unused");

    return std::vector<double>();
}

Eigen::MatrixXd TET10::diffusion_1dof(const double t, const std::vector<double>& A) const{
    (void)t;

    Eigen::Matrix<double, NODES_PER_ELEM, NODES_PER_ELEM> M;
    M.fill(0);
    Eigen::Matrix<double, NODE_DOF, NODE_DOF> Am;
    for(size_t i = 0; i < NODE_DOF; ++i){
        for(size_t j = 0; j < NODE_DOF; ++j){
            Am(i,j) = A[i*NODE_DOF + j];
        }
    }

    const auto& gli = utils::GaussLegendreTet<2>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GL_point(it->a, it->b, it->c, it->d);
        const auto dN(this->dN_mat_1dof(p));
        M += it->w*dN.transpose()*Am*dN;
    }
    M *= this->V;

    return M;
}
Eigen::MatrixXd TET10::advection_1dof(const double t, const std::vector<double>& v) const{
    (void)t;

    Eigen::Matrix<double, NODES_PER_ELEM, NODES_PER_ELEM> M;
    M.fill(0);
    Eigen::Vector<double, NODE_DOF> vv{v[0], v[1], v[2]};

    const auto& gli = utils::GaussLegendreTet<3>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GL_point(it->a, it->b, it->c, it->d);
        const auto dN(this->dN_mat_1dof(p));
        const auto N(this->N_mat_1dof(p));
        M += it->w*N*(vv.transpose()*dN);
    }
    M *= this->V;

    return M;
}
Eigen::MatrixXd TET10::absorption_1dof(const double t) const{
    (void)t;

    Eigen::Matrix<double, NODES_PER_ELEM, NODES_PER_ELEM> M;
    M.fill(0);

    const auto& gli = utils::GaussLegendreTet<4>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GL_point(it->a, it->b, it->c, it->d);
        const auto N(this->N_mat_1dof(p));
        M += it->w*N*N.transpose();;
    }
    M *= this->V;

    return M;
}
Eigen::VectorXd TET10::source_1dof(const double t) const{
    (void)t;

    Eigen::Vector<double, NODES_PER_ELEM> M;
    M.fill(0);

    const auto& gli = utils::GaussLegendreTet<2>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GL_point(it->a, it->b, it->c, it->d);
        const auto N(this->N_mat_1dof(p));
        M += it->w*N;
    }
    M *= this->V;

    return M;
}

Eigen::VectorXd TET10::flow_1dof(const double t, const MeshNode** nodes) const{
    (void)t;

    const gp_Vec v1(nodes[0]->point, nodes[1]->point);
    const gp_Vec v2(nodes[0]->point, nodes[2]->point);

    const double delta = v1.Crossed(v2).Magnitude()/2;
    Eigen::Vector<double, NODES_PER_ELEM> M;
    M.fill(0);

    const std::vector<gp_Pnt> points{nodes[0]->point, nodes[1]->point, nodes[2]->point};

    const auto& gli = utils::GaussLegendreTri<2>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GL_point_2D(it->a, it->b, it->c, points);
        const auto N(this->N_mat_1dof(p));
        M += it->w*N;
    }
    M *= delta;

    return M;
}

}
