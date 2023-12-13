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
#include "boundary_element/BTRI6.hpp"
#include "logger.hpp"
#include "utils/gauss_legendre.hpp"
#include "utils/boundary_nullifier.hpp"

namespace boundary_element{

BTRI6::BTRI6(ElementShape s, const MeshElement* const parent):
    BoundaryMeshElement(s.nodes, parent){
    const size_t N = BTRI6::NODES_PER_ELEM;
    
    std::array<double, N> x, y, z;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
        z[i] = this->nodes[i]->point.Z();
    }
    std::array<double, N*N> M = 
        {1, x[0], y[0], x[0]*y[0], x[0]*x[0], y[0]*y[0],
         1, x[1], y[1], x[1]*y[1], x[1]*x[1], y[1]*y[1],
         1, x[2], y[2], x[2]*y[2], x[2]*x[2], y[2]*y[2],
         1, x[3], y[3], x[3]*y[3], x[3]*x[3], y[3]*y[3],
         1, x[4], y[4], x[4]*y[4], x[4]*x[4], y[4]*y[4],
         1, x[5], y[5], x[5]*y[5], x[5]*x[5], y[5]*y[5]};

    std::array<int, N> ipiv;

    // M*C = I -> C=M^-1
    // C = {a[0], a[1], a[2],
    //      b[0], b[1], b[2],
    //      c[0], c[1], c[2]}
    int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating LU in BTRI6.", info);
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating computing inverse from LU in BTRI6.", info);

    for(size_t i = 0; i < N; ++i){
        this->a[i] = M[i];
        this->b[i] = M[i+N];
        this->c[i] = M[i+2*N];
        this->d[i] = M[i+3*N];
        this->e[i] = M[i+4*N];
        this->f[i] = M[i+5*N];
    }

    gp_Vec v1(this->nodes[0]->point, this->nodes[1]->point);
    gp_Vec v2(this->nodes[0]->point, this->nodes[2]->point);

    this->delta = 0.5*(v1.Crossed(v2)).Magnitude();
}

Eigen::MatrixXd BTRI6::diffusion_1dof(const Eigen::MatrixXd& A) const{
}
Eigen::MatrixXd BTRI6::advection_1dof(const Eigen::VectorXd& v) const{
}
Eigen::MatrixXd BTRI6::absorption_1dof() const{
}
Eigen::VectorXd BTRI6::source_1dof() const{
    const auto& gsi = utils::GaussLegendreTri<2>::get();
    Eigen::Vector<double, 6> result{0,0,0,0,0,0};
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto NN = this->N_mat_1dof(p);
        result += it->w*NN;
    }
    return result*delta;
}
Eigen::VectorXd BTRI6::grad_1dof_upos(const gp_Pnt& p, const std::vector<double>& phi) const{
    (void)p;
    Eigen::Vector<double, 6> phiv{0, 0, 0, 0, 0, 0};
    for(size_t i = 0; i < 6; ++i){
        const auto p = this->nodes[i]->u_pos[0];
        if(p > -1){
            phiv[i] = phi[p];
        }
    }
    return this->dN_mat_1dof(p)*phiv;
};
Eigen::VectorXd BTRI6::grad_1dof_id(const gp_Pnt& p, const std::vector<double>& phi) const{
    (void)p;
    Eigen::Vector<double, 6> phiv{0, 0, 0, 0, 0, 0};
    for(size_t i = 0; i < 6; ++i){
        const auto p = this->nodes[i]->id;
        phiv[i] = phi[p];
    }
    return this->dN_mat_1dof(p)*phiv;
};

Eigen::MatrixXd BTRI6::int_grad_1dof() const{
    //return delta*this->dN_mat_1dof();
}

Eigen::VectorXd BTRI6::source_1dof(const Eigen::Vector<double, 3>& v) const{
}

std::array<Eigen::VectorXd, 2> BTRI6::source_1dof(const double S13, const double Gxy, const double S12, const double Gxz, const gp_Pnt& center) const{
    const auto& gsi = utils::GaussLegendreTri<4>::get();
    std::array<Eigen::VectorXd, 2>  result{
        Eigen::Vector<double, 6>{0,0,0,0,0,0}, 
        Eigen::Vector<double, 6>{0,0,0,0,0,0}};
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        // CENTER
        Eigen::Vector<double, 3> X{p.X() - center.X(), p.Y() - center.Y(), p.Z() - center.Z()};
        const auto NN = this->N_mat_1dof(p);
        result[0] += it->w*NN*(-2*S13*X[1]);
        result[1] += it->w*NN*(-2*S12*X[2]);
    }
    result[0] *= delta;
    result[1] *= delta;

    return result;
}

Eigen::VectorXd BTRI6::source_1dof(double dcurv_v, double dcurv_w, const gp_Pnt& center) const{
    const auto& gsi = utils::GaussLegendreTri<4>::get();
    Eigen::Vector<double, 6> result{0,0,0,0,0,0};
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        // CENTER
        Eigen::Vector<double, 3> X{p.X() - center.X(), p.Y() - center.Y(), p.Z() - center.Z()};
        const auto NN = this->N_mat_1dof(p);
        result += it->w*NN*(dcurv_v*X[1] + dcurv_w*X[2]);
    }
    return result*delta;
}

Eigen::VectorXd BTRI6::flow_1dof(double dcurv_v, double dcurv_w, const gp_Pnt& center, const gp_Dir& n, const std::vector<gp_Pnt>& edges) const{
    const auto& gsi = utils::GaussLegendre<5>::get();
    Eigen::Vector<double, 6> result{0,0,0,0,0,0};
    gp_Vec dist(edges[0], edges[1]);
    const double rnorm = 0.5*dist.Magnitude();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const double s = (it->x + 1)/2;
        const gp_Pnt p = edges[0].Translated(s*dist);
        // CENTER
        Eigen::Vector<double, 3> X{p.X() - center.X(), p.Y() - center.Y(), p.Z() - center.Z()};
        const auto NN = this->N_mat_1dof(p);
        result += it->w*NN*(-dcurv_w*X[2]*X[2]*n.Y() + dcurv_v*X[1]*X[1]*n.Z())/2;
    }
    return result*rnorm;
}

Eigen::VectorXd BTRI6::source_grad_1dof(const Eigen::VectorXd& v) const{
}

Eigen::MatrixXd BTRI6::L4(const Eigen::MatrixXd& B) const{
    const auto ddN = this->ddN_mat_1dof();

    return delta*ddN.transpose()*B*ddN;
}
Eigen::MatrixXd BTRI6::L3(const Eigen::MatrixXd& B) const{
    const auto ddN = this->ddN_mat_1dof();
    const auto dN = this->dN_mat_1dof(this->GS_point(1.0/3.0, 1.0/3.0, 1.0/3.0));

    return delta*ddN.transpose()*B*dN;
}
Eigen::MatrixXd BTRI6::L2(const Eigen::MatrixXd& B) const{
    const auto& gsi = utils::GaussLegendreTri<2>::get();
    Eigen::Matrix<double, 6, 6> result;
    result.fill(0);
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto dN = this->dN_mat_1dof(p);
        result += it->w*dN.transpose()*B*dN;
    }
    return result*delta;
}
}
