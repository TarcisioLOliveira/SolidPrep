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

#include "element/TRI3.hpp"
#include "cblas.h"
#include "logger.hpp"
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS.hxx>
#include <IntTools_EdgeEdge.hxx>

namespace element{

TRI3::TRI3(ElementShape s):
    MeshElementCommon2DTri<TRI3>(s.nodes){
    const size_t N = this->NODES_PER_ELEM;
    
    gp_Pnt p[3];
    for(size_t i = 0; i < N; ++i){
        const auto& n = this->nodes[i];
        p[i] = n->point;
    }
    for(size_t i = 0; i < N; ++i){
        size_t j = (i + 1) % 3;
        size_t k = (i + 2) % 3;

        a[i] = p[j].X()*p[k].Y() - p[k].X()*p[j].Y();
        b[i] = p[j].Y() - p[k].Y();
        c[i] = p[k].X() - p[j].X();
    }

    gp_Mat deltaM(1, p[0].X(), p[0].Y(), 1, p[1].X(), p[1].Y(), 1, p[2].X(), p[2].Y());

    this->delta = 0.5*std::abs(deltaM.Determinant());
}

std::vector<double> TRI3::get_k(const std::vector<double>& D, const double t) const{
    const size_t N = this->NODES_PER_ELEM;

    std::vector<double> B(3*2*N, 0);

    std::vector<gp_Pnt> p;
    for(size_t i = 0; i < N; ++i){
        const auto& n = this->nodes[i];
        p.push_back(n->point);
    }

    gp_Mat deltaM(1, p[0].X(), p[0].Y(), 1, p[1].X(), p[1].Y(), 1, p[2].X(), p[2].Y());

    double Delta = 0.5*std::abs(deltaM.Determinant());

    for(size_t i = 0; i < N; ++i){
        B[i*2 + 0*2*N] = b[i]/(2*Delta);
        B[i*2 + 1*2*N] = 0;
        B[i*2 + 2*2*N] = c[i]/(2*Delta);
        B[i*2 + 0*2*N + 1] = 0;
        B[i*2 + 1*2*N + 1] = c[i]/(2*Delta);
        B[i*2 + 2*2*N + 1] = b[i]/(2*Delta);
    }

    std::vector<double> DB(3*2*N, 0);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 2*N, 3, 1, D.data(), 3, B.data(), 2*N, 0, DB.data(), 2*N);

    std::vector<double> K(2*N*2*N, 0);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 2*N, 2*N, 3, 1, B.data(), 2*N, DB.data(), 2*N, 0, K.data(), 2*N);

    cblas_dscal(K.size(), t*Delta, K.data(), 1);

    return K;
}

std::vector<double> TRI3::get_DB(const std::vector<double>& D, const gp_Pnt& point) const{
    (void)point;
    const size_t N = this->NODES_PER_ELEM;

    std::vector<double> B(3*2*N, 0);

    std::vector<gp_Pnt> p;
    for(size_t i = 0; i < N; ++i){
        const auto& n = this->nodes[i];
        p.push_back(n->point);
    }

    gp_Mat deltaM(1, p[0].X(), p[0].Y(), 1, p[1].X(), p[1].Y(), 1, p[2].X(), p[2].Y());

    double Delta = 0.5*deltaM.Determinant();

    std::vector<double> a, b, c;
    for(size_t i = 0; i < N; ++i){
        size_t j = (i + 1) % 3;
        size_t k = (i + 2) % 3;

        a.push_back(p[j].X()*p[k].Y() - p[k].X()*p[j].Y());
        b.push_back(p[j].Y() - p[k].Y());
        c.push_back(p[k].X() - p[j].X());
    }

    for(size_t i = 0; i < N; ++i){
        B[i*2 + 0*2*N] = b[i]/(2*Delta);
        B[i*2 + 1*2*N] = 0;
        B[i*2 + 2*2*N] = c[i]/(2*Delta);
        B[i*2 + 0*2*N + 1] = 0;
        B[i*2 + 1*2*N + 1] = c[i]/(2*Delta);
        B[i*2 + 2*2*N + 1] = b[i]/(2*Delta);
    }

    std::vector<double> DB(3*2*N, 0);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 2*N, 3, 1, D.data(), 3, B.data(), 2*N, 0, DB.data(), 2*N);

    return DB;
}

std::vector<double> TRI3::get_B(const gp_Pnt& point) const{
    (void)point;
    std::vector<double> B{
    b[0]/(2*delta)
    ,
    0
    ,
    b[1]/(2*delta)
    ,
    0
    ,
    b[2]/(2*delta)
    ,
    0
    ,
    0
    ,
    c[0]/(2*delta)
    ,
    0
    ,
    c[1]/(2*delta)
    ,
    0
    ,
    c[2]/(2*delta)
    ,
    c[0]/(2*delta)
    ,
    b[0]/(2*delta)
    ,
    c[1]/(2*delta)
    ,
    b[1]/(2*delta)
    ,
    c[2]/(2*delta)
    ,
    b[2]/(2*delta)
    };
    return B;
}

std::vector<double> TRI3::get_Nf(const double t, const std::vector<gp_Pnt>& points) const{
    const double x[]{points[0].X(), points[1].X()};
    const double y[]{points[0].Y(), points[1].Y()};

    std::vector<double> Nf({
    t*(2*a[0] + b[0]*x[0] + b[0]*x[1] + c[0]*y[0] + c[0]*y[1])*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])/(4*delta)
    ,
    0
    ,
    0
    ,
    t*(2*a[0] + b[0]*x[0] + b[0]*x[1] + c[0]*y[0] + c[0]*y[1])*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])/(4*delta)
    ,
    t*(2*a[1] + b[1]*x[0] + b[1]*x[1] + c[1]*y[0] + c[1]*y[1])*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])/(4*delta)
    ,
    0
    ,
    0
    ,
    t*(2*a[1] + b[1]*x[0] + b[1]*x[1] + c[1]*y[0] + c[1]*y[1])*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])/(4*delta)
    ,
    t*(2*a[2] + b[2]*x[0] + b[2]*x[1] + c[2]*y[0] + c[2]*y[1])*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])/(4*delta)
    ,
    0
    ,
    0
    ,
    t*(2*a[2] + b[2]*x[0] + b[2]*x[1] + c[2]*y[0] + c[2]*y[1])*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])/(4*delta)
    });
    
    return Nf;
}

std::vector<double> TRI3::helmholtz_tensor(const double t, const double r) const{
    std::vector<double> h{
    t*(3*b[0]*b[0]*r*r + 3*c[0]*c[0]*r*r + 2*delta*delta)/(12*delta)
    ,
    t*(9*b[0]*b[1]*r*r + 9*c[0]*c[1]*r*r + 4*delta*delta)/(36*delta)
    ,
    t*(9*b[0]*b[2]*r*r + 9*c[0]*c[2]*r*r + 4*delta*delta)/(36*delta)
    ,
    t*(9*b[0]*b[1]*r*r + 9*c[0]*c[1]*r*r + 4*delta*delta)/(36*delta)
    ,
    t*(3*b[1]*b[1]*r*r + 3*c[1]*c[1]*r*r + 2*delta*delta)/(12*delta)
    ,
    t*(9*b[1]*b[2]*r*r + 9*c[1]*c[2]*r*r + 4*delta*delta)/(36*delta)
    ,
    t*(9*b[0]*b[2]*r*r + 9*c[0]*c[2]*r*r + 4*delta*delta)/(36*delta)
    ,
    t*(9*b[1]*b[2]*r*r + 9*c[1]*c[2]*r*r + 4*delta*delta)/(36*delta)
    ,
    t*(3*b[2]*b[2]*r*r + 3*c[2]*c[2]*r*r + 2*delta*delta)/(12*delta)
    };

    return h;
}

std::vector<double> TRI3::get_phi_radial(const double t, const double beta, const double vp, const std::vector<double>& axis, const std::vector<double>& center, const double rho) const{
    const size_t N = this->NODES_PER_ELEM;

    std::vector<double> x(3);
    std::vector<double> y(3);
    for(size_t i = 0; i < N; ++i){
        const auto& n = this->nodes[i];
        x[i] = n->point.X();
        y[i] = n->point.Y();
    }

    const Eigen::Vector<double, 2> A(axis[0], axis[1]);
    const Eigen::Vector<double, 2> C(center[0], center[1]);

    Eigen::Matrix<double, 3, 3> phi = get_phi_radial_base(x[0]/3 + x[1]/3 + x[2]/3,
                                                  y[0]/3 + y[1]/3 + y[2]/3,
                                                  A, C,
                                                  t, beta, vp, rho);
    phi *= -27/48;
    double z[3]{0.6, 0.2, 0.2};
    for(long i = 0; i < 3; ++i){
        long j = (i + 1) % 3;
        long k = (i + 2) % 3;

        auto phi_tmp = get_phi_radial_base(z[i]*x[0] + z[j]*x[1] + z[k]*x[2],
                                           z[i]*y[0] + z[j]*y[1] + z[k]*y[2],
                                           A, C,
                                           t, beta, vp, rho);

        phi += 25*phi_tmp/48;
    }
    std::vector<double> phi_vec(9,0);
    phi.transposeInPlace();

    std::copy(phi.data(), phi.data()+9, phi_vec.begin());

    return phi_vec;
}

Eigen::Matrix<double, 3, 3> TRI3::get_phi_radial_base(const double x, const double y, const Eigen::Vector<double, 2>& A, const Eigen::Vector<double, 2>& C, const double t, const double beta, const double vp, const double rho) const{
    const Eigen::Vector<double, 2> P(x, y);
    const Eigen::Vector<double, 2> v = (C - P) + (C - P).dot(A)*A;
    const double vv = v.dot(v);
    const double dv = -2 + A.dot(A);

    const auto NN = this->N_mat_1dof(x, y);
    const auto dNN = this->dN_mat_1dof();

    const auto phi = t*((beta*rho + vp*dv)*NN*NN.transpose() + vv*dNN.transpose()*dNN + vp*NN*v.transpose()*dNN);

    return phi;
}

std::vector<double> TRI3::get_phi_grad(const double t, const double beta) const{
    std::vector<double> phi{
    beta*delta*t/6
    ,
    beta*delta*t/9
    ,
    beta*delta*t/9
    ,
    beta*delta*t/9
    ,
    beta*delta*t/6
    ,
    beta*delta*t/9
    ,
    beta*delta*t/9
    ,
    beta*delta*t/9
    ,
    beta*delta*t/6
    };

    return phi;
}

std::vector<double> TRI3::get_phi_unidirectional(const double t, const double beta, const double l, const std::vector<double>& v, const double vn) const{
    std::vector<double> phi{
    t*(-9*b[0]*b[0]*l*l - 9*c[0]*c[0]*l*l + 2*delta*(3*b[0]*l*v[0]*vn - 2*beta*delta + 3*c[0]*l*v[1]*vn))/(36*delta)
    ,
    t*(-9*b[0]*b[1]*l*l - 9*c[0]*c[1]*l*l + 2*delta*(3*b[1]*l*v[0]*vn - 2*beta*delta + 3*c[1]*l*v[1]*vn))/(36*delta)
    ,
    t*(-9*b[0]*b[2]*l*l - 9*c[0]*c[2]*l*l + 2*delta*(3*b[2]*l*v[0]*vn - 2*beta*delta + 3*c[2]*l*v[1]*vn))/(36*delta)
    ,
    t*(-9*b[0]*b[1]*l*l - 9*c[0]*c[1]*l*l + 2*delta*(3*b[0]*l*v[0]*vn - 2*beta*delta + 3*c[0]*l*v[1]*vn))/(36*delta)
    ,
    t*(-9*b[1]*b[1]*l*l - 9*c[1]*c[1]*l*l + 2*delta*(3*b[1]*l*v[0]*vn - 2*beta*delta + 3*c[1]*l*v[1]*vn))/(36*delta)
    ,
    t*(-9*b[1]*b[2]*l*l - 9*c[1]*c[2]*l*l + 2*delta*(3*b[2]*l*v[0]*vn - 2*beta*delta + 3*c[2]*l*v[1]*vn))/(36*delta)
    ,
    t*(-9*b[0]*b[2]*l*l - 9*c[0]*c[2]*l*l + 2*delta*(3*b[0]*l*v[0]*vn - 2*beta*delta + 3*c[0]*l*v[1]*vn))/(36*delta)
    ,
    t*(-9*b[1]*b[2]*l*l - 9*c[1]*c[2]*l*l + 2*delta*(3*b[1]*l*v[0]*vn - 2*beta*delta + 3*c[1]*l*v[1]*vn))/(36*delta)
    ,
    t*(-9*b[2]*b[2]*l*l - 9*c[2]*c[2]*l*l + 2*delta*(3*b[2]*l*v[0]*vn - 2*beta*delta + 3*c[2]*l*v[1]*vn))/(36*delta)
    };

    return phi;
}

std::vector<double> TRI3::helmholtz_vector(const double t) const{
    const double txdelta = this->get_volume(t);

    double Ni = txdelta/NODES_PER_ELEM;

    std::vector<double> NT{Ni, Ni, Ni};

    return NT;
}

std::vector<double> TRI3::get_nodal_density_gradient(gp_Pnt p) const{
    (void)p;
    
    return std::vector<double>{b[0]/(2*delta), b[1]/(2*delta), b[2]/(2*delta),
                               c[0]/(2*delta), c[1]/(2*delta), c[2]/(2*delta)};
}

}
