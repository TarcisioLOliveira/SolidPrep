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
    MeshElementCommon2DTri<TRI3>(s.nodes, s.nodes_sorted){}

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

std::vector<double> TRI3::get_Nf(const double t, const std::vector<gp_Pnt>& points) const{
    const size_t N = this->NODES_PER_ELEM;

    std::vector<gp_Pnt> p;
    for(size_t i = 0; i < N; ++i){
        const auto& n = this->nodes[i];
        p.push_back(n->point);
    }

    gp_Mat deltaM(1, p[0].X(), p[0].Y(), 1, p[1].X(), p[1].Y(), 1, p[2].X(), p[2].Y());

    double delta = 0.5*deltaM.Determinant();

    std::vector<double> a, b, c;
    for(size_t i = 0; i < N; ++i){
        size_t j = (i + 1) % 3;
        size_t k = (i + 2) % 3;

        a.push_back(p[j].X()*p[k].Y() - p[k].X()*p[j].Y());
        b.push_back(p[j].Y() - p[k].Y());
        c.push_back(p[k].X() - p[j].X());
    }

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

}
