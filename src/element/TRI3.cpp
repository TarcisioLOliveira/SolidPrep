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
#include "project_data.hpp"
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS.hxx>
#include <IntTools_EdgeEdge.hxx>

namespace element{

TRI3::TRI3(ElementShape s, ProjectData* data):
    MeshElementCommon2DTri<TRI3>(s.nodes, data->material.get(), data->thickness){}

std::vector<double> TRI3::get_k() const{
    size_t N = this->nodes.size();

    std::vector<double> B(3*2*N, 0);

    std::vector<gp_Pnt> p;
    for(auto n:this->nodes){
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
    auto D = this->mat->stiffness_2D();

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 2*N, 3, 1, D.data(), 3, B.data(), 2*N, 0, DB.data(), 2*N);

    std::vector<double> K(2*N*2*N, 0);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 2*N, 2*N, 3, 1, B.data(), 2*N, DB.data(), 2*N, 0, K.data(), 2*N);

    cblas_dscal(K.size(), this->t*Delta, K.data(), 1);

    return K;
}

std::vector<double> TRI3::get_DB(const gp_Pnt& point) const{
    (void)point;
    size_t N = this->nodes.size();

    std::vector<double> B(3*2*N, 0);

    std::vector<gp_Pnt> p;
    for(auto n:this->nodes){
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
    auto D = this->mat->stiffness_2D();

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 2*N, 3, 1, D.data(), 3, B.data(), 2*N, 0, DB.data(), 2*N);

    return DB;
}

std::vector<double> TRI3::get_Nf(const std::vector<gp_Pnt>& points) const{

    std::vector<gp_Pnt> p;
    for(auto n:this->nodes){
        p.push_back(n->point);
    }

    gp_Mat deltaM(1, p[0].X(), p[0].Y(), 1, p[1].X(), p[1].Y(), 1, p[2].X(), p[2].Y());

    double delta = 0.5*deltaM.Determinant();

    std::vector<double> a, b, c;
    for(size_t i = 0; i < this->nodes.size(); ++i){
        size_t j = (i + 1) % 3;
        size_t k = (i + 2) % 3;

        a.push_back(p[j].X()*p[k].Y() - p[k].X()*p[j].Y());
        b.push_back(p[j].Y() - p[k].Y());
        c.push_back(p[k].X() - p[j].X());
    }

    double a0 = a[0];
    double a1 = a[1];
    double a2 = a[2];
    double b0 = b[0];
    double b1 = b[1];
    double b2 = b[2];
    double c0 = c[0];
    double c1 = c[1];
    double c2 = c[2];

    double x1 = points[0].X();
    double y1 = points[0].Y();
    double x2 = points[1].X();
    double y2 = points[1].Y();

    std::vector<double> Nf({
        t*(2*a0 + b0*x1 + b0*x2 + c0*y1 + c0*y2)*sqrt(x1*x1 - 2*x1*x2 + x2*x2 + y1*y1 - 2*y1*y2 + y2*y2)/(4*delta)
        ,
        0
        ,
        0
        ,
        t*(2*a0 + b0*x1 + b0*x2 + c0*y1 + c0*y2)*sqrt(x1*x1 - 2*x1*x2 + x2*x2 + y1*y1 - 2*y1*y2 + y2*y2)/(4*delta)
        ,
        t*(2*a1 + b1*x1 + b1*x2 + c1*y1 + c1*y2)*sqrt(x1*x1 - 2*x1*x2 + x2*x2 + y1*y1 - 2*y1*y2 + y2*y2)/(4*delta)
        ,
        0
        ,
        0
        ,
        t*(2*a1 + b1*x1 + b1*x2 + c1*y1 + c1*y2)*sqrt(x1*x1 - 2*x1*x2 + x2*x2 + y1*y1 - 2*y1*y2 + y2*y2)/(4*delta)
        ,
        t*(2*a2 + b2*x1 + b2*x2 + c2*y1 + c2*y2)*sqrt(x1*x1 - 2*x1*x2 + x2*x2 + y1*y1 - 2*y1*y2 + y2*y2)/(4*delta)
        ,
        0
        ,
        0
        ,
        t*(2*a2 + b2*x1 + b2*x2 + c2*y1 + c2*y2)*sqrt(x1*x1 - 2*x1*x2 + x2*x2 + y1*y1 - 2*y1*y2 + y2*y2)/(4*delta)
    });
    
    return Nf;
}

}
