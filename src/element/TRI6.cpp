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

#include "element/TRI6.hpp"
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

TRI6::TRI6(ElementShape s):
    MeshElementCommon2DTri<TRI6>(s.nodes){
    
    gp_Pnt p[3];
    for(size_t i = 0; i < 3; ++i){
        const auto& n = this->nodes[i];
        p[i] = n->point;
    }
    for(size_t i = 0; i < 3; ++i){
        size_t j = (i + 1) % 3;
        size_t k = (i + 2) % 3;

        a[i] = p[j].X()*p[k].Y() - p[k].X()*p[j].Y();
        b[i] = p[j].Y() - p[k].Y();
        c[i] = p[k].X() - p[j].X();
    }

    gp_Mat deltaM(1, p[0].X(), p[0].Y(), 1, p[1].X(), p[1].Y(), 1, p[2].X(), p[2].Y());

    this->delta = 0.5*std::abs(deltaM.Determinant());
}

std::vector<double> TRI6::get_k(const std::vector<double>& D, const double t) const{
    std::vector<double> k{
    t*(D[0]*b[0]*b[0] + D[2]*b[0]*c[0] + D[6]*b[0]*c[0] + D[8]*c[0]*c[0])/(4*delta)
    ,
    t*(D[1]*b[0]*c[0] + D[2]*b[0]*b[0] + D[7]*c[0]*c[0] + D[8]*b[0]*c[0])/(4*delta)
    ,
    t*(-D[0]*b[0]*b[1] - D[2]*b[0]*c[1] - D[6]*b[1]*c[0] - D[8]*c[0]*c[1])/(12*delta)
    ,
    t*(-D[1]*b[0]*c[1] - D[2]*b[0]*b[1] - D[7]*c[0]*c[1] - D[8]*b[1]*c[0])/(12*delta)
    ,
    t*(-D[0]*b[0]*b[2] - D[2]*b[0]*c[2] - D[6]*b[2]*c[0] - D[8]*c[0]*c[2])/(12*delta)
    ,
    t*(-D[1]*b[0]*c[2] - D[2]*b[0]*b[2] - D[7]*c[0]*c[2] - D[8]*b[2]*c[0])/(12*delta)
    ,
    t*(D[0]*b[0]*b[1] + D[2]*b[0]*c[1] + D[6]*b[1]*c[0] + D[8]*c[0]*c[1])/(3*delta)
    ,
    t*(D[1]*b[0]*c[1] + D[2]*b[0]*b[1] + D[7]*c[0]*c[1] + D[8]*b[1]*c[0])/(3*delta)
    ,
    0
    ,
    0
    ,
    t*(D[0]*b[0]*b[2] + D[2]*b[0]*c[2] + D[6]*b[2]*c[0] + D[8]*c[0]*c[2])/(3*delta)
    ,
    t*(D[1]*b[0]*c[2] + D[2]*b[0]*b[2] + D[7]*c[0]*c[2] + D[8]*b[2]*c[0])/(3*delta)
    ,
    t*(D[3]*b[0]*c[0] + D[5]*c[0]*c[0] + D[6]*b[0]*b[0] + D[8]*b[0]*c[0])/(4*delta)
    ,
    t*(D[4]*c[0]*c[0] + D[5]*b[0]*c[0] + D[7]*b[0]*c[0] + D[8]*b[0]*b[0])/(4*delta)
    ,
    t*(-D[3]*b[1]*c[0] - D[5]*c[0]*c[1] - D[6]*b[0]*b[1] - D[8]*b[0]*c[1])/(12*delta)
    ,
    t*(-D[4]*c[0]*c[1] - D[5]*b[1]*c[0] - D[7]*b[0]*c[1] - D[8]*b[0]*b[1])/(12*delta)
    ,
    t*(-D[3]*b[2]*c[0] - D[5]*c[0]*c[2] - D[6]*b[0]*b[2] - D[8]*b[0]*c[2])/(12*delta)
    ,
    t*(-D[4]*c[0]*c[2] - D[5]*b[2]*c[0] - D[7]*b[0]*c[2] - D[8]*b[0]*b[2])/(12*delta)
    ,
    t*(D[3]*b[1]*c[0] + D[5]*c[0]*c[1] + D[6]*b[0]*b[1] + D[8]*b[0]*c[1])/(3*delta)
    ,
    t*(D[4]*c[0]*c[1] + D[5]*b[1]*c[0] + D[7]*b[0]*c[1] + D[8]*b[0]*b[1])/(3*delta)
    ,
    0
    ,
    0
    ,
    t*(D[3]*b[2]*c[0] + D[5]*c[0]*c[2] + D[6]*b[0]*b[2] + D[8]*b[0]*c[2])/(3*delta)
    ,
    t*(D[4]*c[0]*c[2] + D[5]*b[2]*c[0] + D[7]*b[0]*c[2] + D[8]*b[0]*b[2])/(3*delta)
    ,
    t*(-D[0]*b[0]*b[1] - D[2]*b[1]*c[0] - D[6]*b[0]*c[1] - D[8]*c[0]*c[1])/(12*delta)
    ,
    t*(-D[1]*b[1]*c[0] - D[2]*b[0]*b[1] - D[7]*c[0]*c[1] - D[8]*b[0]*c[1])/(12*delta)
    ,
    t*(D[0]*b[1]*b[1] + D[2]*b[1]*c[1] + D[6]*b[1]*c[1] + D[8]*c[1]*c[1])/(4*delta)
    ,
    t*(D[1]*b[1]*c[1] + D[2]*b[1]*b[1] + D[7]*c[1]*c[1] + D[8]*b[1]*c[1])/(4*delta)
    ,
    t*(-D[0]*b[1]*b[2] - D[2]*b[1]*c[2] - D[6]*b[2]*c[1] - D[8]*c[1]*c[2])/(12*delta)
    ,
    t*(-D[1]*b[1]*c[2] - D[2]*b[1]*b[2] - D[7]*c[1]*c[2] - D[8]*b[2]*c[1])/(12*delta)
    ,
    t*(D[0]*b[0]*b[1] + D[2]*b[1]*c[0] + D[6]*b[0]*c[1] + D[8]*c[0]*c[1])/(3*delta)
    ,
    t*(D[1]*b[1]*c[0] + D[2]*b[0]*b[1] + D[7]*c[0]*c[1] + D[8]*b[0]*c[1])/(3*delta)
    ,
    t*(D[0]*b[1]*b[2] + D[2]*b[1]*c[2] + D[6]*b[2]*c[1] + D[8]*c[1]*c[2])/(3*delta)
    ,
    t*(D[1]*b[1]*c[2] + D[2]*b[1]*b[2] + D[7]*c[1]*c[2] + D[8]*b[2]*c[1])/(3*delta)
    ,
    0
    ,
    0
    ,
    t*(-D[3]*b[0]*c[1] - D[5]*c[0]*c[1] - D[6]*b[0]*b[1] - D[8]*b[1]*c[0])/(12*delta)
    ,
    t*(-D[4]*c[0]*c[1] - D[5]*b[0]*c[1] - D[7]*b[1]*c[0] - D[8]*b[0]*b[1])/(12*delta)
    ,
    t*(D[3]*b[1]*c[1] + D[5]*c[1]*c[1] + D[6]*b[1]*b[1] + D[8]*b[1]*c[1])/(4*delta)
    ,
    t*(D[4]*c[1]*c[1] + D[5]*b[1]*c[1] + D[7]*b[1]*c[1] + D[8]*b[1]*b[1])/(4*delta)
    ,
    t*(-D[3]*b[2]*c[1] - D[5]*c[1]*c[2] - D[6]*b[1]*b[2] - D[8]*b[1]*c[2])/(12*delta)
    ,
    t*(-D[4]*c[1]*c[2] - D[5]*b[2]*c[1] - D[7]*b[1]*c[2] - D[8]*b[1]*b[2])/(12*delta)
    ,
    t*(D[3]*b[0]*c[1] + D[5]*c[0]*c[1] + D[6]*b[0]*b[1] + D[8]*b[1]*c[0])/(3*delta)
    ,
    t*(D[4]*c[0]*c[1] + D[5]*b[0]*c[1] + D[7]*b[1]*c[0] + D[8]*b[0]*b[1])/(3*delta)
    ,
    t*(D[3]*b[2]*c[1] + D[5]*c[1]*c[2] + D[6]*b[1]*b[2] + D[8]*b[1]*c[2])/(3*delta)
    ,
    t*(D[4]*c[1]*c[2] + D[5]*b[2]*c[1] + D[7]*b[1]*c[2] + D[8]*b[1]*b[2])/(3*delta)
    ,
    0
    ,
    0
    ,
    t*(-D[0]*b[0]*b[2] - D[2]*b[2]*c[0] - D[6]*b[0]*c[2] - D[8]*c[0]*c[2])/(12*delta)
    ,
    t*(-D[1]*b[2]*c[0] - D[2]*b[0]*b[2] - D[7]*c[0]*c[2] - D[8]*b[0]*c[2])/(12*delta)
    ,
    t*(-D[0]*b[1]*b[2] - D[2]*b[2]*c[1] - D[6]*b[1]*c[2] - D[8]*c[1]*c[2])/(12*delta)
    ,
    t*(-D[1]*b[2]*c[1] - D[2]*b[1]*b[2] - D[7]*c[1]*c[2] - D[8]*b[1]*c[2])/(12*delta)
    ,
    t*(D[0]*b[2]*b[2] + D[2]*b[2]*c[2] + D[6]*b[2]*c[2] + D[8]*c[2]*c[2])/(4*delta)
    ,
    t*(D[1]*b[2]*c[2] + D[2]*b[2]*b[2] + D[7]*c[2]*c[2] + D[8]*b[2]*c[2])/(4*delta)
    ,
    0
    ,
    0
    ,
    t*(D[0]*b[1]*b[2] + D[2]*b[2]*c[1] + D[6]*b[1]*c[2] + D[8]*c[1]*c[2])/(3*delta)
    ,
    t*(D[1]*b[2]*c[1] + D[2]*b[1]*b[2] + D[7]*c[1]*c[2] + D[8]*b[1]*c[2])/(3*delta)
    ,
    t*(D[0]*b[0]*b[2] + D[2]*b[2]*c[0] + D[6]*b[0]*c[2] + D[8]*c[0]*c[2])/(3*delta)
    ,
    t*(D[1]*b[2]*c[0] + D[2]*b[0]*b[2] + D[7]*c[0]*c[2] + D[8]*b[0]*c[2])/(3*delta)
    ,
    t*(-D[3]*b[0]*c[2] - D[5]*c[0]*c[2] - D[6]*b[0]*b[2] - D[8]*b[2]*c[0])/(12*delta)
    ,
    t*(-D[4]*c[0]*c[2] - D[5]*b[0]*c[2] - D[7]*b[2]*c[0] - D[8]*b[0]*b[2])/(12*delta)
    ,
    t*(-D[3]*b[1]*c[2] - D[5]*c[1]*c[2] - D[6]*b[1]*b[2] - D[8]*b[2]*c[1])/(12*delta)
    ,
    t*(-D[4]*c[1]*c[2] - D[5]*b[1]*c[2] - D[7]*b[2]*c[1] - D[8]*b[1]*b[2])/(12*delta)
    ,
    t*(D[3]*b[2]*c[2] + D[5]*c[2]*c[2] + D[6]*b[2]*b[2] + D[8]*b[2]*c[2])/(4*delta)
    ,
    t*(D[4]*c[2]*c[2] + D[5]*b[2]*c[2] + D[7]*b[2]*c[2] + D[8]*b[2]*b[2])/(4*delta)
    ,
    0
    ,
    0
    ,
    t*(D[3]*b[1]*c[2] + D[5]*c[1]*c[2] + D[6]*b[1]*b[2] + D[8]*b[2]*c[1])/(3*delta)
    ,
    t*(D[4]*c[1]*c[2] + D[5]*b[1]*c[2] + D[7]*b[2]*c[1] + D[8]*b[1]*b[2])/(3*delta)
    ,
    t*(D[3]*b[0]*c[2] + D[5]*c[0]*c[2] + D[6]*b[0]*b[2] + D[8]*b[2]*c[0])/(3*delta)
    ,
    t*(D[4]*c[0]*c[2] + D[5]*b[0]*c[2] + D[7]*b[2]*c[0] + D[8]*b[0]*b[2])/(3*delta)
    ,
    t*(D[0]*b[0]*b[1] + D[2]*b[1]*c[0] + D[6]*b[0]*c[1] + D[8]*c[0]*c[1])/(3*delta)
    ,
    t*(D[1]*b[1]*c[0] + D[2]*b[0]*b[1] + D[7]*c[0]*c[1] + D[8]*b[0]*c[1])/(3*delta)
    ,
    t*(D[0]*b[0]*b[1] + D[2]*b[0]*c[1] + D[6]*b[1]*c[0] + D[8]*c[0]*c[1])/(3*delta)
    ,
    t*(D[1]*b[0]*c[1] + D[2]*b[0]*b[1] + D[7]*c[0]*c[1] + D[8]*b[1]*c[0])/(3*delta)
    ,
    0
    ,
    0
    ,
    t*(2*D[0]*b[0]*b[0] + 2*D[0]*b[0]*b[1] + 2*D[0]*b[1]*b[1] + 2*D[2]*b[0]*c[0] + D[2]*b[0]*c[1] + D[2]*b[1]*c[0] + 2*D[2]*b[1]*c[1] + 2*D[6]*b[0]*c[0] + D[6]*b[0]*c[1] + D[6]*b[1]*c[0] + 2*D[6]*b[1]*c[1] + 2*D[8]*c[0]*c[0] + 2*D[8]*c[0]*c[1] + 2*D[8]*c[1]*c[1])/(3*delta)
    ,
    t*(2*D[1]*b[0]*c[0] + D[1]*b[0]*c[1] + D[1]*b[1]*c[0] + 2*D[1]*b[1]*c[1] + 2*D[2]*b[0]*b[0] + 2*D[2]*b[0]*b[1] + 2*D[2]*b[1]*b[1] + 2*D[7]*c[0]*c[0] + 2*D[7]*c[0]*c[1] + 2*D[7]*c[1]*c[1] +2*D[8]*b[0]*c[0] + D[8]*b[0]*c[1] + D[8]*b[1]*c[0] + 2*D[8]*b[1]*c[1])/(3*delta)
    ,
    t*(D[0]*b[0]*b[1] + 2*D[0]*b[0]*b[2] + D[0]*b[1]*b[1] + D[0]*b[1]*b[2] + D[2]*b[0]*c[1] + 2*D[2]*b[0]*c[2] + D[2]*b[1]*c[1] + D[2]*b[1]*c[2] + D[6]*b[1]*c[0] + D[6]*b[1]*c[1] + 2*D[6]*b[2]*c[0] + D[6]*b[2]*c[1] + D[8]*c[0]*c[1] + 2*D[8]*c[0]*c[2] + D[8]*c[1]*c[1] + D[8]*c[1]*c[2])/(3*delta)
    ,
    t*(D[1]*b[0]*c[1] + 2*D[1]*b[0]*c[2] + D[1]*b[1]*c[1] + D[1]*b[1]*c[2] + D[2]*b[0]*b[1] + 2*D[2]*b[0]*b[2] + D[2]*b[1]*b[1] + D[2]*b[1]*b[2] + D[7]*c[0]*c[1] + 2*D[7]*c[0]*c[2] + D[7]*c[1]*c[1] + D[7]*c[1]*c[2] + D[8]*b[1]*c[0] + D[8]*b[1]*c[1] + 2*D[8]*b[2]*c[0] + D[8]*b[2]*c[1])/(3*delta)
    ,
    t*(D[0]*b[0]*b[0] + D[0]*b[0]*b[1] + D[0]*b[0]*b[2] + 2*D[0]*b[1]*b[2] + D[2]*b[0]*c[0] + D[2]*b[0]*c[2] + D[2]*b[1]*c[0] + 2*D[2]*b[1]*c[2] + D[6]*b[0]*c[0] + D[6]*b[0]*c[1] + D[6]*b[2]*c[0] + 2*D[6]*b[2]*c[1] + D[8]*c[0]*c[0] + D[8]*c[0]*c[1] + D[8]*c[0]*c[2] + 2*D[8]*c[1]*c[2])/(3*delta)
    ,
    t*(D[1]*b[0]*c[0] + D[1]*b[0]*c[2] + D[1]*b[1]*c[0] + 2*D[1]*b[1]*c[2] + D[2]*b[0]*b[0] + D[2]*b[0]*b[1] + D[2]*b[0]*b[2] + 2*D[2]*b[1]*b[2] + D[7]*c[0]*c[0] + D[7]*c[0]*c[1] + D[7]*c[0]*c[2] + 2*D[7]*c[1]*c[2] + D[8]*b[0]*c[0] + D[8]*b[0]*c[1] + D[8]*b[2]*c[0] + 2*D[8]*b[2]*c[1])/(3*delta)
    ,
    t*(D[3]*b[0]*c[1] + D[5]*c[0]*c[1] + D[6]*b[0]*b[1] + D[8]*b[1]*c[0])/(3*delta)
    ,
    t*(D[4]*c[0]*c[1] + D[5]*b[0]*c[1] + D[7]*b[1]*c[0] + D[8]*b[0]*b[1])/(3*delta)
    ,
    t*(D[3]*b[1]*c[0] + D[5]*c[0]*c[1] + D[6]*b[0]*b[1] + D[8]*b[0]*c[1])/(3*delta)
    ,
    t*(D[4]*c[0]*c[1] + D[5]*b[1]*c[0] + D[7]*b[0]*c[1] + D[8]*b[0]*b[1])/(3*delta)
    ,
    0
    ,
    0
    ,
    t*(2*D[3]*b[0]*c[0] + D[3]*b[0]*c[1] + D[3]*b[1]*c[0] + 2*D[3]*b[1]*c[1] + 2*D[5]*c[0]*c[0] + 2*D[5]*c[0]*c[1] + 2*D[5]*c[1]*c[1] + 2*D[6]*b[0]*b[0] + 2*D[6]*b[0]*b[1] + 2*D[6]*b[1]*b[1] +2*D[8]*b[0]*c[0] + D[8]*b[0]*c[1] + D[8]*b[1]*c[0] + 2*D[8]*b[1]*c[1])/(3*delta)
    ,
    t*(2*D[4]*c[0]*c[0] + 2*D[4]*c[0]*c[1] + 2*D[4]*c[1]*c[1] + 2*D[5]*b[0]*c[0] + D[5]*b[0]*c[1] + D[5]*b[1]*c[0] + 2*D[5]*b[1]*c[1] + 2*D[7]*b[0]*c[0] + D[7]*b[0]*c[1] + D[7]*b[1]*c[0] + 2*D[7]*b[1]*c[1] + 2*D[8]*b[0]*b[0] + 2*D[8]*b[0]*b[1] + 2*D[8]*b[1]*b[1])/(3*delta)
    ,
    t*(D[3]*b[1]*c[0] + D[3]*b[1]*c[1] + 2*D[3]*b[2]*c[0] + D[3]*b[2]*c[1] + D[5]*c[0]*c[1] + 2*D[5]*c[0]*c[2] + D[5]*c[1]*c[1] + D[5]*c[1]*c[2] + D[6]*b[0]*b[1] + 2*D[6]*b[0]*b[2] + D[6]*b[1]*b[1] + D[6]*b[1]*b[2] + D[8]*b[0]*c[1] + 2*D[8]*b[0]*c[2] + D[8]*b[1]*c[1] + D[8]*b[1]*c[2])/(3*delta)
    ,
    t*(D[4]*c[0]*c[1] + 2*D[4]*c[0]*c[2] + D[4]*c[1]*c[1] + D[4]*c[1]*c[2] + D[5]*b[1]*c[0] + D[5]*b[1]*c[1] + 2*D[5]*b[2]*c[0] + D[5]*b[2]*c[1] + D[7]*b[0]*c[1] + 2*D[7]*b[0]*c[2] + D[7]*b[1]*c[1] + D[7]*b[1]*c[2] + D[8]*b[0]*b[1] + 2*D[8]*b[0]*b[2] + D[8]*b[1]*b[1] + D[8]*b[1]*b[2])/(3*delta)
    ,
    t*(D[3]*b[0]*c[0] + D[3]*b[0]*c[1] + D[3]*b[2]*c[0] + 2*D[3]*b[2]*c[1] + D[5]*c[0]*c[0] + D[5]*c[0]*c[1] + D[5]*c[0]*c[2] + 2*D[5]*c[1]*c[2] + D[6]*b[0]*b[0] + D[6]*b[0]*b[1] + D[6]*b[0]*b[2] + 2*D[6]*b[1]*b[2] + D[8]*b[0]*c[0] + D[8]*b[0]*c[2] + D[8]*b[1]*c[0] + 2*D[8]*b[1]*c[2])/(3*delta)
    ,
    t*(D[4]*c[0]*c[0] + D[4]*c[0]*c[1] + D[4]*c[0]*c[2] + 2*D[4]*c[1]*c[2] + D[5]*b[0]*c[0] + D[5]*b[0]*c[1] + D[5]*b[2]*c[0] + 2*D[5]*b[2]*c[1] + D[7]*b[0]*c[0] + D[7]*b[0]*c[2] + D[7]*b[1]*c[0] + 2*D[7]*b[1]*c[2] + D[8]*b[0]*b[0] + D[8]*b[0]*b[1] + D[8]*b[0]*b[2] + 2*D[8]*b[1]*b[2])/(3*delta)
    ,
    0
    ,
    0
    ,
    t*(D[0]*b[1]*b[2] + D[2]*b[2]*c[1] + D[6]*b[1]*c[2] + D[8]*c[1]*c[2])/(3*delta)
    ,
    t*(D[1]*b[2]*c[1] + D[2]*b[1]*b[2] + D[7]*c[1]*c[2] + D[8]*b[1]*c[2])/(3*delta)
    ,
    t*(D[0]*b[1]*b[2] + D[2]*b[1]*c[2] + D[6]*b[2]*c[1] + D[8]*c[1]*c[2])/(3*delta)
    ,
    t*(D[1]*b[1]*c[2] + D[2]*b[1]*b[2] + D[7]*c[1]*c[2] + D[8]*b[2]*c[1])/(3*delta)
    ,
    t*(D[0]*b[0]*b[1] + 2*D[0]*b[0]*b[2] + D[0]*b[1]*b[1] + D[0]*b[1]*b[2] + D[2]*b[1]*c[0] + D[2]*b[1]*c[1] + 2*D[2]*b[2]*c[0] + D[2]*b[2]*c[1] + D[6]*b[0]*c[1] + 2*D[6]*b[0]*c[2] + D[6]*b[1]*c[1] + D[6]*b[1]*c[2] + D[8]*c[0]*c[1] + 2*D[8]*c[0]*c[2] + D[8]*c[1]*c[1] + D[8]*c[1]*c[2])/(3*delta)
    ,
    t*(D[1]*b[1]*c[0] + D[1]*b[1]*c[1] + 2*D[1]*b[2]*c[0] + D[1]*b[2]*c[1] + D[2]*b[0]*b[1] + 2*D[2]*b[0]*b[2] + D[2]*b[1]*b[1] + D[2]*b[1]*b[2] + D[7]*c[0]*c[1] + 2*D[7]*c[0]*c[2] + D[7]*c[1]*c[1] + D[7]*c[1]*c[2] + D[8]*b[0]*c[1] + 2*D[8]*b[0]*c[2] + D[8]*b[1]*c[1] + D[8]*b[1]*c[2])/(3*delta)
    ,
    t*(2*D[0]*b[1]*b[1] + 2*D[0]*b[1]*b[2] + 2*D[0]*b[2]*b[2] + 2*D[2]*b[1]*c[1] + D[2]*b[1]*c[2] + D[2]*b[2]*c[1] + 2*D[2]*b[2]*c[2] + 2*D[6]*b[1]*c[1] + D[6]*b[1]*c[2] + D[6]*b[2]*c[1] + 2*D[6]*b[2]*c[2] + 2*D[8]*c[1]*c[1] + 2*D[8]*c[1]*c[2] + 2*D[8]*c[2]*c[2])/(3*delta)
    ,
    t*(2*D[1]*b[1]*c[1] + D[1]*b[1]*c[2] + D[1]*b[2]*c[1] + 2*D[1]*b[2]*c[2] + 2*D[2]*b[1]*b[1] + 2*D[2]*b[1]*b[2] + 2*D[2]*b[2]*b[2] + 2*D[7]*c[1]*c[1] + 2*D[7]*c[1]*c[2] + 2*D[7]*c[2]*c[2] +2*D[8]*b[1]*c[1] + D[8]*b[1]*c[2] + D[8]*b[2]*c[1] + 2*D[8]*b[2]*c[2])/(3*delta)
    ,
    t*(2*D[0]*b[0]*b[1] + D[0]*b[0]*b[2] + D[0]*b[1]*b[2] + D[0]*b[2]*b[2] + 2*D[2]*b[1]*c[0] + D[2]*b[1]*c[2] + D[2]*b[2]*c[0] + D[2]*b[2]*c[2] + 2*D[6]*b[0]*c[1] + D[6]*b[0]*c[2] + D[6]*b[2]*c[1] + D[6]*b[2]*c[2] + 2*D[8]*c[0]*c[1] + D[8]*c[0]*c[2] + D[8]*c[1]*c[2] + D[8]*c[2]*c[2])/(3*delta)
    ,
    t*(2*D[1]*b[1]*c[0] + D[1]*b[1]*c[2] + D[1]*b[2]*c[0] + D[1]*b[2]*c[2] + 2*D[2]*b[0]*b[1] + D[2]*b[0]*b[2] + D[2]*b[1]*b[2] + D[2]*b[2]*b[2] + 2*D[7]*c[0]*c[1] + D[7]*c[0]*c[2] + D[7]*c[1]*c[2] + D[7]*c[2]*c[2] + 2*D[8]*b[0]*c[1] + D[8]*b[0]*c[2] + D[8]*b[2]*c[1] + D[8]*b[2]*c[2])/(3*delta)
    ,
    0
    ,
    0
    ,
    t*(D[3]*b[1]*c[2] + D[5]*c[1]*c[2] + D[6]*b[1]*b[2] + D[8]*b[2]*c[1])/(3*delta)
    ,
    t*(D[4]*c[1]*c[2] + D[5]*b[1]*c[2] + D[7]*b[2]*c[1] + D[8]*b[1]*b[2])/(3*delta)
    ,
    t*(D[3]*b[2]*c[1] + D[5]*c[1]*c[2] + D[6]*b[1]*b[2] + D[8]*b[1]*c[2])/(3*delta)
    ,
    t*(D[4]*c[1]*c[2] + D[5]*b[2]*c[1] + D[7]*b[1]*c[2] + D[8]*b[1]*b[2])/(3*delta)
    ,
    t*(D[3]*b[0]*c[1] + 2*D[3]*b[0]*c[2] + D[3]*b[1]*c[1] + D[3]*b[1]*c[2] + D[5]*c[0]*c[1] + 2*D[5]*c[0]*c[2] + D[5]*c[1]*c[1] + D[5]*c[1]*c[2] + D[6]*b[0]*b[1] + 2*D[6]*b[0]*b[2] + D[6]*b[1]*b[1] + D[6]*b[1]*b[2] + D[8]*b[1]*c[0] + D[8]*b[1]*c[1] + 2*D[8]*b[2]*c[0] + D[8]*b[2]*c[1])/(3*delta)
    ,
    t*(D[4]*c[0]*c[1] + 2*D[4]*c[0]*c[2] + D[4]*c[1]*c[1] + D[4]*c[1]*c[2] + D[5]*b[0]*c[1] + 2*D[5]*b[0]*c[2] + D[5]*b[1]*c[1] + D[5]*b[1]*c[2] + D[7]*b[1]*c[0] + D[7]*b[1]*c[1] + 2*D[7]*b[2]*c[0] + D[7]*b[2]*c[1] + D[8]*b[0]*b[1] + 2*D[8]*b[0]*b[2] + D[8]*b[1]*b[1] + D[8]*b[1]*b[2])/(3*delta)
    ,
    t*(2*D[3]*b[1]*c[1] + D[3]*b[1]*c[2] + D[3]*b[2]*c[1] + 2*D[3]*b[2]*c[2] + 2*D[5]*c[1]*c[1] + 2*D[5]*c[1]*c[2] + 2*D[5]*c[2]*c[2] + 2*D[6]*b[1]*b[1] + 2*D[6]*b[1]*b[2] + 2*D[6]*b[2]*b[2] +2*D[8]*b[1]*c[1] + D[8]*b[1]*c[2] + D[8]*b[2]*c[1] + 2*D[8]*b[2]*c[2])/(3*delta)
    ,
    t*(2*D[4]*c[1]*c[1] + 2*D[4]*c[1]*c[2] + 2*D[4]*c[2]*c[2] + 2*D[5]*b[1]*c[1] + D[5]*b[1]*c[2] + D[5]*b[2]*c[1] + 2*D[5]*b[2]*c[2] + 2*D[7]*b[1]*c[1] + D[7]*b[1]*c[2] + D[7]*b[2]*c[1] + 2*D[7]*b[2]*c[2] + 2*D[8]*b[1]*b[1] + 2*D[8]*b[1]*b[2] + 2*D[8]*b[2]*b[2])/(3*delta)
    ,
    t*(2*D[3]*b[0]*c[1] + D[3]*b[0]*c[2] + D[3]*b[2]*c[1] + D[3]*b[2]*c[2] + 2*D[5]*c[0]*c[1] + D[5]*c[0]*c[2] + D[5]*c[1]*c[2] + D[5]*c[2]*c[2] + 2*D[6]*b[0]*b[1] + D[6]*b[0]*b[2] + D[6]*b[1]*b[2] + D[6]*b[2]*b[2] + 2*D[8]*b[1]*c[0] + D[8]*b[1]*c[2] + D[8]*b[2]*c[0] + D[8]*b[2]*c[2])/(3*delta)
    ,
    t*(2*D[4]*c[0]*c[1] + D[4]*c[0]*c[2] + D[4]*c[1]*c[2] + D[4]*c[2]*c[2] + 2*D[5]*b[0]*c[1] + D[5]*b[0]*c[2] + D[5]*b[2]*c[1] + D[5]*b[2]*c[2] + 2*D[7]*b[1]*c[0] + D[7]*b[1]*c[2] + D[7]*b[2]*c[0] + D[7]*b[2]*c[2] + 2*D[8]*b[0]*b[1] + D[8]*b[0]*b[2] + D[8]*b[1]*b[2] + D[8]*b[2]*b[2])/(3*delta)
    ,
    t*(D[0]*b[0]*b[2] + D[2]*b[2]*c[0] + D[6]*b[0]*c[2] + D[8]*c[0]*c[2])/(3*delta)
    ,
    t*(D[1]*b[2]*c[0] + D[2]*b[0]*b[2] + D[7]*c[0]*c[2] + D[8]*b[0]*c[2])/(3*delta)
    ,
    0
    ,
    0
    ,
    t*(D[0]*b[0]*b[2] + D[2]*b[0]*c[2] + D[6]*b[2]*c[0] + D[8]*c[0]*c[2])/(3*delta)
    ,
    t*(D[1]*b[0]*c[2] + D[2]*b[0]*b[2] + D[7]*c[0]*c[2] + D[8]*b[2]*c[0])/(3*delta)
    ,
    t*(D[0]*b[0]*b[0] + D[0]*b[0]*b[1] + D[0]*b[0]*b[2] + 2*D[0]*b[1]*b[2] + D[2]*b[0]*c[0] + D[2]*b[0]*c[1] + D[2]*b[2]*c[0] + 2*D[2]*b[2]*c[1] + D[6]*b[0]*c[0] + D[6]*b[0]*c[2] + D[6]*b[1]*c[0] + 2*D[6]*b[1]*c[2] + D[8]*c[0]*c[0] + D[8]*c[0]*c[1] + D[8]*c[0]*c[2] + 2*D[8]*c[1]*c[2])/(3*delta)
    ,
    t*(D[1]*b[0]*c[0] + D[1]*b[0]*c[1] + D[1]*b[2]*c[0] + 2*D[1]*b[2]*c[1] + D[2]*b[0]*b[0] + D[2]*b[0]*b[1] + D[2]*b[0]*b[2] + 2*D[2]*b[1]*b[2] + D[7]*c[0]*c[0] + D[7]*c[0]*c[1] + D[7]*c[0]*c[2] + 2*D[7]*c[1]*c[2] + D[8]*b[0]*c[0] + D[8]*b[0]*c[2] + D[8]*b[1]*c[0] + 2*D[8]*b[1]*c[2])/(3*delta)
    ,
    t*(2*D[0]*b[0]*b[1] + D[0]*b[0]*b[2] + D[0]*b[1]*b[2] + D[0]*b[2]*b[2] + 2*D[2]*b[0]*c[1] + D[2]*b[0]*c[2] + D[2]*b[2]*c[1] + D[2]*b[2]*c[2] + 2*D[6]*b[1]*c[0] + D[6]*b[1]*c[2] + D[6]*b[2]*c[0] + D[6]*b[2]*c[2] + 2*D[8]*c[0]*c[1] + D[8]*c[0]*c[2] + D[8]*c[1]*c[2] + D[8]*c[2]*c[2])/(3*delta)
    ,
    t*(2*D[1]*b[0]*c[1] + D[1]*b[0]*c[2] + D[1]*b[2]*c[1] + D[1]*b[2]*c[2] + 2*D[2]*b[0]*b[1] + D[2]*b[0]*b[2] + D[2]*b[1]*b[2] + D[2]*b[2]*b[2] + 2*D[7]*c[0]*c[1] + D[7]*c[0]*c[2] + D[7]*c[1]*c[2] + D[7]*c[2]*c[2] + 2*D[8]*b[1]*c[0] + D[8]*b[1]*c[2] + D[8]*b[2]*c[0] + D[8]*b[2]*c[2])/(3*delta)
    ,
    t*(2*D[0]*b[0]*b[0] + 2*D[0]*b[0]*b[2] + 2*D[0]*b[2]*b[2] + 2*D[2]*b[0]*c[0] + D[2]*b[0]*c[2] + D[2]*b[2]*c[0] + 2*D[2]*b[2]*c[2] + 2*D[6]*b[0]*c[0] + D[6]*b[0]*c[2] + D[6]*b[2]*c[0] + 2*D[6]*b[2]*c[2] + 2*D[8]*c[0]*c[0] + 2*D[8]*c[0]*c[2] + 2*D[8]*c[2]*c[2])/(3*delta)
    ,
    t*(2*D[1]*b[0]*c[0] + D[1]*b[0]*c[2] + D[1]*b[2]*c[0] + 2*D[1]*b[2]*c[2] + 2*D[2]*b[0]*b[0] + 2*D[2]*b[0]*b[2] + 2*D[2]*b[2]*b[2] + 2*D[7]*c[0]*c[0] + 2*D[7]*c[0]*c[2] + 2*D[7]*c[2]*c[2] +2*D[8]*b[0]*c[0] + D[8]*b[0]*c[2] + D[8]*b[2]*c[0] + 2*D[8]*b[2]*c[2])/(3*delta)
    ,
    t*(D[3]*b[0]*c[2] + D[5]*c[0]*c[2] + D[6]*b[0]*b[2] + D[8]*b[2]*c[0])/(3*delta)
    ,
    t*(D[4]*c[0]*c[2] + D[5]*b[0]*c[2] + D[7]*b[2]*c[0] + D[8]*b[0]*b[2])/(3*delta)
    ,
    0
    ,
    0
    ,
    t*(D[3]*b[2]*c[0] + D[5]*c[0]*c[2] + D[6]*b[0]*b[2] + D[8]*b[0]*c[2])/(3*delta)
    ,
    t*(D[4]*c[0]*c[2] + D[5]*b[2]*c[0] + D[7]*b[0]*c[2] + D[8]*b[0]*b[2])/(3*delta)
    ,
    t*(D[3]*b[0]*c[0] + D[3]*b[0]*c[2] + D[3]*b[1]*c[0] + 2*D[3]*b[1]*c[2] + D[5]*c[0]*c[0] + D[5]*c[0]*c[1] + D[5]*c[0]*c[2] + 2*D[5]*c[1]*c[2] + D[6]*b[0]*b[0] + D[6]*b[0]*b[1] + D[6]*b[0]*b[2] + 2*D[6]*b[1]*b[2] + D[8]*b[0]*c[0] + D[8]*b[0]*c[1] + D[8]*b[2]*c[0] + 2*D[8]*b[2]*c[1])/(3*delta)
    ,
    t*(D[4]*c[0]*c[0] + D[4]*c[0]*c[1] + D[4]*c[0]*c[2] + 2*D[4]*c[1]*c[2] + D[5]*b[0]*c[0] + D[5]*b[0]*c[2] + D[5]*b[1]*c[0] + 2*D[5]*b[1]*c[2] + D[7]*b[0]*c[0] + D[7]*b[0]*c[1] + D[7]*b[2]*c[0] + 2*D[7]*b[2]*c[1] + D[8]*b[0]*b[0] + D[8]*b[0]*b[1] + D[8]*b[0]*b[2] + 2*D[8]*b[1]*b[2])/(3*delta)
    ,
    t*(2*D[3]*b[1]*c[0] + D[3]*b[1]*c[2] + D[3]*b[2]*c[0] + D[3]*b[2]*c[2] + 2*D[5]*c[0]*c[1] + D[5]*c[0]*c[2] + D[5]*c[1]*c[2] + D[5]*c[2]*c[2] + 2*D[6]*b[0]*b[1] + D[6]*b[0]*b[2] + D[6]*b[1]*b[2] + D[6]*b[2]*b[2] + 2*D[8]*b[0]*c[1] + D[8]*b[0]*c[2] + D[8]*b[2]*c[1] + D[8]*b[2]*c[2])/(3*delta)
    ,
    t*(2*D[4]*c[0]*c[1] + D[4]*c[0]*c[2] + D[4]*c[1]*c[2] + D[4]*c[2]*c[2] + 2*D[5]*b[1]*c[0] + D[5]*b[1]*c[2] + D[5]*b[2]*c[0] + D[5]*b[2]*c[2] + 2*D[7]*b[0]*c[1] + D[7]*b[0]*c[2] + D[7]*b[2]*c[1] + D[7]*b[2]*c[2] + 2*D[8]*b[0]*b[1] + D[8]*b[0]*b[2] + D[8]*b[1]*b[2] + D[8]*b[2]*b[2])/(3*delta)
    ,
    t*(2*D[3]*b[0]*c[0] + D[3]*b[0]*c[2] + D[3]*b[2]*c[0] + 2*D[3]*b[2]*c[2] + 2*D[5]*c[0]*c[0] + 2*D[5]*c[0]*c[2] + 2*D[5]*c[2]*c[2] + 2*D[6]*b[0]*b[0] + 2*D[6]*b[0]*b[2] + 2*D[6]*b[2]*b[2] +2*D[8]*b[0]*c[0] + D[8]*b[0]*c[2] + D[8]*b[2]*c[0] + 2*D[8]*b[2]*c[2])/(3*delta)
    ,
    t*(2*D[4]*c[0]*c[0] + 2*D[4]*c[0]*c[2] + 2*D[4]*c[2]*c[2] + 2*D[5]*b[0]*c[0] + D[5]*b[0]*c[2] + D[5]*b[2]*c[0] + 2*D[5]*b[2]*c[2] + 2*D[7]*b[0]*c[0] + D[7]*b[0]*c[2] + D[7]*b[2]*c[0] + 2*D[7]*b[2]*c[2] + 2*D[8]*b[0]*b[0] + 2*D[8]*b[0]*b[2] + 2*D[8]*b[2]*b[2])/(3*delta)
    };
    return k;
}

std::vector<double> TRI6::get_DB(const std::vector<double>& D, const gp_Pnt& point) const{
    const double x = point.X();
    const double y = point.Y();

    const double L0 = this->L(x, y, 0);
    const double L1 = this->L(x, y, 1);
    const double L2 = this->L(x, y, 2);

    std::vector<double> DB{
    (-D[0]*b[0] - D[2]*c[0] + 4*(D[0]*b[0] + D[2]*c[0])*L0)/(2*delta)
    ,
    (-D[1]*c[0] - D[2]*b[0] + 4*(D[1]*c[0] + D[2]*b[0])*L0)/(2*delta)
    ,
    (-D[0]*b[1] - D[2]*c[1] + 4*(D[0]*b[1] + D[2]*c[1])*L1)/(2*delta)
    ,
    (-D[1]*c[1] - D[2]*b[1] + 4*(D[1]*c[1] + D[2]*b[1])*L1)/(2*delta)
    ,
    (-D[0]*b[2] - D[2]*c[2] + 4*(D[0]*b[2] + D[2]*c[2])*L2)/(2*delta)
    ,
    (-D[1]*c[2] - D[2]*b[2] + 4*(D[1]*c[2] + D[2]*b[2])*L2)/(2*delta)
    ,
    2*((D[0]*b[0] + D[2]*c[0])*L1 + (D[0]*b[1] + D[2]*c[1])*L0)/delta
    ,
    2*((D[1]*c[0] + D[2]*b[0])*L1 + (D[1]*c[1] + D[2]*b[1])*L0)/delta
    ,
    2*((D[0]*b[1] + D[2]*c[1])*L2 + (D[0]*b[2] + D[2]*c[2])*L1)/delta
    ,
    2*((D[1]*c[1] + D[2]*b[1])*L2 + (D[1]*c[2] + D[2]*b[2])*L1)/delta
    ,
    2*((D[0]*b[0] + D[2]*c[0])*L2 + (D[0]*b[2] + D[2]*c[2])*L0)/delta
    ,
    2*((D[1]*c[0] + D[2]*b[0])*L2 + (D[1]*c[2] + D[2]*b[2])*L0)/delta
    ,
    (-D[3]*b[0] - D[5]*c[0] + 4*(D[3]*b[0] + D[5]*c[0])*L0)/(2*delta)
    ,
    (-D[4]*c[0] - D[5]*b[0] + 4*(D[4]*c[0] + D[5]*b[0])*L0)/(2*delta)
    ,
    (-D[3]*b[1] - D[5]*c[1] + 4*(D[3]*b[1] + D[5]*c[1])*L1)/(2*delta)
    ,
    (-D[4]*c[1] - D[5]*b[1] + 4*(D[4]*c[1] + D[5]*b[1])*L1)/(2*delta)
    ,
    (-D[3]*b[2] - D[5]*c[2] + 4*(D[3]*b[2] + D[5]*c[2])*L2)/(2*delta)
    ,
    (-D[4]*c[2] - D[5]*b[2] + 4*(D[4]*c[2] + D[5]*b[2])*L2)/(2*delta)
    ,
    2*((D[3]*b[0] + D[5]*c[0])*L1 + (D[3]*b[1] + D[5]*c[1])*L0)/delta
    ,
    2*((D[4]*c[0] + D[5]*b[0])*L1 + (D[4]*c[1] + D[5]*b[1])*L0)/delta
    ,
    2*((D[3]*b[1] + D[5]*c[1])*L2 + (D[3]*b[2] + D[5]*c[2])*L1)/delta
    ,
    2*((D[4]*c[1] + D[5]*b[1])*L2 + (D[4]*c[2] + D[5]*b[2])*L1)/delta
    ,
    2*((D[3]*b[0] + D[5]*c[0])*L2 + (D[3]*b[2] + D[5]*c[2])*L0)/delta
    ,
    2*((D[4]*c[0] + D[5]*b[0])*L2 + (D[4]*c[2] + D[5]*b[2])*L0)/delta
    ,
    (-D[6]*b[0] - D[8]*c[0] + 4*(D[6]*b[0] + D[8]*c[0])*L0)/(2*delta)
    ,
    (-D[7]*c[0] - D[8]*b[0] + 4*(D[7]*c[0] + D[8]*b[0])*L0)/(2*delta)
    ,
    (-D[6]*b[1] - D[8]*c[1] + 4*(D[6]*b[1] + D[8]*c[1])*L1)/(2*delta)
    ,
    (-D[7]*c[1] - D[8]*b[1] + 4*(D[7]*c[1] + D[8]*b[1])*L1)/(2*delta)
    ,
    (-D[6]*b[2] - D[8]*c[2] + 4*(D[6]*b[2] + D[8]*c[2])*L2)/(2*delta)
    ,
    (-D[7]*c[2] - D[8]*b[2] + 4*(D[7]*c[2] + D[8]*b[2])*L2)/(2*delta)
    ,
    2*((D[6]*b[0] + D[8]*c[0])*L1 + (D[6]*b[1] + D[8]*c[1])*L0)/delta
    ,
    2*((D[7]*c[0] + D[8]*b[0])*L1 + (D[7]*c[1] + D[8]*b[1])*L0)/delta
    ,
    2*((D[6]*b[1] + D[8]*c[1])*L2 + (D[6]*b[2] + D[8]*c[2])*L1)/delta
    ,
    2*((D[7]*c[1] + D[8]*b[1])*L2 + (D[7]*c[2] + D[8]*b[2])*L1)/delta
    ,
    2*((D[6]*b[0] + D[8]*c[0])*L2 + (D[6]*b[2] + D[8]*c[2])*L0)/delta
    ,
    2*((D[7]*c[0] + D[8]*b[0])*L2 + (D[7]*c[2] + D[8]*b[2])*L0)/delta
    };

    return DB;
}
std::vector<double> TRI6::get_Nf(const double t, const std::vector<gp_Pnt>& points) const{
    const double x[]{points[0].X(), points[1].X()};
    const double y[]{points[0].Y(), points[1].Y()};

    std::vector<double> Nf{
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[0]*a[0] + 6*a[0]*b[0]*x[0] + 6*a[0]*b[0]*x[1] + 6*a[0]*c[0]*y[0] + 6*a[0]*c[0]*y[1] + 2*b[0]*b[0]*x[0]*x[0] + 2*b[0]*b[0]*x[0]*x[1] + 2*b[0]*b[0]*x[1]*x[1] + 4*b[0]*c[0]*x[0]*y[0] + 2*b[0]*c[0]*x[0]*y[1] + 2*b[0]*c[0]*x[1]*y[0] + 4*b[0]*c[0]*x[1]*y[1] + 2*c[0]*c[0]*y[0]*y[0] + 2*c[0]*c[0]*y[0]*y[1] + 2*c[0]*c[0]*y[1]*y[1] - 3*delta*(2*a[0] + b[0]*x[0] + b[0]*x[1] + c[0]*y[0] + c[0]*y[1]))/(12*delta*delta)
    ,
    0
    ,
    0
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[0]*a[0] + 6*a[0]*b[0]*x[0] + 6*a[0]*b[0]*x[1] + 6*a[0]*c[0]*y[0] + 6*a[0]*c[0]*y[1] + 2*b[0]*b[0]*x[0]*x[0] + 2*b[0]*b[0]*x[0]*x[1] + 2*b[0]*b[0]*x[1]*x[1] + 4*b[0]*c[0]*x[0]*y[0] + 2*b[0]*c[0]*x[0]*y[1] + 2*b[0]*c[0]*x[1]*y[0] + 4*b[0]*c[0]*x[1]*y[1] + 2*c[0]*c[0]*y[0]*y[0] + 2*c[0]*c[0]*y[0]*y[1] + 2*c[0]*c[0]*y[1]*y[1] - 3*delta*(2*a[0] + b[0]*x[0] + b[0]*x[1] + c[0]*y[0] + c[0]*y[1]))/(12*delta*delta)
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[1]*a[1] + 6*a[1]*b[1]*x[0] + 6*a[1]*b[1]*x[1] + 6*a[1]*c[1]*y[0] + 6*a[1]*c[1]*y[1] + 2*b[1]*b[1]*x[0]*x[0] + 2*b[1]*b[1]*x[0]*x[1] + 2*b[1]*b[1]*x[1]*x[1] + 4*b[1]*c[1]*x[0]*y[0] + 2*b[1]*c[1]*x[0]*y[1] + 2*b[1]*c[1]*x[1]*y[0] + 4*b[1]*c[1]*x[1]*y[1] + 2*c[1]*c[1]*y[0]*y[0] + 2*c[1]*c[1]*y[0]*y[1] + 2*c[1]*c[1]*y[1]*y[1] - 3*delta*(2*a[1] + b[1]*x[0] + b[1]*x[1] + c[1]*y[0] + c[1]*y[1]))/(12*delta*delta)
    ,
    0
    ,
    0
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[1]*a[1] + 6*a[1]*b[1]*x[0] + 6*a[1]*b[1]*x[1] + 6*a[1]*c[1]*y[0] + 6*a[1]*c[1]*y[1] + 2*b[1]*b[1]*x[0]*x[0] + 2*b[1]*b[1]*x[0]*x[1] + 2*b[1]*b[1]*x[1]*x[1] + 4*b[1]*c[1]*x[0]*y[0] + 2*b[1]*c[1]*x[0]*y[1] + 2*b[1]*c[1]*x[1]*y[0] + 4*b[1]*c[1]*x[1]*y[1] + 2*c[1]*c[1]*y[0]*y[0] + 2*c[1]*c[1]*y[0]*y[1] + 2*c[1]*c[1]*y[1]*y[1] - 3*delta*(2*a[1] + b[1]*x[0] + b[1]*x[1] + c[1]*y[0] + c[1]*y[1]))/(12*delta*delta)
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[2]*a[2] + 6*a[2]*b[2]*x[0] + 6*a[2]*b[2]*x[1] + 6*a[2]*c[2]*y[0] + 6*a[2]*c[2]*y[1] + 2*b[2]*b[2]*x[0]*x[0] + 2*b[2]*b[2]*x[0]*x[1] + 2*b[2]*b[2]*x[1]*x[1] + 4*b[2]*c[2]*x[0]*y[0] + 2*b[2]*c[2]*x[0]*y[1] + 2*b[2]*c[2]*x[1]*y[0] + 4*b[2]*c[2]*x[1]*y[1] + 2*c[2]*c[2]*y[0]*y[0] + 2*c[2]*c[2]*y[0]*y[1] + 2*c[2]*c[2]*y[1]*y[1] - 3*delta*(2*a[2] + b[2]*x[0] + b[2]*x[1] + c[2]*y[0] + c[2]*y[1]))/(12*delta*delta)
    ,
    0
    ,
    0
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[2]*a[2] + 6*a[2]*b[2]*x[0] + 6*a[2]*b[2]*x[1] + 6*a[2]*c[2]*y[0] + 6*a[2]*c[2]*y[1] + 2*b[2]*b[2]*x[0]*x[0] + 2*b[2]*b[2]*x[0]*x[1] + 2*b[2]*b[2]*x[1]*x[1] + 4*b[2]*c[2]*x[0]*y[0] + 2*b[2]*c[2]*x[0]*y[1] + 2*b[2]*c[2]*x[1]*y[0] + 4*b[2]*c[2]*x[1]*y[1] + 2*c[2]*c[2]*y[0]*y[0] + 2*c[2]*c[2]*y[0]*y[1] + 2*c[2]*c[2]*y[1]*y[1] - 3*delta*(2*a[2] + b[2]*x[0] + b[2]*x[1] + c[2]*y[0] + c[2]*y[1]))/(12*delta*delta)
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[0]*a[1] + 3*a[0]*b[1]*x[0] + 3*a[0]*b[1]*x[1] + 3*a[0]*c[1]*y[0] + 3*a[0]*c[1]*y[1] + 3*a[1]*b[0]*x[0] + 3*a[1]*b[0]*x[1] + 3*a[1]*c[0]*y[0] + 3*a[1]*c[0]*y[1] + 2*b[0]*b[1]*x[0]*x[0] + 2*b[0]*b[1]*x[0]*x[1] + 2*b[0]*b[1]*x[1]*x[1] + 2*b[0]*c[1]*x[0]*y[0] + b[0]*c[1]*x[0]*y[1] + b[0]*c[1]*x[1]*y[0] + 2*b[0]*c[1]*x[1]*y[1] + 2*b[1]*c[0]*x[0]*y[0] + b[1]*c[0]*x[0]*y[1] + b[1]*c[0]*x[1]*y[0] + 2*b[1]*c[0]*x[1]*y[1] + 2*c[0]*c[1]*y[0]*y[0] + 2*c[0]*c[1]*y[0]*y[1] + 2*c[0]*c[1]*y[1]*y[1])/(6*delta*delta)
    ,
    0
    ,
    0
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[0]*a[1] + 3*a[0]*b[1]*x[0] + 3*a[0]*b[1]*x[1] + 3*a[0]*c[1]*y[0] + 3*a[0]*c[1]*y[1] + 3*a[1]*b[0]*x[0] + 3*a[1]*b[0]*x[1] + 3*a[1]*c[0]*y[0] + 3*a[1]*c[0]*y[1] + 2*b[0]*b[1]*x[0]*x[0] + 2*b[0]*b[1]*x[0]*x[1] + 2*b[0]*b[1]*x[1]*x[1] + 2*b[0]*c[1]*x[0]*y[0] + b[0]*c[1]*x[0]*y[1] + b[0]*c[1]*x[1]*y[0] + 2*b[0]*c[1]*x[1]*y[1] + 2*b[1]*c[0]*x[0]*y[0] + b[1]*c[0]*x[0]*y[1] + b[1]*c[0]*x[1]*y[0] + 2*b[1]*c[0]*x[1]*y[1] + 2*c[0]*c[1]*y[0]*y[0] + 2*c[0]*c[1]*y[0]*y[1] + 2*c[0]*c[1]*y[1]*y[1])/(6*delta*delta)
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[1]*a[2] + 3*a[1]*b[2]*x[0] + 3*a[1]*b[2]*x[1] + 3*a[1]*c[2]*y[0] + 3*a[1]*c[2]*y[1] + 3*a[2]*b[1]*x[0] + 3*a[2]*b[1]*x[1] + 3*a[2]*c[1]*y[0] + 3*a[2]*c[1]*y[1] + 2*b[1]*b[2]*x[0]*x[0] + 2*b[1]*b[2]*x[0]*x[1] + 2*b[1]*b[2]*x[1]*x[1] + 2*b[1]*c[2]*x[0]*y[0] + b[1]*c[2]*x[0]*y[1] + b[1]*c[2]*x[1]*y[0] + 2*b[1]*c[2]*x[1]*y[1] + 2*b[2]*c[1]*x[0]*y[0] + b[2]*c[1]*x[0]*y[1] + b[2]*c[1]*x[1]*y[0] + 2*b[2]*c[1]*x[1]*y[1] + 2*c[1]*c[2]*y[0]*y[0] + 2*c[1]*c[2]*y[0]*y[1] + 2*c[1]*c[2]*y[1]*y[1])/(6*delta*delta)
    ,
    0
    ,
    0
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[1]*a[2] + 3*a[1]*b[2]*x[0] + 3*a[1]*b[2]*x[1] + 3*a[1]*c[2]*y[0] + 3*a[1]*c[2]*y[1] + 3*a[2]*b[1]*x[0] + 3*a[2]*b[1]*x[1] + 3*a[2]*c[1]*y[0] + 3*a[2]*c[1]*y[1] + 2*b[1]*b[2]*x[0]*x[0] + 2*b[1]*b[2]*x[0]*x[1] + 2*b[1]*b[2]*x[1]*x[1] + 2*b[1]*c[2]*x[0]*y[0] + b[1]*c[2]*x[0]*y[1] + b[1]*c[2]*x[1]*y[0] + 2*b[1]*c[2]*x[1]*y[1] + 2*b[2]*c[1]*x[0]*y[0] + b[2]*c[1]*x[0]*y[1] + b[2]*c[1]*x[1]*y[0] + 2*b[2]*c[1]*x[1]*y[1] + 2*c[1]*c[2]*y[0]*y[0] + 2*c[1]*c[2]*y[0]*y[1] + 2*c[1]*c[2]*y[1]*y[1])/(6*delta*delta)
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[0]*a[2] + 3*a[0]*b[2]*x[0] + 3*a[0]*b[2]*x[1] + 3*a[0]*c[2]*y[0] + 3*a[0]*c[2]*y[1] + 3*a[2]*b[0]*x[0] + 3*a[2]*b[0]*x[1] + 3*a[2]*c[0]*y[0] + 3*a[2]*c[0]*y[1] + 2*b[0]*b[2]*x[0]*x[0] + 2*b[0]*b[2]*x[0]*x[1] + 2*b[0]*b[2]*x[1]*x[1] + 2*b[0]*c[2]*x[0]*y[0] + b[0]*c[2]*x[0]*y[1] + b[0]*c[2]*x[1]*y[0] + 2*b[0]*c[2]*x[1]*y[1] + 2*b[2]*c[0]*x[0]*y[0] + b[2]*c[0]*x[0]*y[1] + b[2]*c[0]*x[1]*y[0] + 2*b[2]*c[0]*x[1]*y[1] + 2*c[0]*c[2]*y[0]*y[0] + 2*c[0]*c[2]*y[0]*y[1] + 2*c[0]*c[2]*y[1]*y[1])/(6*delta*delta)
    ,
    0
    ,
    0
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[0]*a[2] + 3*a[0]*b[2]*x[0] + 3*a[0]*b[2]*x[1] + 3*a[0]*c[2]*y[0] + 3*a[0]*c[2]*y[1] + 3*a[2]*b[0]*x[0] + 3*a[2]*b[0]*x[1] + 3*a[2]*c[0]*y[0] + 3*a[2]*c[0]*y[1] + 2*b[0]*b[2]*x[0]*x[0] + 2*b[0]*b[2]*x[0]*x[1] + 2*b[0]*b[2]*x[1]*x[1] + 2*b[0]*c[2]*x[0]*y[0] + b[0]*c[2]*x[0]*y[1] + b[0]*c[2]*x[1]*y[0] + 2*b[0]*c[2]*x[1]*y[1] + 2*b[2]*c[0]*x[0]*y[0] + b[2]*c[0]*x[0]*y[1] + b[2]*c[0]*x[1]*y[0] + 2*b[2]*c[0]*x[1]*y[1] + 2*c[0]*c[2]*y[0]*y[0] + 2*c[0]*c[2]*y[0]*y[1] + 2*c[0]*c[2]*y[1]*y[1])/(6*delta*delta)
    };

    return Nf;
}

}
