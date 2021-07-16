/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   Foobar is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Foobar is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with SolidPrep.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "element/GT9.hpp"
#include "cblas.h"
#include "logger.hpp"
#include "project_data.hpp"
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>

namespace element{

GT9::GT9(ElementShape s, ProjectData* data):
    MeshElement(s.nodes), mat(data->material.get()), t(data->thickness){}

std::vector<double> GT9::get_k() const{

    size_t N = this->nodes.size();

    std::vector<gp_Pnt> p;
    for(auto n:this->nodes){
        p.push_back(n->point);
    }

    gp_Mat deltaM(1, p[0].X(), p[0].Y(), 1, p[1].X(), p[1].Y(), 1, p[2].X(), p[2].Y());

    double delta = 0.5*std::abs(deltaM.Determinant());

    std::vector<double> a, b, c;
    for(size_t i = 0; i < N; ++i){
        size_t j = (i + 1) % 3;
        size_t k = (i + 2) % 3;

        a.push_back(p[j].X()*p[k].Y() - p[k].X()*p[j].Y());
        b.push_back(p[j].Y() - p[k].Y());
        c.push_back(p[k].X() - p[j].X());
    }

    double b0 = b[0];
    double b1 = b[1];
    double b2 = b[2];
    double c0 = c[0];
    double c1 = c[1];
    double c2 = c[2];
    
    auto D = this->mat->stiffness_2D();

    double d0 = D[0];
    double d1 = D[1];
    double d2 = D[2];
    double d3 = D[3];
    double d4 = D[4];
    double d5 = D[5];
    double d6 = D[6];
    double d7 = D[7];
    double d8 = D[8];

    std::vector<double> K{
this->t*(b0*b0*d0 + b0*c0*d2 + b0*c0*d6 + c0*c0*d8)/(4*delta)
,
this->t*(b0*b0*d2 + b0*c0*d1 + b0*c0*d8 + c0*c0*d7)/(4*delta)
,
this->t*0.125f*(-(b0*b0*b1*d0 + b0*b0*c1*d2 + b0*b1*c0*d2 + b0*b1*c0*d6 + b0*c0*c1*d1 + b0*c0*c1*d8 + b1*c0*c0*d8 + c0*c0*c1*d7)/3 + (b0*b0*b2*d0 + b0*b0*c2*d2 + b0*b2*c0*d2 + b0*b2*c0*d6 + b0*c0*c2*d1 + b0*c0*c2*d8 + b2*c0*c0*d8 + c0*c0*c2*d7)/3)/delta
,
this->t*(b0*b1*d0 + b0*c1*d2 + b1*c0*d6 + c0*c1*d8)/(4*delta)
,
this->t*(b0*b1*d2 + b0*c1*d1 + b1*c0*d8 + c0*c1*d7)/(4*delta)
,
this->t*0.125f*((b0*b0*b1*d0 + b0*b0*c1*d2 + b0*b1*c0*d2 + b0*b1*c0*d6 + b0*c0*c1*d1 + b0*c0*c1*d8 + b1*c0*c0*d8 + c0*c0*c1*d7)/3 - (b0*b1*b2*d0 + b0*b1*c2*d2 + b0*b2*c1*d2 + b0*c1*c2*d1 + b1*b2*c0*d6 + b1*c0*c2*d8 + b2*c0*c1*d8 + c0*c1*c2*d7)/3)/delta
,
this->t*(b0*b2*d0 + b0*c2*d2 + b2*c0*d6 + c0*c2*d8)/(4*delta)
,
this->t*(b0*b2*d2 + b0*c2*d1 + b2*c0*d8 + c0*c2*d7)/(4*delta)
,
this->t*0.125f*(-(b0*b0*b2*d0 + b0*b0*c2*d2 + b0*b2*c0*d2 + b0*b2*c0*d6 + b0*c0*c2*d1 + b0*c0*c2*d8 + b2*c0*c0*d8 + c0*c0*c2*d7)/3 + (b0*b1*b2*d0 + b0*b1*c2*d2 + b0*b2*c1*d2 + b0*c1*c2*d1 + b1*b2*c0*d6 + b1*c0*c2*d8 + b2*c0*c1*d8 + c0*c1*c2*d7)/3)/delta
,
this->t*(b0*b0*d6 + b0*c0*d3 + b0*c0*d8 + c0*c0*d5)/(4*delta)
,
this->t*(b0*b0*d8 + b0*c0*d5 + b0*c0*d7 + c0*c0*d4)/(4*delta)
,
this->t*0.125f*(-(b0*b0*b1*d6 + b0*b0*c1*d8 + b0*b1*c0*d3 + b0*b1*c0*d8 + b0*c0*c1*d5 + b0*c0*c1*d7 + b1*c0*c0*d5 + c0*c0*c1*d4)/3 + (b0*b0*b2*d6 + b0*b0*c2*d8 + b0*b2*c0*d3 + b0*b2*c0*d8 + b0*c0*c2*d5 + b0*c0*c2*d7 + b2*c0*c0*d5 + c0*c0*c2*d4)/3)/delta
,
this->t*(b0*b1*d6 + b0*c1*d8 + b1*c0*d3 + c0*c1*d5)/(4*delta)
,
this->t*(b0*b1*d8 + b0*c1*d7 + b1*c0*d5 + c0*c1*d4)/(4*delta)
,
this->t*0.125f*((b0*b0*b1*d6 + b0*b0*c1*d8 + b0*b1*c0*d3 + b0*b1*c0*d8 + b0*c0*c1*d5 + b0*c0*c1*d7 + b1*c0*c0*d5 + c0*c0*c1*d4)/3 - (b0*b1*b2*d6 + b0*b1*c2*d8 + b0*b2*c1*d8 + b0*c1*c2*d7 + b1*b2*c0*d3 + b1*c0*c2*d5 + b2*c0*c1*d5 + c0*c1*c2*d4)/3)/delta
,
this->t*(b0*b2*d6 + b0*c2*d8 + b2*c0*d3 + c0*c2*d5)/(4*delta)
,
this->t*(b0*b2*d8 + b0*c2*d7 + b2*c0*d5 + c0*c2*d4)/(4*delta)
,
this->t*0.125f*(-(b0*b0*b2*d6 + b0*b0*c2*d8 + b0*b2*c0*d3 + b0*b2*c0*d8 + b0*c0*c2*d5 + b0*c0*c2*d7 + b2*c0*c0*d5 + c0*c0*c2*d4)/3 + (b0*b1*b2*d6 + b0*b1*c2*d8 + b0*b2*c1*d8 + b0*c1*c2*d7 + b1*b2*c0*d3 + b1*c0*c2*d5 + b2*c0*c1*d5 + c0*c1*c2*d4)/3)/delta
,
this->t*0.125f*(-(b0*b0*b1*d0 + b0*b0*c1*d6 + b0*b1*c0*d2 + b0*b1*c0*d6 + b0*c0*c1*d3 + b0*c0*c1*d8 + b1*c0*c0*d8 + c0*c0*c1*d5)/3 + (b0*b0*b2*d0 + b0*b0*c2*d6 + b0*b2*c0*d2 + b0*b2*c0*d6 + b0*c0*c2*d3 + b0*c0*c2*d8 + b2*c0*c0*d8 + c0*c0*c2*d5)/3)/delta
,
this->t*0.125f*(-(b0*b0*b1*d2 + b0*b0*c1*d8 + b0*b1*c0*d1 + b0*b1*c0*d8 + b0*c0*c1*d5 + b0*c0*c1*d7 + b1*c0*c0*d7 + c0*c0*c1*d4)/3 + (b0*b0*b2*d2 + b0*b0*c2*d8 + b0*b2*c0*d1 + b0*b2*c0*d8 + b0*c0*c2*d5 + b0*c0*c2*d7 + b2*c0*c0*d7 + c0*c0*c2*d4)/3)/delta
,
this->t*((0.0625f*b0*b0*b1*b1*d0 + 0.0625f*b0*b0*b1*c1*d2 + 0.0625f*b0*b0*b1*c1*d6 + 0.0625f*b0*b0*c1*c1*d8 + 0.0625f*b0*b1*b1*c0*d2 + 0.0625f*b0*b1*b1*c0*d6 + 0.0625f*b0*b1*c0*c1*d1 + 0.0625f*b0*b1*c0*c1*d3 + 0.125f*b0*b1*c0*c1*d8 + 0.0625f*b0*c0*c1*c1*d5 + 0.0625f*b0*c0*c1*c1*d7 + 0.0625f*b1*b1*c0*c0*d8 + 0.0625f*b1*c0*c0*c1*d5 + 0.0625f*b1*c0*c0*c1*d7 + 0.0625f*c0*c0*c1*c1*d4)/6 + (0.0625f*b0*b0*b2*b2*d0 + 0.0625f*b0*b0*b2*c2*d2 + 0.0625f*b0*b0*b2*c2*d6 + 0.0625f*b0*b0*c2*c2*d8 + 0.0625f*b0*b2*b2*c0*d2 + 0.0625f*b0*b2*b2*c0*d6 + 0.0625f*b0*b2*c0*c2*d1 + 0.0625f*b0*b2*c0*c2*d3 + 0.125f*b0*b2*c0*c2*d8 + 0.0625f*b0*c0*c2*c2*d5 + 0.0625f*b0*c0*c2*c2*d7 + 0.0625f*b2*b2*c0*c0*d8 + 0.0625f*b2*c0*c0*c2*d5 + 0.0625f*b2*c0*c0*c2*d7 + 0.0625f*c0*c0*c2*c2*d4)/6 - (0.125f*b0*b0*b1*b2*d0 + 0.0625f*b0*b0*b1*c2*d2 + 0.0625f*b0*b0*b1*c2*d6 + 0.0625f*b0*b0*b2*c1*d2 + 0.0625f*b0*b0*b2*c1*d6 + 0.125f*b0*b0*c1*c2*d8 + 0.125f*b0*b1*b2*c0*d2 + 0.125f*b0*b1*b2*c0*d6 + 0.0625f*b0*b1*c0*c2*d1 + 0.0625f*b0*b1*c0*c2*d3 + 0.125f*b0*b1*c0*c2*d8 + 0.0625f*b0*b2*c0*c1*d1 + 0.0625f*b0*b2*c0*c1*d3 + 0.125f*b0*b2*c0*c1*d8 + 0.125f*b0*c0*c1*c2*d5 + 0.125f*b0*c0*c1*c2*d7 + 0.125f*b1*b2*c0*c0*d8 + 0.0625f*b1*c0*c0*c2*d5 + 0.0625f*b1*c0*c0*c2*d7 + 0.0625f*b2*c0*c0*c1*d5 + 0.0625f*b2*c0*c0*c1*d7 + 0.125f*c0*c0*c1*c2*d4)/12)/delta
,
this->t*0.125f*(-(b0*b1*b1*d0 + b0*b1*c1*d2 + b0*b1*c1*d6 + b0*c1*c1*d8 + b1*b1*c0*d6 + b1*c0*c1*d3 + b1*c0*c1*d8 + c0*c1*c1*d5)/3 + (b0*b1*b2*d0 + b0*b1*c2*d6 + b0*b2*c1*d2 + b0*c1*c2*d8 + b1*b2*c0*d6 + b1*c0*c2*d3 + b2*c0*c1*d8 + c0*c1*c2*d5)/3)/delta
,
this->t*0.125f*(-(b0*b1*b1*d2 + b0*b1*c1*d1 + b0*b1*c1*d8 + b0*c1*c1*d7 + b1*b1*c0*d8 + b1*c0*c1*d5 + b1*c0*c1*d7 + c0*c1*c1*d4)/3 + (b0*b1*b2*d2 + b0*b1*c2*d8 + b0*b2*c1*d1 + b0*c1*c2*d7 + b1*b2*c0*d8 + b1*c0*c2*d5 + b2*c0*c1*d7 + c0*c1*c2*d4)/3)/delta
,
this->t*(-(0.0625f*b0*b0*b1*b1*d0 + 0.0625f*b0*b0*b1*c1*d2 + 0.0625f*b0*b0*b1*c1*d6 + 0.0625f*b0*b0*c1*c1*d8 + 0.0625f*b0*b1*b1*c0*d2 + 0.0625f*b0*b1*b1*c0*d6 + 0.0625f*b0*b1*c0*c1*d1 + 0.0625f*b0*b1*c0*c1*d3 + 0.125f*b0*b1*c0*c1*d8 + 0.0625f*b0*c0*c1*c1*d5 + 0.0625f*b0*c0*c1*c1*d7 + 0.0625f*b1*b1*c0*c0*d8 + 0.0625f*b1*c0*c0*c1*d5 + 0.0625f*b1*c0*c0*c1*d7 + 0.0625f*c0*c0*c1*c1*d4)/6 - 0.0625f*(b0*b1*b2*b2*d0 + b0*b1*b2*c2*d2 + b0*b1*b2*c2*d6 + b0*b1*c2*c2*d8 + b0*b2*b2*c1*d2 + b0*b2*c1*c2*d1 + b0*b2*c1*c2*d8 + b0*c1*c2*c2*d7 + b1*b2*b2*c0*d6 + b1*b2*c0*c2*d3 + b1*b2*c0*c2*d8 + b1*c0*c2*c2*d5 + b2*b2*c0*c1*d8 + b2*c0*c1*c2*d5 + b2*c0*c1*c2*d7 + c0*c1*c2*c2*d4)/12 + 0.0625f*(b0*b1*b1*b2*d0 + b0*b1*b1*c2*d2 + b0*b1*b2*c1*d2 + b0*b1*b2*c1*d6 + b0*b1*c1*c2*d1 + b0*b1*c1*c2*d8 + b0*b2*c1*c1*d8 + b0*c1*c1*c2*d7 + b1*b1*b2*c0*d6 + b1*b1*c0*c2*d8 + b1*b2*c0*c1*d3 + b1*b2*c0*c1*d8 + b1*c0*c1*c2*d5 + b1*c0*c1*c2*d7 + b2*c0*c1*c1*d5 + c0*c1*c1*c2*d4)/12 + 0.0625f*(b0*b0*b1*b2*d0 + b0*b0*b1*c2*d6 + b0*b0*b2*c1*d2 + b0*b0*c1*c2*d8 + b0*b1*b2*c0*d2 + b0*b1*b2*c0*d6 + b0*b1*c0*c2*d3 + b0*b1*c0*c2*d8 + b0*b2*c0*c1*d1 + b0*b2*c0*c1*d8 + b0*c0*c1*c2*d5 + b0*c0*c1*c2*d7 + b1*b2*c0*c0*d8 + b1*c0*c0*c2*d5 + b2*c0*c0*c1*d7 + c0*c0*c1*c2*d4)/12)/delta
,
this->t*0.125f*((b0*b2*b2*d0 + b0*b2*c2*d2 + b0*b2*c2*d6 + b0*c2*c2*d8 + b2*b2*c0*d6 + b2*c0*c2*d3 + b2*c0*c2*d8 + c0*c2*c2*d5)/3 - (b0*b1*b2*d0 + b0*b1*c2*d2 + b0*b2*c1*d6 + b0*c1*c2*d8 + b1*b2*c0*d6 + b1*c0*c2*d8 + b2*c0*c1*d3 + c0*c1*c2*d5)/3)/delta
,
this->t*0.125f*((b0*b2*b2*d2 + b0*b2*c2*d1 + b0*b2*c2*d8 + b0*c2*c2*d7 + b2*b2*c0*d8 + b2*c0*c2*d5 + b2*c0*c2*d7 + c0*c2*c2*d4)/3 - (b0*b1*b2*d2 + b0*b1*c2*d1 + b0*b2*c1*d8 + b0*c1*c2*d7 + b1*b2*c0*d8 + b1*c0*c2*d7 + b2*c0*c1*d5 + c0*c1*c2*d4)/3)/delta
,
this->t*(-(0.0625f*b0*b0*b2*b2*d0 + 0.0625f*b0*b0*b2*c2*d2 + 0.0625f*b0*b0*b2*c2*d6 + 0.0625f*b0*b0*c2*c2*d8 + 0.0625f*b0*b2*b2*c0*d2 + 0.0625f*b0*b2*b2*c0*d6 + 0.0625f*b0*b2*c0*c2*d1 + 0.0625f*b0*b2*c0*c2*d3 + 0.125f*b0*b2*c0*c2*d8 + 0.0625f*b0*c0*c2*c2*d5 + 0.0625f*b0*c0*c2*c2*d7 + 0.0625f*b2*b2*c0*c0*d8 + 0.0625f*b2*c0*c0*c2*d5 + 0.0625f*b2*c0*c0*c2*d7 + 0.0625f*c0*c0*c2*c2*d4)/6 + 0.0625f*(b0*b1*b2*b2*d0 + b0*b1*b2*c2*d2 + b0*b1*b2*c2*d6 + b0*b1*c2*c2*d8 + b0*b2*b2*c1*d2 + b0*b2*c1*c2*d1 + b0*b2*c1*c2*d8 + b0*c1*c2*c2*d7 + b1*b2*b2*c0*d6 + b1*b2*c0*c2*d3 + b1*b2*c0*c2*d8 + b1*c0*c2*c2*d5 + b2*b2*c0*c1*d8 + b2*c0*c1*c2*d5 + b2*c0*c1*c2*d7 + c0*c1*c2*c2*d4)/12 - 0.0625f*(b0*b1*b1*b2*d0 + b0*b1*b1*c2*d2 + b0*b1*b2*c1*d2 + b0*b1*b2*c1*d6 + b0*b1*c1*c2*d1 + b0*b1*c1*c2*d8 + b0*b2*c1*c1*d8 + b0*c1*c1*c2*d7 + b1*b1*b2*c0*d6 + b1*b1*c0*c2*d8 + b1*b2*c0*c1*d3 + b1*b2*c0*c1*d8 + b1*c0*c1*c2*d5 + b1*c0*c1*c2*d7 + b2*c0*c1*c1*d5 + c0*c1*c1*c2*d4)/12 + 0.0625f*(b0*b0*b1*b2*d0 + b0*b0*b1*c2*d2 + b0*b0*b2*c1*d6 + b0*b0*c1*c2*d8 + b0*b1*b2*c0*d2 + b0*b1*b2*c0*d6 + b0*b1*c0*c2*d1 + b0*b1*c0*c2*d8 + b0*b2*c0*c1*d3 + b0*b2*c0*c1*d8 + b0*c0*c1*c2*d5 + b0*c0*c1*c2*d7 + b1*b2*c0*c0*d8 + b1*c0*c0*c2*d7 + b2*c0*c0*c1*d5 + c0*c0*c1*c2*d4)/12)/delta
,
this->t*(b0*b1*d0 + b0*c1*d6 + b1*c0*d2 + c0*c1*d8)/(4*delta)
,
this->t*(b0*b1*d2 + b0*c1*d8 + b1*c0*d1 + c0*c1*d7)/(4*delta)
,
this->t*0.125f*(-(b0*b1*b1*d0 + b0*b1*c1*d2 + b0*b1*c1*d6 + b0*c1*c1*d8 + b1*b1*c0*d2 + b1*c0*c1*d1 + b1*c0*c1*d8 + c0*c1*c1*d7)/3 + (b0*b1*b2*d0 + b0*b1*c2*d2 + b0*b2*c1*d6 + b0*c1*c2*d8 + b1*b2*c0*d2 + b1*c0*c2*d1 + b2*c0*c1*d8 + c0*c1*c2*d7)/3)/delta
,
this->t*(b1*b1*d0 + b1*c1*d2 + b1*c1*d6 + c1*c1*d8)/(4*delta)
,
this->t*(b1*b1*d2 + b1*c1*d1 + b1*c1*d8 + c1*c1*d7)/(4*delta)
,
this->t*0.125f*((b0*b1*b1*d0 + b0*b1*c1*d2 + b0*b1*c1*d6 + b0*c1*c1*d8 + b1*b1*c0*d2 + b1*c0*c1*d1 + b1*c0*c1*d8 + c0*c1*c1*d7)/3 - (b1*b1*b2*d0 + b1*b1*c2*d2 + b1*b2*c1*d2 + b1*b2*c1*d6 + b1*c1*c2*d1 + b1*c1*c2*d8 + b2*c1*c1*d8 + c1*c1*c2*d7)/3)/delta
,
this->t*(b1*b2*d0 + b1*c2*d2 + b2*c1*d6 + c1*c2*d8)/(4*delta)
,
this->t*(b1*b2*d2 + b1*c2*d1 + b2*c1*d8 + c1*c2*d7)/(4*delta)
,
this->t*0.125f*((b1*b1*b2*d0 + b1*b1*c2*d2 + b1*b2*c1*d2 + b1*b2*c1*d6 + b1*c1*c2*d1 + b1*c1*c2*d8 + b2*c1*c1*d8 + c1*c1*c2*d7)/3 - (b0*b1*b2*d0 + b0*b1*c2*d2 + b0*b2*c1*d6 + b0*c1*c2*d8 + b1*b2*c0*d2 + b1*c0*c2*d1 + b2*c0*c1*d8 + c0*c1*c2*d7)/3)/delta
,
this->t*(b0*b1*d6 + b0*c1*d3 + b1*c0*d8 + c0*c1*d5)/(4*delta)
,
this->t*(b0*b1*d8 + b0*c1*d5 + b1*c0*d7 + c0*c1*d4)/(4*delta)
,
this->t*0.125f*(-(b0*b1*b1*d6 + b0*b1*c1*d3 + b0*b1*c1*d8 + b0*c1*c1*d5 + b1*b1*c0*d8 + b1*c0*c1*d5 + b1*c0*c1*d7 + c0*c1*c1*d4)/3 + (b0*b1*b2*d6 + b0*b1*c2*d8 + b0*b2*c1*d3 + b0*c1*c2*d5 + b1*b2*c0*d8 + b1*c0*c2*d7 + b2*c0*c1*d5 + c0*c1*c2*d4)/3)/delta
,
this->t*(b1*b1*d6 + b1*c1*d3 + b1*c1*d8 + c1*c1*d5)/(4*delta)
,
this->t*(b1*b1*d8 + b1*c1*d5 + b1*c1*d7 + c1*c1*d4)/(4*delta)
,
this->t*0.125f*((b0*b1*b1*d6 + b0*b1*c1*d3 + b0*b1*c1*d8 + b0*c1*c1*d5 + b1*b1*c0*d8 + b1*c0*c1*d5 + b1*c0*c1*d7 + c0*c1*c1*d4)/3 - (b1*b1*b2*d6 + b1*b1*c2*d8 + b1*b2*c1*d3 + b1*b2*c1*d8 + b1*c1*c2*d5 + b1*c1*c2*d7 + b2*c1*c1*d5 + c1*c1*c2*d4)/3)/delta
,
this->t*(b1*b2*d6 + b1*c2*d8 + b2*c1*d3 + c1*c2*d5)/(4*delta)
,
this->t*(b1*b2*d8 + b1*c2*d7 + b2*c1*d5 + c1*c2*d4)/(4*delta)
,
this->t*0.125f*((b1*b1*b2*d6 + b1*b1*c2*d8 + b1*b2*c1*d3 + b1*b2*c1*d8 + b1*c1*c2*d5 + b1*c1*c2*d7 + b2*c1*c1*d5 + c1*c1*c2*d4)/3 - (b0*b1*b2*d6 + b0*b1*c2*d8 + b0*b2*c1*d3 + b0*c1*c2*d5 + b1*b2*c0*d8 + b1*c0*c2*d7 + b2*c0*c1*d5 + c0*c1*c2*d4)/3)/delta
,
this->t*0.125f*((b0*b0*b1*d0 + b0*b0*c1*d6 + b0*b1*c0*d2 + b0*b1*c0*d6 + b0*c0*c1*d3 + b0*c0*c1*d8 + b1*c0*c0*d8 + c0*c0*c1*d5)/3 - (b0*b1*b2*d0 + b0*b1*c2*d6 + b0*b2*c1*d6 + b0*c1*c2*d3 + b1*b2*c0*d2 + b1*c0*c2*d8 + b2*c0*c1*d8 + c0*c1*c2*d5)/3)/delta
,
this->t*0.125f*((b0*b0*b1*d2 + b0*b0*c1*d8 + b0*b1*c0*d1 + b0*b1*c0*d8 + b0*c0*c1*d5 + b0*c0*c1*d7 + b1*c0*c0*d7 + c0*c0*c1*d4)/3 - (b0*b1*b2*d2 + b0*b1*c2*d8 + b0*b2*c1*d8 + b0*c1*c2*d5 + b1*b2*c0*d1 + b1*c0*c2*d7 + b2*c0*c1*d7 + c0*c1*c2*d4)/3)/delta
,
this->t*(-(0.0625f*b0*b0*b1*b1*d0 + 0.0625f*b0*b0*b1*c1*d2 + 0.0625f*b0*b0*b1*c1*d6 + 0.0625f*b0*b0*c1*c1*d8 + 0.0625f*b0*b1*b1*c0*d2 + 0.0625f*b0*b1*b1*c0*d6 + 0.0625f*b0*b1*c0*c1*d1 + 0.0625f*b0*b1*c0*c1*d3 + 0.125f*b0*b1*c0*c1*d8 + 0.0625f*b0*c0*c1*c1*d5 + 0.0625f*b0*c0*c1*c1*d7 + 0.0625f*b1*b1*c0*c0*d8 + 0.0625f*b1*c0*c0*c1*d5 + 0.0625f*b1*c0*c0*c1*d7 + 0.0625f*c0*c0*c1*c1*d4)/6 - 0.0625f*(b0*b1*b2*b2*d0 + b0*b1*b2*c2*d2 + b0*b1*b2*c2*d6 + b0*b1*c2*c2*d8 + b0*b2*b2*c1*d6 + b0*b2*c1*c2*d3 + b0*b2*c1*c2*d8 + b0*c1*c2*c2*d5 + b1*b2*b2*c0*d2 + b1*b2*c0*c2*d1 + b1*b2*c0*c2*d8 + b1*c0*c2*c2*d7 + b2*b2*c0*c1*d8 + b2*c0*c1*c2*d5 + b2*c0*c1*c2*d7 + c0*c1*c2*c2*d4)/12 + 0.0625f*(b0*b1*b1*b2*d0 + b0*b1*b1*c2*d6 + b0*b1*b2*c1*d2 + b0*b1*b2*c1*d6 + b0*b1*c1*c2*d3 + b0*b1*c1*c2*d8 + b0*b2*c1*c1*d8 + b0*c1*c1*c2*d5 + b1*b1*b2*c0*d2 + b1*b1*c0*c2*d8 + b1*b2*c0*c1*d1 + b1*b2*c0*c1*d8 + b1*c0*c1*c2*d5 + b1*c0*c1*c2*d7 + b2*c0*c1*c1*d7 + c0*c1*c1*c2*d4)/12 + 0.0625f*(b0*b0*b1*b2*d0 + b0*b0*b1*c2*d2 + b0*b0*b2*c1*d6 + b0*b0*c1*c2*d8 + b0*b1*b2*c0*d2 + b0*b1*b2*c0*d6 + b0*b1*c0*c2*d1 + b0*b1*c0*c2*d8 + b0*b2*c0*c1*d3 + b0*b2*c0*c1*d8 + b0*c0*c1*c2*d5 + b0*c0*c1*c2*d7 + b1*b2*c0*c0*d8 + b1*c0*c0*c2*d7 + b2*c0*c0*c1*d5 + c0*c0*c1*c2*d4)/12)/delta
,
this->t*0.125f*((b0*b1*b1*d0 + b0*b1*c1*d2 + b0*b1*c1*d6 + b0*c1*c1*d8 + b1*b1*c0*d6 + b1*c0*c1*d3 + b1*c0*c1*d8 + c0*c1*c1*d5)/3 - (b1*b1*b2*d0 + b1*b1*c2*d6 + b1*b2*c1*d2 + b1*b2*c1*d6 + b1*c1*c2*d3 + b1*c1*c2*d8 + b2*c1*c1*d8 + c1*c1*c2*d5)/3)/delta
,
this->t*0.125f*((b0*b1*b1*d2 + b0*b1*c1*d1 + b0*b1*c1*d8 + b0*c1*c1*d7 + b1*b1*c0*d8 + b1*c0*c1*d5 + b1*c0*c1*d7 + c0*c1*c1*d4)/3 - (b1*b1*b2*d2 + b1*b1*c2*d8 + b1*b2*c1*d1 + b1*b2*c1*d8 + b1*c1*c2*d5 + b1*c1*c2*d7 + b2*c1*c1*d7 + c1*c1*c2*d4)/3)/delta
,
this->t*((0.0625f*b0*b0*b1*b1*d0 + 0.0625f*b0*b0*b1*c1*d2 + 0.0625f*b0*b0*b1*c1*d6 + 0.0625f*b0*b0*c1*c1*d8 + 0.0625f*b0*b1*b1*c0*d2 + 0.0625f*b0*b1*b1*c0*d6 + 0.0625f*b0*b1*c0*c1*d1 + 0.0625f*b0*b1*c0*c1*d3 + 0.125f*b0*b1*c0*c1*d8 + 0.0625f*b0*c0*c1*c1*d5 + 0.0625f*b0*c0*c1*c1*d7 + 0.0625f*b1*b1*c0*c0*d8 + 0.0625f*b1*c0*c0*c1*d5 + 0.0625f*b1*c0*c0*c1*d7 + 0.0625f*c0*c0*c1*c1*d4)/6 + (0.0625f*b1*b1*b2*b2*d0 + 0.0625f*b1*b1*b2*c2*d2 + 0.0625f*b1*b1*b2*c2*d6 + 0.0625f*b1*b1*c2*c2*d8 + 0.0625f*b1*b2*b2*c1*d2 + 0.0625f*b1*b2*b2*c1*d6 + 0.0625f*b1*b2*c1*c2*d1 + 0.0625f*b1*b2*c1*c2*d3 + 0.125f*b1*b2*c1*c2*d8 + 0.0625f*b1*c1*c2*c2*d5 + 0.0625f*b1*c1*c2*c2*d7 + 0.0625f*b2*b2*c1*c1*d8 + 0.0625f*b2*c1*c1*c2*d5 + 0.0625f*b2*c1*c1*c2*d7 + 0.0625f*c1*c1*c2*c2*d4)/6 - (0.125f*b0*b1*b1*b2*d0 + 0.0625f*b0*b1*b1*c2*d2 + 0.0625f*b0*b1*b1*c2*d6 + 0.125f*b0*b1*b2*c1*d2 + 0.125f*b0*b1*b2*c1*d6 + 0.0625f*b0*b1*c1*c2*d1 + 0.0625f*b0*b1*c1*c2*d3 + 0.125f*b0*b1*c1*c2*d8 + 0.125f*b0*b2*c1*c1*d8 + 0.0625f*b0*c1*c1*c2*d5 + 0.0625f*b0*c1*c1*c2*d7 + 0.0625f*b1*b1*b2*c0*d2 + 0.0625f*b1*b1*b2*c0*d6 + 0.125f*b1*b1*c0*c2*d8 + 0.0625f*b1*b2*c0*c1*d1 + 0.0625f*b1*b2*c0*c1*d3 + 0.125f*b1*b2*c0*c1*d8 + 0.125f*b1*c0*c1*c2*d5 + 0.125f*b1*c0*c1*c2*d7 + 0.0625f*b2*c0*c1*c1*d5 + 0.0625f*b2*c0*c1*c1*d7 + 0.125f*c0*c1*c1*c2*d4)/12)/delta
,
this->t*0.125f*(-(b1*b2*b2*d0 + b1*b2*c2*d2 + b1*b2*c2*d6 + b1*c2*c2*d8 + b2*b2*c1*d6 + b2*c1*c2*d3 + b2*c1*c2*d8 + c1*c2*c2*d5)/3 + (b0*b1*b2*d0 + b0*b1*c2*d2 + b0*b2*c1*d6 + b0*c1*c2*d8 + b1*b2*c0*d6 + b1*c0*c2*d8 + b2*c0*c1*d3 + c0*c1*c2*d5)/3)/delta
,
this->t*0.125f*(-(b1*b2*b2*d2 + b1*b2*c2*d1 + b1*b2*c2*d8 + b1*c2*c2*d7 + b2*b2*c1*d8 + b2*c1*c2*d5 + b2*c1*c2*d7 + c1*c2*c2*d4)/3 + (b0*b1*b2*d2 + b0*b1*c2*d1 + b0*b2*c1*d8 + b0*c1*c2*d7 + b1*b2*c0*d8 + b1*c0*c2*d7 + b2*c0*c1*d5 + c0*c1*c2*d4)/3)/delta
,
this->t*(-(0.0625f*b1*b1*b2*b2*d0 + 0.0625f*b1*b1*b2*c2*d2 + 0.0625f*b1*b1*b2*c2*d6 + 0.0625f*b1*b1*c2*c2*d8 + 0.0625f*b1*b2*b2*c1*d2 + 0.0625f*b1*b2*b2*c1*d6 + 0.0625f*b1*b2*c1*c2*d1 + 0.0625f*b1*b2*c1*c2*d3 + 0.125f*b1*b2*c1*c2*d8 + 0.0625f*b1*c1*c2*c2*d5 + 0.0625f*b1*c1*c2*c2*d7 + 0.0625f*b2*b2*c1*c1*d8 + 0.0625f*b2*c1*c1*c2*d5 + 0.0625f*b2*c1*c1*c2*d7 + 0.0625f*c1*c1*c2*c2*d4)/6 + 0.0625f*(b0*b1*b2*b2*d0 + b0*b1*b2*c2*d2 + b0*b1*b2*c2*d6 + b0*b1*c2*c2*d8 + b0*b2*b2*c1*d6 + b0*b2*c1*c2*d3 + b0*b2*c1*c2*d8 + b0*c1*c2*c2*d5 + b1*b2*b2*c0*d2 + b1*b2*c0*c2*d1 + b1*b2*c0*c2*d8 + b1*c0*c2*c2*d7 + b2*b2*c0*c1*d8 + b2*c0*c1*c2*d5 + b2*c0*c1*c2*d7 + c0*c1*c2*c2*d4)/12 + 0.0625f*(b0*b1*b1*b2*d0 + b0*b1*b1*c2*d2 + b0*b1*b2*c1*d2 + b0*b1*b2*c1*d6 + b0*b1*c1*c2*d1 + b0*b1*c1*c2*d8 + b0*b2*c1*c1*d8 + b0*c1*c1*c2*d7 + b1*b1*b2*c0*d6 + b1*b1*c0*c2*d8 + b1*b2*c0*c1*d3 + b1*b2*c0*c1*d8 + b1*c0*c1*c2*d5 + b1*c0*c1*c2*d7 + b2*c0*c1*c1*d5 + c0*c1*c1*c2*d4)/12 - 0.0625f*(b0*b0*b1*b2*d0 + b0*b0*b1*c2*d2 + b0*b0*b2*c1*d6 + b0*b0*c1*c2*d8 + b0*b1*b2*c0*d2 + b0*b1*b2*c0*d6 + b0*b1*c0*c2*d1 + b0*b1*c0*c2*d8 + b0*b2*c0*c1*d3 + b0*b2*c0*c1*d8 + b0*c0*c1*c2*d5 + b0*c0*c1*c2*d7 + b1*b2*c0*c0*d8 + b1*c0*c0*c2*d7 + b2*c0*c0*c1*d5 + c0*c0*c1*c2*d4)/12)/delta
,
this->t*(b0*b2*d0 + b0*c2*d6 + b2*c0*d2 + c0*c2*d8)/(4*delta)
,
this->t*(b0*b2*d2 + b0*c2*d8 + b2*c0*d1 + c0*c2*d7)/(4*delta)
,
this->t*0.125f*((b0*b2*b2*d0 + b0*b2*c2*d2 + b0*b2*c2*d6 + b0*c2*c2*d8 + b2*b2*c0*d2 + b2*c0*c2*d1 + b2*c0*c2*d8 + c0*c2*c2*d7)/3 - (b0*b1*b2*d0 + b0*b1*c2*d6 + b0*b2*c1*d2 + b0*c1*c2*d8 + b1*b2*c0*d2 + b1*c0*c2*d8 + b2*c0*c1*d1 + c0*c1*c2*d7)/3)/delta
,
this->t*(b1*b2*d0 + b1*c2*d6 + b2*c1*d2 + c1*c2*d8)/(4*delta)
,
this->t*(b1*b2*d2 + b1*c2*d8 + b2*c1*d1 + c1*c2*d7)/(4*delta)
,
this->t*0.125f*(-(b1*b2*b2*d0 + b1*b2*c2*d2 + b1*b2*c2*d6 + b1*c2*c2*d8 + b2*b2*c1*d2 + b2*c1*c2*d1 + b2*c1*c2*d8 + c1*c2*c2*d7)/3 + (b0*b1*b2*d0 + b0*b1*c2*d6 + b0*b2*c1*d2 + b0*c1*c2*d8 + b1*b2*c0*d2 + b1*c0*c2*d8 + b2*c0*c1*d1 + c0*c1*c2*d7)/3)/delta
,
this->t*(b2*b2*d0 + b2*c2*d2 + b2*c2*d6 + c2*c2*d8)/(4*delta)
,
this->t*(b2*b2*d2 + b2*c2*d1 + b2*c2*d8 + c2*c2*d7)/(4*delta)
,
this->t*0.125f*(-(b0*b2*b2*d0 + b0*b2*c2*d2 + b0*b2*c2*d6 + b0*c2*c2*d8 + b2*b2*c0*d2 + b2*c0*c2*d1 + b2*c0*c2*d8 + c0*c2*c2*d7)/3 + (b1*b2*b2*d0 + b1*b2*c2*d2 + b1*b2*c2*d6 + b1*c2*c2*d8 + b2*b2*c1*d2 + b2*c1*c2*d1 + b2*c1*c2*d8 + c1*c2*c2*d7)/3)/delta
,
this->t*(b0*b2*d6 + b0*c2*d3 + b2*c0*d8 + c0*c2*d5)/(4*delta)
,
this->t*(b0*b2*d8 + b0*c2*d5 + b2*c0*d7 + c0*c2*d4)/(4*delta)
,
this->t*0.125f*((b0*b2*b2*d6 + b0*b2*c2*d3 + b0*b2*c2*d8 + b0*c2*c2*d5 + b2*b2*c0*d8 + b2*c0*c2*d5 + b2*c0*c2*d7 + c0*c2*c2*d4)/3 - (b0*b1*b2*d6 + b0*b1*c2*d3 + b0*b2*c1*d8 + b0*c1*c2*d5 + b1*b2*c0*d8 + b1*c0*c2*d5 + b2*c0*c1*d7 + c0*c1*c2*d4)/3)/delta
,
this->t*(b1*b2*d6 + b1*c2*d3 + b2*c1*d8 + c1*c2*d5)/(4*delta)
,
this->t*(b1*b2*d8 + b1*c2*d5 + b2*c1*d7 + c1*c2*d4)/(4*delta)
,
this->t*0.125f*(-(b1*b2*b2*d6 + b1*b2*c2*d3 + b1*b2*c2*d8 + b1*c2*c2*d5 + b2*b2*c1*d8 + b2*c1*c2*d5 + b2*c1*c2*d7 + c1*c2*c2*d4)/3 + (b0*b1*b2*d6 + b0*b1*c2*d3 + b0*b2*c1*d8 + b0*c1*c2*d5 + b1*b2*c0*d8 + b1*c0*c2*d5 + b2*c0*c1*d7 + c0*c1*c2*d4)/3)/delta
,
this->t*(b2*b2*d6 + b2*c2*d3 + b2*c2*d8 + c2*c2*d5)/(4*delta)
,
this->t*(b2*b2*d8 + b2*c2*d5 + b2*c2*d7 + c2*c2*d4)/(4*delta)
,
this->t*0.125f*(-(b0*b2*b2*d6 + b0*b2*c2*d3 + b0*b2*c2*d8 + b0*c2*c2*d5 + b2*b2*c0*d8 + b2*c0*c2*d5 + b2*c0*c2*d7 + c0*c2*c2*d4)/3 + (b1*b2*b2*d6 + b1*b2*c2*d3 + b1*b2*c2*d8 + b1*c2*c2*d5 + b2*b2*c1*d8 + b2*c1*c2*d5 + b2*c1*c2*d7 + c1*c2*c2*d4)/3)/delta
,
this->t*0.125f*(-(b0*b0*b2*d0 + b0*b0*c2*d6 + b0*b2*c0*d2 + b0*b2*c0*d6 + b0*c0*c2*d3 + b0*c0*c2*d8 + b2*c0*c0*d8 + c0*c0*c2*d5)/3 + (b0*b1*b2*d0 + b0*b1*c2*d6 + b0*b2*c1*d6 + b0*c1*c2*d3 + b1*b2*c0*d2 + b1*c0*c2*d8 + b2*c0*c1*d8 + c0*c1*c2*d5)/3)/delta
,
this->t*0.125f*(-(b0*b0*b2*d2 + b0*b0*c2*d8 + b0*b2*c0*d1 + b0*b2*c0*d8 + b0*c0*c2*d5 + b0*c0*c2*d7 + b2*c0*c0*d7 + c0*c0*c2*d4)/3 + (b0*b1*b2*d2 + b0*b1*c2*d8 + b0*b2*c1*d8 + b0*c1*c2*d5 + b1*b2*c0*d1 + b1*c0*c2*d7 + b2*c0*c1*d7 + c0*c1*c2*d4)/3)/delta
,
this->t*(-(0.0625f*b0*b0*b2*b2*d0 + 0.0625f*b0*b0*b2*c2*d2 + 0.0625f*b0*b0*b2*c2*d6 + 0.0625f*b0*b0*c2*c2*d8 + 0.0625f*b0*b2*b2*c0*d2 + 0.0625f*b0*b2*b2*c0*d6 + 0.0625f*b0*b2*c0*c2*d1 + 0.0625f*b0*b2*c0*c2*d3 + 0.125f*b0*b2*c0*c2*d8 + 0.0625f*b0*c0*c2*c2*d5 + 0.0625f*b0*c0*c2*c2*d7 + 0.0625f*b2*b2*c0*c0*d8 + 0.0625f*b2*c0*c0*c2*d5 + 0.0625f*b2*c0*c0*c2*d7 + 0.0625f*c0*c0*c2*c2*d4)/6 + 0.0625f*(b0*b1*b2*b2*d0 + b0*b1*b2*c2*d2 + b0*b1*b2*c2*d6 + b0*b1*c2*c2*d8 + b0*b2*b2*c1*d6 + b0*b2*c1*c2*d3 + b0*b2*c1*c2*d8 + b0*c1*c2*c2*d5 + b1*b2*b2*c0*d2 + b1*b2*c0*c2*d1 + b1*b2*c0*c2*d8 + b1*c0*c2*c2*d7 + b2*b2*c0*c1*d8 + b2*c0*c1*c2*d5 + b2*c0*c1*c2*d7 + c0*c1*c2*c2*d4)/12 - 0.0625f*(b0*b1*b1*b2*d0 + b0*b1*b1*c2*d6 + b0*b1*b2*c1*d2 + b0*b1*b2*c1*d6 + b0*b1*c1*c2*d3 + b0*b1*c1*c2*d8 + b0*b2*c1*c1*d8 + b0*c1*c1*c2*d5 + b1*b1*b2*c0*d2 + b1*b1*c0*c2*d8 + b1*b2*c0*c1*d1 + b1*b2*c0*c1*d8 + b1*c0*c1*c2*d5 + b1*c0*c1*c2*d7 + b2*c0*c1*c1*d7 + c0*c1*c1*c2*d4)/12 + 0.0625f*(b0*b0*b1*b2*d0 + b0*b0*b1*c2*d6 + b0*b0*b2*c1*d2 + b0*b0*c1*c2*d8 + b0*b1*b2*c0*d2 + b0*b1*b2*c0*d6 + b0*b1*c0*c2*d3 + b0*b1*c0*c2*d8 + b0*b2*c0*c1*d1 + b0*b2*c0*c1*d8 + b0*c0*c1*c2*d5 + b0*c0*c1*c2*d7 + b1*b2*c0*c0*d8 + b1*c0*c0*c2*d5 + b2*c0*c0*c1*d7 + c0*c0*c1*c2*d4)/12)/delta
,
this->t*0.125f*((b1*b1*b2*d0 + b1*b1*c2*d6 + b1*b2*c1*d2 + b1*b2*c1*d6 + b1*c1*c2*d3 + b1*c1*c2*d8 + b2*c1*c1*d8 + c1*c1*c2*d5)/3 - (b0*b1*b2*d0 + b0*b1*c2*d6 + b0*b2*c1*d2 + b0*c1*c2*d8 + b1*b2*c0*d6 + b1*c0*c2*d3 + b2*c0*c1*d8 + c0*c1*c2*d5)/3)/delta
,
this->t*0.125f*((b1*b1*b2*d2 + b1*b1*c2*d8 + b1*b2*c1*d1 + b1*b2*c1*d8 + b1*c1*c2*d5 + b1*c1*c2*d7 + b2*c1*c1*d7 + c1*c1*c2*d4)/3 - (b0*b1*b2*d2 + b0*b1*c2*d8 + b0*b2*c1*d1 + b0*c1*c2*d7 + b1*b2*c0*d8 + b1*c0*c2*d5 + b2*c0*c1*d7 + c0*c1*c2*d4)/3)/delta
,
this->t*(-(0.0625f*b1*b1*b2*b2*d0 + 0.0625f*b1*b1*b2*c2*d2 + 0.0625f*b1*b1*b2*c2*d6 + 0.0625f*b1*b1*c2*c2*d8 + 0.0625f*b1*b2*b2*c1*d2 + 0.0625f*b1*b2*b2*c1*d6 + 0.0625f*b1*b2*c1*c2*d1 + 0.0625f*b1*b2*c1*c2*d3 + 0.125f*b1*b2*c1*c2*d8 + 0.0625f*b1*c1*c2*c2*d5 + 0.0625f*b1*c1*c2*c2*d7 + 0.0625f*b2*b2*c1*c1*d8 + 0.0625f*b2*c1*c1*c2*d5 + 0.0625f*b2*c1*c1*c2*d7 + 0.0625f*c1*c1*c2*c2*d4)/6 + 0.0625f*(b0*b1*b2*b2*d0 + b0*b1*b2*c2*d2 + b0*b1*b2*c2*d6 + b0*b1*c2*c2*d8 + b0*b2*b2*c1*d2 + b0*b2*c1*c2*d1 + b0*b2*c1*c2*d8 + b0*c1*c2*c2*d7 + b1*b2*b2*c0*d6 + b1*b2*c0*c2*d3 + b1*b2*c0*c2*d8 + b1*c0*c2*c2*d5 + b2*b2*c0*c1*d8 + b2*c0*c1*c2*d5 + b2*c0*c1*c2*d7 + c0*c1*c2*c2*d4)/12 + 0.0625f*(b0*b1*b1*b2*d0 + b0*b1*b1*c2*d6 + b0*b1*b2*c1*d2 + b0*b1*b2*c1*d6 + b0*b1*c1*c2*d3 + b0*b1*c1*c2*d8 + b0*b2*c1*c1*d8 + b0*c1*c1*c2*d5 + b1*b1*b2*c0*d2 + b1*b1*c0*c2*d8 + b1*b2*c0*c1*d1 + b1*b2*c0*c1*d8 + b1*c0*c1*c2*d5 + b1*c0*c1*c2*d7 + b2*c0*c1*c1*d7 + c0*c1*c1*c2*d4)/12 - 0.0625f*(b0*b0*b1*b2*d0 + b0*b0*b1*c2*d6 + b0*b0*b2*c1*d2 + b0*b0*c1*c2*d8 + b0*b1*b2*c0*d2 + b0*b1*b2*c0*d6 + b0*b1*c0*c2*d3 + b0*b1*c0*c2*d8 + b0*b2*c0*c1*d1 + b0*b2*c0*c1*d8 + b0*c0*c1*c2*d5 + b0*c0*c1*c2*d7 + b1*b2*c0*c0*d8 + b1*c0*c0*c2*d5 + b2*c0*c0*c1*d7 + c0*c0*c1*c2*d4)/12)/delta
,
this->t*0.125f*(-(b0*b2*b2*d0 + b0*b2*c2*d2 + b0*b2*c2*d6 + b0*c2*c2*d8 + b2*b2*c0*d6 + b2*c0*c2*d3 + b2*c0*c2*d8 + c0*c2*c2*d5)/3 + (b1*b2*b2*d0 + b1*b2*c2*d2 + b1*b2*c2*d6 + b1*c2*c2*d8 + b2*b2*c1*d6 + b2*c1*c2*d3 + b2*c1*c2*d8 + c1*c2*c2*d5)/3)/delta
,
this->t*0.125f*(-(b0*b2*b2*d2 + b0*b2*c2*d1 + b0*b2*c2*d8 + b0*c2*c2*d7 + b2*b2*c0*d8 + b2*c0*c2*d5 + b2*c0*c2*d7 + c0*c2*c2*d4)/3 + (b1*b2*b2*d2 + b1*b2*c2*d1 + b1*b2*c2*d8 + b1*c2*c2*d7 + b2*b2*c1*d8 + b2*c1*c2*d5 + b2*c1*c2*d7 + c1*c2*c2*d4)/3)/delta
,
this->t*((0.0625f*b0*b0*b2*b2*d0 + 0.0625f*b0*b0*b2*c2*d2 + 0.0625f*b0*b0*b2*c2*d6 + 0.0625f*b0*b0*c2*c2*d8 + 0.0625f*b0*b2*b2*c0*d2 + 0.0625f*b0*b2*b2*c0*d6 + 0.0625f*b0*b2*c0*c2*d1 + 0.0625f*b0*b2*c0*c2*d3 + 0.125f*b0*b2*c0*c2*d8 + 0.0625f*b0*c0*c2*c2*d5 + 0.0625f*b0*c0*c2*c2*d7 + 0.0625f*b2*b2*c0*c0*d8 + 0.0625f*b2*c0*c0*c2*d5 + 0.0625f*b2*c0*c0*c2*d7 + 0.0625f*c0*c0*c2*c2*d4)/6 + (0.0625f*b1*b1*b2*b2*d0 + 0.0625f*b1*b1*b2*c2*d2 + 0.0625f*b1*b1*b2*c2*d6 + 0.0625f*b1*b1*c2*c2*d8 + 0.0625f*b1*b2*b2*c1*d2 + 0.0625f*b1*b2*b2*c1*d6 + 0.0625f*b1*b2*c1*c2*d1 + 0.0625f*b1*b2*c1*c2*d3 + 0.125f*b1*b2*c1*c2*d8 + 0.0625f*b1*c1*c2*c2*d5 + 0.0625f*b1*c1*c2*c2*d7 + 0.0625f*b2*b2*c1*c1*d8 + 0.0625f*b2*c1*c1*c2*d5 + 0.0625f*b2*c1*c1*c2*d7 + 0.0625f*c1*c1*c2*c2*d4)/6 - (0.125f*b0*b1*b2*b2*d0 + 0.125f*b0*b1*b2*c2*d2 + 0.125f*b0*b1*b2*c2*d6 + 0.125f*b0*b1*c2*c2*d8 + 0.0625f*b0*b2*b2*c1*d2 + 0.0625f*b0*b2*b2*c1*d6 + 0.0625f*b0*b2*c1*c2*d1 + 0.0625f*b0*b2*c1*c2*d3 + 0.125f*b0*b2*c1*c2*d8 + 0.0625f*b0*c1*c2*c2*d5 + 0.0625f*b0*c1*c2*c2*d7 + 0.0625f*b1*b2*b2*c0*d2 + 0.0625f*b1*b2*b2*c0*d6 + 0.0625f*b1*b2*c0*c2*d1 + 0.0625f*b1*b2*c0*c2*d3 + 0.125f*b1*b2*c0*c2*d8 + 0.0625f*b1*c0*c2*c2*d5 + 0.0625f*b1*c0*c2*c2*d7 + 0.125f*b2*b2*c0*c1*d8 + 0.125f*b2*c0*c1*c2*d5 + 0.125f*b2*c0*c1*c2*d7 + 0.125f*c0*c1*c2*c2*d4)/12)/delta
};

    return K;
}

MeshNode* GT9::get_stresses(size_t node, const std::vector<double>& u, double density) const{
    logger::log_assert(node >= 0 && node <= 3, logger::ERROR, "wrong value for BeamLinear2D node, must be either 0 or 1.");
    
    size_t N = this->nodes.size();

    std::vector<double> DB = this->get_DB(this->nodes[node]->point);

    MeshNode2D* n = static_cast<MeshNode2D*>(this->nodes[node]);
    for(size_t i = 0; i < 3; ++i){
        n->results[i] = 0;
        for(size_t l = 0; l < 3; ++l){
            for(size_t j = 0; j < 3; ++j){
                if(this->nodes[l]->u_pos[j] > -1){
                    n->results[i] += DB[3*N*i + 3*l + j]*u[this->nodes[l]->u_pos[j]];
                }
            }
        }
        n->results[i] = density*std::abs(n->results[i]);
    }

    return this->get_node(node);
}

double GT9::get_stress_at(gp_Pnt point, const std::vector<double>& u) const{
    size_t N = this->nodes.size();

    std::vector<double> DB = this->get_DB(point);

    std::vector<double> results(3, 0);
    for(size_t i = 0; i < 3; ++i){
        results[i] = 0;
        for(size_t l = 0; l < 3; ++l){
            for(size_t j = 0; j < 3; ++j){
                if(this->nodes[l]->u_pos[j] > -1){
                    results[i] += DB[3*N*i + 3*l + j]*u[this->nodes[l]->u_pos[j]];
                }
            }
        }
        results[i] = std::abs(results[i]);
    }

    return std::sqrt(std::pow(results[0], 2) - results[0]*results[1] + std::pow(results[1], 2) + 3*std::pow(results[2], 2));
}

MeshNode* GT9::get_internal_loads(size_t node, const std::vector<double>& u) const{
    logger::log_assert(node >= 0 && node <= 3, logger::ERROR, "wrong value for BeamLinear2D node, must be either 0 or 1.");

    std::vector<double> k = this->get_k();

    MeshNode2D* n = static_cast<MeshNode2D*>(this->nodes[node]);
    for(int i = 0; i < 3; ++i){
        n->results[i] = 0;
        for(size_t l = 0; l < 3; ++l){
            for(int j = 0; j < 3; ++j){
                if(this->nodes[l]->u_pos[j] > -1){
                    n->results[i] += k[node*3*3+l*3+j]*u[this->nodes[l]->u_pos[j]];
                }
            }
        }
        n->results[i] = std::abs(n->results[i]);
    }

    return this->get_node(node);
}

double GT9::get_compliance(const std::vector<double>& u, const std::vector<double>& l) const{
    auto k = this->get_k();
    std::vector<double> u_vec(9, 0);
    for(size_t i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            if(this->nodes[i]->u_pos[j] > -1){
                u_vec[i*3+j] = u[this->nodes[i]->u_pos[j]];
            }
        }
    }

    std::vector<double> f_vec(9, 0);

    if(l.size() > 0){
        std::vector<double> l_vec(9, 0);
        for(size_t k = 0; k < 3; ++k){
            for(int j = 0; j < 3; ++j){
                if(this->nodes[k]->u_pos[j] > -1){
                    l_vec[k*3+j] = l[this->nodes[k]->u_pos[j]];
                }
            }
        }
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 9, 1, 9, 1, k.data(), 9, u_vec.data(), 1, 0, f_vec.data(), 1);
        return cblas_ddot(9, l_vec.data(), 1, f_vec.data(), 1);
    } else {
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 9, 1, 9, 1, k.data(), 9, u_vec.data(), 1, 0, f_vec.data(), 1);
        return cblas_ddot(9, u_vec.data(), 1, f_vec.data(), 1);
    }


}

double GT9::get_volume() const{
    gp_Mat deltaM(1, this->nodes[0]->point.X(), this->nodes[0]->point.Y(),
                  1, this->nodes[1]->point.X(), this->nodes[1]->point.Y(),
                  1, this->nodes[2]->point.X(), this->nodes[2]->point.Y());

    return 0.5*std::abs(deltaM.Determinant())*this->t;
}

void GT9::get_virtual_load(double mult, gp_Pnt point, const std::vector<double>& u, std::vector<double>& l) const{
    std::vector<double> DB = this->get_DB(point);
    std::vector<double> V{1, -0.5, 0,
                         -0.5, 1, 0,
                         0,   0, 1.5};

    std::vector<double> u_vec(9, 0);
    for(size_t k = 0; k < 3; ++k){
        for(int j = 0; j < 3; ++j){
            if(this->nodes[k]->u_pos[j] > -1){
                u_vec[k*3+j] = u[this->nodes[k]->u_pos[j]];
            }
        }
    }

    std::vector<double> f_vec(9, 0);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 1, 9, 1, DB.data(), 9, u_vec.data(), 1, 0, f_vec.data(), 1);

    std::vector<double> res(f_vec);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 1, 3, 1, V.data(), 3, res.data(), 1, 0, f_vec.data(), 1);
    res = f_vec;

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 9, 1, 3, 1, DB.data(), 9, res.data(), 1, 0, f_vec.data(), 1);
    cblas_dscal(9, mult, f_vec.data(), 1);

    for(size_t k = 0; k < 3; ++k){
        for(int j = 0; j < 3; ++j){
            if(this->nodes[k]->u_pos[j] > -1){
                l[this->nodes[k]->u_pos[j]] += f_vec[k*3+j];
            }
        }
    }
}

TopoDS_Shape GT9::get_shape() const{
    TopoDS_Vertex v1 = BRepBuilderAPI_MakeVertex(this->nodes[0]->point);
    TopoDS_Vertex v2 = BRepBuilderAPI_MakeVertex(this->nodes[1]->point);
    TopoDS_Vertex v3 = BRepBuilderAPI_MakeVertex(this->nodes[2]->point);

    TopoDS_Edge e1 = BRepBuilderAPI_MakeEdge(v1, v2);
    TopoDS_Edge e2 = BRepBuilderAPI_MakeEdge(v2, v3);
    TopoDS_Edge e3 = BRepBuilderAPI_MakeEdge(v3, v1);

    TopoDS_Wire w = BRepBuilderAPI_MakeWire(e1, e2, e3);

    TopoDS_Face f = BRepBuilderAPI_MakeFace(w);

    return f;
}

gp_Pnt GT9::get_centroid() const{
    double x = 0;
    double y = 0;
    for(auto& n : this->nodes){
        x += n->point.X();
        y += n->point.Y();
    }

    return gp_Pnt(x/this->nodes.size(), y/this->nodes.size(), 0);
}


std::vector<double> GT9::get_DB(gp_Pnt point) const{
    size_t N = this->nodes.size();

    std::vector<gp_Pnt> p;
    for(auto n:this->nodes){
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

    double L0 = (a[0] + b[0]*point.X() + c[0]*point.Y())/(2*delta);
    double L1 = (a[1] + b[1]*point.X() + c[1]*point.Y())/(2*delta);
    double L2 = (a[2] + b[2]*point.X() + c[2]*point.Y())/(2*delta);

    double b0 = b[0];
    double b1 = b[1];
    double b2 = b[2];
    double c0 = c[0];
    double c1 = c[1];
    double c2 = c[2];
    
    auto D = this->mat->stiffness_2D();

    double d0 = D[0];
    double d1 = D[1];
    double d2 = D[2];
    double d3 = D[3];
    double d4 = D[4];
    double d5 = D[5];
    double d6 = D[6];
    double d7 = D[7];
    double d8 = D[8];

    std::vector<double> DB{
(b0*d0 + c0*d2)/(2*delta)
,
(b0*d2 + c0*d1)/(2*delta)
,
0.25f*(-(b0*b1*d0 + b0*c1*d2 + b1*c0*d2 + c0*c1*d1)*L2 + (b0*b2*d0 + b0*c2*d2 + b2*c0*d2 + c0*c2*d1)*L1)/delta
,
(b1*d0 + c1*d2)/(2*delta)
,
(b1*d2 + c1*d1)/(2*delta)
,
0.25f*((b0*b1*d0 + b0*c1*d2 + b1*c0*d2 + c0*c1*d1)*L2 - (b1*b2*d0 + b1*c2*d2 + b2*c1*d2 + c1*c2*d1)*L0)/delta
,
(b2*d0 + c2*d2)/(2*delta)
,
(b2*d2 + c2*d1)/(2*delta)
,
0.25f*(-(b0*b2*d0 + b0*c2*d2 + b2*c0*d2 + c0*c2*d1)*L1 + (b1*b2*d0 + b1*c2*d2 + b2*c1*d2 + c1*c2*d1)*L0)/delta
,
(b0*d3 + c0*d5)/(2*delta)
,
(b0*d5 + c0*d4)/(2*delta)
,
0.25f*(-(b0*b1*d3 + b0*c1*d5 + b1*c0*d5 + c0*c1*d4)*L2 + (b0*b2*d3 + b0*c2*d5 + b2*c0*d5 + c0*c2*d4)*L1)/delta
,
(b1*d3 + c1*d5)/(2*delta)
,
(b1*d5 + c1*d4)/(2*delta)
,
0.25f*((b0*b1*d3 + b0*c1*d5 + b1*c0*d5 + c0*c1*d4)*L2 - (b1*b2*d3 + b1*c2*d5 + b2*c1*d5 + c1*c2*d4)*L0)/delta
,
(b2*d3 + c2*d5)/(2*delta)
,
(b2*d5 + c2*d4)/(2*delta)
,
0.25f*(-(b0*b2*d3 + b0*c2*d5 + b2*c0*d5 + c0*c2*d4)*L1 + (b1*b2*d3 + b1*c2*d5 + b2*c1*d5 + c1*c2*d4)*L0)/delta
,
(b0*d6 + c0*d8)/(2*delta)
,
(b0*d8 + c0*d7)/(2*delta)
,
0.25f*(-(b0*b1*d6 + b0*c1*d8 + b1*c0*d8 + c0*c1*d7)*L2 + (b0*b2*d6 + b0*c2*d8 + b2*c0*d8 + c0*c2*d7)*L1)/delta
,
(b1*d6 + c1*d8)/(2*delta)
,
(b1*d8 + c1*d7)/(2*delta)
,
0.25f*((b0*b1*d6 + b0*c1*d8 + b1*c0*d8 + c0*c1*d7)*L2 - (b1*b2*d6 + b1*c2*d8 + b2*c1*d8 + c1*c2*d7)*L0)/delta
,
(b2*d6 + c2*d8)/(2*delta)
,
(b2*d8 + c2*d7)/(2*delta)
,
0.25f*(-(b0*b2*d6 + b0*c2*d8 + b2*c0*d8 + c0*c2*d7)*L1 + (b1*b2*d6 + b1*c2*d8 + b2*c1*d8 + c1*c2*d7)*L0)/delta
    };

    return DB;
}

}
