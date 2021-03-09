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

#include "support.hpp"

#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <TopoDS_Wire.hxx>

Support::Support(bool X, bool Y, std::vector<std::array<double, 2>> vertices){
    this->X = X;
    this->Y = Y;

    gp_Pnt p(vertices[0][0], vertices[0][1], 0.0);
    gp_Pnt pi(vertices[1][0], vertices[1][1], 0.0);
    TopoDS_Edge e = BRepBuilderAPI_MakeEdge(p, pi);
    TopoDS_Wire w2 = BRepBuilderAPI_MakeWire(e);
    p = pi;

    for(auto i = vertices.begin()+2; i < vertices.end(); ++i){
        pi = gp_Pnt((*i)[0], (*i)[1], 0);
        e = BRepBuilderAPI_MakeEdge(p, pi);
        w2 = BRepBuilderAPI_MakeWire(w2, e);
        p = pi;
    }

    this->shape = w2;
}

Support::Support(bool X, bool Y, bool Z, std::vector<std::array<double, 3>> vertices){
    // TODO
}


bool Support::is_inside(gp_Pnt p) const{
    BRepClass3d_SolidClassifier insider(this->shape);
    insider.Perform(p, 0);
    return insider.State() == TopAbs_ON;
}


double Support::get_distance(gp_Pnt p) const{
    TopoDS_Vertex v = BRepBuilderAPI_MakeVertex(p);

    BRepExtrema_DistShapeShape d(v, this->shape);
    d.Perform();
    return d.Value();
}
