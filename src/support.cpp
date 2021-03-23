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
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Shell.hxx>
#include <BRep_Builder.hxx>
#include <STEPCAFControl_Writer.hxx>
#include <STEPControl_StepModelType.hxx>

Support::Support(bool X, bool Y, double thickness, std::vector<std::array<double, 2>> vertices){
    this->X = X;
    this->Y = Y;

    // THIS IS FAULTY. Need a better way to create the shape, otherwise
    // SolidClassifier and IntCurvesFace_ShapeIntersector won't work properly.
    // This workaround is sufficient for pathfinding but may cause problems
    // later on.
    // thickness = 0.002;
    
    TopoDS_Shell sh;
    BRep_Builder builder;
    builder.MakeShell(sh);

    for(size_t i = 1; i < vertices.size(); ++i){
        gp_Pnt p1(vertices[i-1][0], vertices[i-1][1], -thickness/2);
        gp_Pnt p2(vertices[i-1][0], vertices[i-1][1], thickness/2);
        gp_Pnt p3(vertices[i][0], vertices[i][1], thickness/2);
        gp_Pnt p4(vertices[i][0], vertices[i][1], -thickness/2);

        TopoDS_Edge e1 = BRepBuilderAPI_MakeEdge(p1, p2);
        TopoDS_Edge e2 = BRepBuilderAPI_MakeEdge(p2, p3);
        TopoDS_Edge e3 = BRepBuilderAPI_MakeEdge(p3, p4);
        TopoDS_Edge e4 = BRepBuilderAPI_MakeEdge(p4, p1);
        TopoDS_Wire w = BRepBuilderAPI_MakeWire(e1, e2, e3, e4);
        TopoDS_Face f = BRepBuilderAPI_MakeFace(w);
        builder.Add(sh, f);
    }

    this->shape = sh;
}

Support::Support(bool X, bool Y, bool Z, std::vector<std::array<double, 3>> vertices){
    // TODO
}


bool Support::is_inside(gp_Pnt p) const{
    BRepClass3d_SolidClassifier insider(this->shape);
    insider.Perform(p, 0.01);
    return insider.State() == TopAbs_ON;
}


double Support::get_distance(gp_Pnt p) const{
    TopoDS_Vertex v = BRepBuilderAPI_MakeVertex(p);

    BRepExtrema_DistShapeShape d(v, this->shape);
    d.Perform();
    return d.Value();
}
