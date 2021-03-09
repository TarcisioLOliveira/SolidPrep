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

#include "force.hpp"
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRep_Tool.hxx>
#include <BOPTools_AlgoTools3D.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include "logger.hpp"


Force::Force(std::vector<std::array<double, 2>> vertices, double thickness, std::array<double, 2> force){
    this->force = gp_Vec(force[0], force[1], 0);

    // Calculation of properties
    size_t size = vertices.size();
    logger::log_assert(size > 1, logger::ERROR, "Forces must have two or more vertices.");
    logger::log_assert(vertices[0] != vertices[size-1], logger::ERROR, "Force's vertices must not form a closed polygon.");

    gp_Pnt p1(vertices[0][0], vertices[0][1], 0);
    gp_Pnt p2(vertices[size-1][0], vertices[size-1][1], 0);
    this->max_dim = p1.Distance(p2);
    this->centroid = p1;
    this->centroid.BaryCenter(1, p2, 1);

    gp_Pnt p0(p1.X(), p1.Y(), -thickness/2);
    p1.SetZ(thickness/2);
    p2.SetZ(thickness/2);
    gp_Pnt p3(p2.X(), p2.Y(), -thickness/2);

    TopoDS_Edge e0 = BRepBuilderAPI_MakeEdge(p0, p1);
    TopoDS_Edge e1 = BRepBuilderAPI_MakeEdge(p1, p2);
    TopoDS_Edge e2 = BRepBuilderAPI_MakeEdge(p2, p3);
    TopoDS_Edge e3 = BRepBuilderAPI_MakeEdge(p3, p0);
    TopoDS_Wire w = BRepBuilderAPI_MakeWire(e0, e1, e2, e3);
    TopoDS_Face face = BRepBuilderAPI_MakeFace(w);

    GProp_GProps props;
    BRepGProp::SurfaceProperties(face, props);

    this->inertia = props.MatrixOfInertia();
    gp_Pnt cm = props.CentreOfMass();

    Handle(Geom_Surface) surf = BRep_Tool::Surface(face);
    BOPTools_AlgoTools3D::GetNormalToSurface(surf, cm.X(), cm.Y(), this->normal);

    // Shape creation.
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


Force::Force(std::vector<std::array<double, 3>> vertices,  std::array<double, 3> force){
    // TODO
}

bool Force::is_inside(gp_Pnt p) const{
    BRepClass3d_SolidClassifier insider(this->shape);
    insider.Perform(p, 0);
    return insider.State() == TopAbs_ON;
}
