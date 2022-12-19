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

#include "cross_section.hpp"
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
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <TopoDS_Shell.hxx>
#include <BRep_Builder.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include "utils.hpp"
#include <gp_Circ.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>

CrossSection::CrossSection(std::vector<gp_Pnt> vertices, double thickness){
    size_t size = vertices.size();
    logger::log_assert(size > 1, logger::ERROR, "Bi-dimensional cross sections must have two or more vertices.");
    logger::log_assert(!utils::equal(vertices[0], vertices[size-1]), logger::ERROR, "The vertices of a bi-dimensional cross section must not form a closed polygon.");

    // Creates a 2D cross-section so that OCCT functions for line/face 
    // intersection can work correctly.
    gp_Pnt p1(vertices[0]);
    gp_Pnt p2(vertices[size-1]);
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

    Handle(Geom_Surface) surf = BRep_Tool::Surface(face);
    BOPTools_AlgoTools3D::GetNormalToSurface(surf, this->centroid.X(), this->centroid.Y(), this->normal);
    this->area = props.Mass();

    gp_Ax1 axis(this->centroid, gp_Dir(0, 0, 1));
    double ang = this->normal.AngleWithRef(gp_Dir(1,0,0), gp_Dir(0,0,1));
    gp_Trsf t;
    t.SetRotation(axis, ang);
    BRepBuilderAPI_Transform transf(face, t, true);
    face = TopoDS::Face(transf.Shape());
    BRepGProp::SurfaceProperties(face, props);

    this->inertia = props.MatrixOfInertia();

    // Shape creation.
    TopoDS_Shell sh;
    BRep_Builder builder;
    builder.MakeShell(sh);

    for(size_t i = 1; i < vertices.size(); ++i){
        gp_Pnt p1(vertices[i-1]);
        p1.SetZ(-thickness/2);
        gp_Pnt p2(vertices[i-1]);
        p2.SetZ(thickness/2);
        gp_Pnt p3(vertices[i]);
        p3.SetZ(thickness/2);
        gp_Pnt p4(vertices[i]);
        p4.SetZ(-thickness/2);

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
CrossSection::CrossSection(gp_Pnt p):
    centroid(p), inertia(), normal(), max_dim(0), shape(BRepBuilderAPI_MakeVertex(p)), area(0){

}

CrossSection::CrossSection(gp_Pnt p, utils::ProblemType type, double radius):
    centroid(p), inertia(), normal(), max_dim(radius), shape(), area(){
    if(type == utils::PROBLEM_TYPE_2D){
        // Creates a circle
        gp_Ax2 axis(p, gp_Dir(0,0,1));
        gp_Circ circ(axis, radius);
        TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(circ);
        TopoDS_Wire wire = BRepBuilderAPI_MakeWire(edge);
        TopoDS_Face face = BRepBuilderAPI_MakeFace(wire);
        this->shape = face;
    } else if(type == utils::PROBLEM_TYPE_3D){
        // Creates a sphere
        this->shape = BRepPrimAPI_MakeSphere(radius);
        gp_Trsf t;
        t.SetTranslation(gp_Pnt(0,0,0), p);
        BRepBuilderAPI_Transform transf(this->shape, t, true);
        this->shape = transf.Shape();
    } 
}

CrossSection::CrossSection(Rectangle r):
    centroid(r.center), inertia(), normal(r.normal), max_dim(std::max(r.w, r.h)), 
    shape(), area(r.w*r.h){
 
    gp_Ax1 ax1(r.center, r.normal.Crossed(gp_Dir(0,0,1))); 
    double ang1 = r.normal.AngleWithRef(gp_Dir(0,0,1), ax1.Direction());
    gp_Ax1 ax2(r.center, r.normal);
    gp_Pnt p1(r.center.X() - r.w/2, r.center.Y() - r.h/2, r.center.Z());
    gp_Pnt p2(r.center.X() + r.w/2, r.center.Y() - r.h/2, r.center.Z());
    gp_Pnt p3(r.center.X() + r.w/2, r.center.Y() + r.h/2, r.center.Z());
    gp_Pnt p4(r.center.X() - r.w/2, r.center.Y() + r.h/2, r.center.Z());
    p1.Rotate(ax1, ang1);
    p2.Rotate(ax1, ang1);
    p3.Rotate(ax1, ang1);
    p4.Rotate(ax1, ang1);
    p1.Rotate(ax2, r.rot_ang);
    p2.Rotate(ax2, r.rot_ang);
    p3.Rotate(ax2, r.rot_ang);
    p4.Rotate(ax2, r.rot_ang);

    const TopoDS_Vertex v1 = BRepBuilderAPI_MakeVertex(p1);
    const TopoDS_Vertex v2 = BRepBuilderAPI_MakeVertex(p2);
    const TopoDS_Vertex v3 = BRepBuilderAPI_MakeVertex(p3);
    const TopoDS_Vertex v4 = BRepBuilderAPI_MakeVertex(p4);

    const TopoDS_Edge e1 = BRepBuilderAPI_MakeEdge(v1, v2);
    const TopoDS_Edge e2 = BRepBuilderAPI_MakeEdge(v2, v3);
    const TopoDS_Edge e3 = BRepBuilderAPI_MakeEdge(v3, v4);
    const TopoDS_Edge e4 = BRepBuilderAPI_MakeEdge(v4, v1);

    const TopoDS_Wire w = BRepBuilderAPI_MakeWire(e1, e2, e3, e4);

    TopoDS_Face f = BRepBuilderAPI_MakeFace(w);

    GProp_GProps props;
    BRepGProp::SurfaceProperties(f, props);

    this->inertia = props.MatrixOfInertia();

    this->shape = std::move(f);
}

CrossSection::CrossSection(std::vector<gp_Pnt> vertices){
    (void)vertices;
    // TODO
}

CrossSection::CrossSection(const std::string& s):
    centroid(), inertia(), normal(), max_dim(),
    shape(utils::load_shape(s, 1.0)), area(){

    GProp_GProps props;
    BRepGProp::SurfaceProperties(this->shape, props);
    this->area = props.Mass();
    this->inertia = props.MatrixOfInertia();
    this->centroid = props.CentreOfMass();

    // TODO: this->normal, this->max_dim
    // maybe allow normal to be user-defined?
}

bool CrossSection::is_inside(gp_Pnt p) const{
    BRepClass3d_SolidClassifier insider(this->shape);
    insider.Perform(p, Precision::Confusion());
    return insider.State() == TopAbs_ON;
}

double CrossSection::get_distance(gp_Pnt p) const{
    TopoDS_Vertex v = BRepBuilderAPI_MakeVertex(p);

    BRepExtrema_DistShapeShape d(v, this->shape);
    d.Perform();
    return d.Value();
}

void CrossSection::set_centroid(gp_Pnt p){
    gp_Trsf t;
    t.SetTranslation(this->centroid, p);
    BRepBuilderAPI_Transform transf(this->shape, t, true);
    this->shape = transf.Shape();
    this->centroid = p;
}
