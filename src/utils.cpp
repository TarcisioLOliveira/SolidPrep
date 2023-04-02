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

#include "utils.hpp"
#include "logger.hpp"
#include "Interface_Static.hxx"
#include <TopExp_Explorer.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Line.hxx>
#include <BRep_Tool.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <vector>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <TopTools_ListOfShape.hxx>
#include <gp_Circ.hxx>
#include <ShapeFix_Wireframe.hxx>
#include <ShapeFix_Shape.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <lapacke.h>
#include <STEPCAFControl_Reader.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BOPAlgo_BOP.hxx>

namespace utils{

void shape_to_file(const std::string& s, const TopoDS_Shape& t){
    STEPControl_Writer writer;
    STEPControl_StepModelType mode = STEPControl_AsIs;
    Interface_Static::SetIVal("write.surfacecurve.mode",0);
    IFSelect_ReturnStatus stat = writer.Transfer(t,mode);
    if(stat != IFSelect_RetDone){
        writer.PrintStatsTransfer(4, 0);
        std::cout << "ERROR: failed to translate model to STEP. Filename: " << s << std::endl;
        exit(EXIT_FAILURE);
    }
    IFSelect_ReturnStatus stat2 = writer.Write(s.c_str());
    if(stat2 != IFSelect_RetDone){
        writer.PrintStatsTransfer(4, 0);
        std::cout << "ERROR: failed to write translated model. Filename: " << s << std::endl;
        exit(EXIT_FAILURE);
    }
}


TopoDS_Shape sweep_surface(const std::vector<gp_Pnt>& spine, const TopoDS_Shape& surface, const TopoDS_Shape& base){
    TopoDS_Shape csec = BRepBuilderAPI_Copy(surface);
    TopoDS_Shape result = BRepBuilderAPI_Copy(base);
    gp_Vec init(spine[0], spine[1]);
    init.Scale(0.5);

    result = cut_shape(result, BRepPrimAPI_MakePrism(csec, init));
    gp_Pnt prev_mid = spine[0];
    prev_mid.BaryCenter(1, spine[1], 1);
    gp_Dir prev_dir(gp_Vec(spine[0], spine[1]));
    gp_Trsf trsf;
    trsf.SetTranslation(init);
    csec = BRepBuilderAPI_Transform(csec, trsf);
    trsf.SetTranslation(gp_Vec(0,0,0));
    for(size_t i = 2; i < spine.size(); ++i){
        gp_Pnt cur_mid = spine[i-1];
        cur_mid.BaryCenter(1, spine[i], 1);
        gp_Dir cur_dir(gp_Vec(spine[i-1], spine[i]));
        if(cur_dir.IsEqual(prev_dir, Precision::Confusion())){
            gp_Vec transl(prev_mid, cur_mid);
            TopoDS_Shape prism = BRepPrimAPI_MakePrism(csec, transl);
            result = cut_shape(result, prism);

            trsf.SetTranslation(transl);
            csec = BRepBuilderAPI_Transform(csec, trsf);
            trsf.SetTranslation(gp_Vec(0,0,0));
        } else {
            gp_Dir axis = prev_dir.Crossed(cur_dir);

            // Find center of rotation
            double d1 = axis.Dot(gp_Vec(gp_Pnt(0,0,0), spine[i-1]));
            double d2 = gp_Vec(prev_dir).Dot(gp_Vec(gp_Pnt(0,0,0), prev_mid));
            double d3 = gp_Vec(cur_dir).Dot(gp_Vec(gp_Pnt(0,0,0), cur_mid));

            gp_Mat M(axis.XYZ(), prev_dir.XYZ(), cur_dir.XYZ());
            gp_Pnt center(0,0,0);
            center.Translate((d1*prev_dir.Crossed(cur_dir)+d2*cur_dir.Crossed(axis)+d3*axis.Crossed(prev_dir))/M.Determinant());

            // std::vector<double> d({d1, d2, d3});
            // std::vector<double> M({axis.X(), axis.Y(), axis.Z(),
            //                        prev_dir.X(), prev_dir.Y(), prev_dir.Z(),
            //                        cur_dir.X(), cur_dir.Y(), cur_dir.Z()});
            // std::vector<int> ipiv(3);
            // LAPACKE_dgesv(LAPACK_ROW_MAJOR, 3, 1, M.data(), 3, ipiv.data(), d.data(), 1);
            // gp_Pnt center(d[0], d[1], d[2]);
            //
            // logger::quick_log(d1, d2, d3);
            // logger::quick_log(axis.X(), axis.Y(), axis.Z(), prev_dir.X(), prev_dir.Y(), cur_dir.X(), cur_dir.Y());
            // logger::quick_log(center.X(), center.Y(), prev_mid.X(), prev_mid.Y(), cur_mid.X(), cur_mid.Y());
            gp_Ax1 ax(center, axis);
            gp_Vec v1(center, cur_mid);
            gp_Vec v2(center, prev_mid);
            double ang = v2.AngleWithRef(v1, axis);
            TopoDS_Shape revol = BRepPrimAPI_MakeRevol(csec, ax, ang, true);
            result = cut_shape(result, revol);

            trsf.SetRotation(ax, ang);
            csec = BRepBuilderAPI_Transform(csec, trsf);
            trsf.SetRotation(ax, 0);
        }

        prev_mid = std::move(cur_mid);
        prev_dir = std::move(cur_dir);
    }
    gp_Vec fin(prev_mid, spine.back());
    TopoDS_Shape prism = BRepPrimAPI_MakePrism(csec, fin);
    result = cut_shape(result, prism);
    result = cut_shape(base, result);

    return result;
}
TopoDS_Shape cut_shape(const TopoDS_Shape& base, const TopoDS_Shape& cutter){
    TopTools_ListOfShape shapes1;
    shapes1.Append(base);
    TopTools_ListOfShape shapes2;
    shapes2.Append(cutter);
    BRepAlgoAPI_Cut cut;
    cut.SetTools(shapes2);
    cut.SetArguments(shapes1);
    cut.SetNonDestructive(false);
    cut.SimplifyResult();
    cut.SetRunParallel(true);
    cut.Build();
    //BRepAlgoAPI_Cut cut(base, cutter);
    return cut.Shape();
}

TopoDS_Shape fast_make_2D_beam(const std::vector<gp_Pnt>& spine, double diameter, const TopoDS_Shape& surface){
    TopoDS_Shape result = BRepBuilderAPI_Copy(surface);
    TopTools_ListOfShape shapes1;
    shapes1.Append(result);
    TopTools_ListOfShape shapes2;
    for(auto& n:spine){
        gp_Ax2 axis(n, gp_Dir(0,0,1));
        gp_Circ circ(axis, diameter/2);
        TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(circ);
        TopoDS_Wire wire = BRepBuilderAPI_MakeWire(edge);
        TopoDS_Face face = BRepBuilderAPI_MakeFace(wire);

        // Fastest method I managed to find, cut from original shape then
        // inverse cut back into original shape.
        // Still a bit slow though.
        shapes2.Clear();
        shapes2.Append(face);
        BRepAlgoAPI_Cut cut;
        cut.SetTools(shapes2);
        cut.SetArguments(shapes1);
        cut.SetNonDestructive(false);
        //cut.SimplifyResult();
        cut.SetRunParallel(true);
        cut.Build();
        result = cut;
        shapes1.Clear();
        shapes1.Append(result);
    }
    result = BRepAlgoAPI_Cut(surface, result);

    // Simplify
    double tol = diameter;
    double prec = tol*2;

    Handle(ShapeFix_Shape) sfs = new ShapeFix_Shape(result);
    sfs->SetPrecision(prec);
    sfs->SetMaxTolerance(2*tol);
    sfs->SetMinTolerance(tol);
    auto sfw = sfs->FixWireTool();
    sfw->ModifyGeometryMode() = true;
    sfw->ModifyTopologyMode() = true;
    sfw->FixSmallMode() = true;
    sfw->FixSmall(false, prec);
    sfs->Perform();
    result = sfs->Shape();

    Handle(ShapeFix_Wireframe) SFWF = new ShapeFix_Wireframe(result);
    SFWF->SetPrecision(prec);
    SFWF->SetMaxTolerance(2*tol);
    SFWF->SetMinTolerance(tol);
    SFWF->ModeDropSmallEdges() = Standard_True;
    SFWF->FixSmallEdges();
    SFWF->FixWireGaps();
    result = SFWF->Shape();

    return result;
}

TopoDS_Shape load_shape(const std::string& path, double scale) {
    TopoDS_Shape s;
    
    STEPControl_Reader reader;
    IFSelect_ReturnStatus stat = reader.ReadFile(path.c_str());
    if(stat != IFSelect_RetDone){
        reader.PrintCheckLoad(false, IFSelect_ItemsByEntity);
        exit(EXIT_FAILURE);
    }

    Standard_Integer NbRoots = reader.NbRootsForTransfer();
    Standard_Integer num = reader.TransferRoots();
    (void) NbRoots;
    (void) num;
    s = reader.OneShape();
    if(s.IsNull()){
        reader.PrintCheckTransfer(true, IFSelect_ItemsByEntity);
        exit(EXIT_FAILURE);
    }

    if(scale != 1){
        gp_Trsf t;
        t.SetScale(gp_Pnt(0, 0, 0), scale);
        BRepBuilderAPI_Transform transf(s, t, true);
        s = transf.Shape();
    }

    return s;
}

}
