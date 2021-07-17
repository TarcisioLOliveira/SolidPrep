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

#include "sizing/beam_sizing.hpp"
#include "project_data.hpp"
#include "utils.hpp"
#include <lapacke.h>
#include <vector>
#include <cmath>
#include <cstring>
#include <logger.hpp>
#include <algorithm>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <gp_Circ.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Shape.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS.hxx>
#include <BOPAlgo_BOP.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Wireframe.hxx>

namespace sizing{

BeamSizing::BeamSizing(ProjectData* data, BeamElementFactory::BeamElementType t):
    Sizing(data), type(t){

}

TopoDS_Shape BeamSizing::run(){
    BeamGraph graph(this->data, this->type);
    graph.run();
    if(this->data->type == utils::PROBLEM_TYPE_2D){
        TopoDS_Shape copy = BRepBuilderAPI_Copy(this->data->ground_structure->shape);

        size_t graph_size = graph.size();
        for(size_t i = 0; i < graph_size; ++i){
            BeamNode* n = graph.get(i);
            double Fx = n->results[0];
            double Fy = n->results[1];
            double Mz = n->results[2];

            gp_Vec normal(n->normal);
            gp_Vec F(Fx, Fy, 0);
            double t = this->data->thickness;

            std::vector<double> S = this->data->material->get_max_stresses(normal);

            double S_f = std::min(S[0], S[1]);
            double S_n = (normal.Dot(F) < 0) ? S[0] : S[1];
            double S_c = S[2];

            // Bending
            double h_f = std::sqrt(6*Mz/(t*S_f));
            // Normal
            double h_n = std::abs(normal.Dot(F))/(t*S_n);
            // Shear
            double h_c = (F - normal.Dot(F)*normal).Magnitude()*(3/(2*t*S_c));

            // debug
            std::cout << h_f << " " << h_n << " " << h_c << " " << n->dim << std::endl;

            double h = 2*std::max({h_f, h_n, h_c, n->dim}); 
            gp_Ax2 axis(n->point, gp_Dir(0,0,1));
            gp_Circ circ(axis, h/2);
            TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(circ);
            TopoDS_Wire wire = BRepBuilderAPI_MakeWire(edge);
            TopoDS_Face face = BRepBuilderAPI_MakeFace(wire);

            // Fastest method I managed to find, cut from original shape then
            // inverse cut back into original shape.
            // Still a bit slow though.
            copy = BRepAlgoAPI_Cut(copy, face);
        }

        copy = BRepAlgoAPI_Cut(this->data->ground_structure->shape, copy);

        // Simplify geometry

        double tol = 10;
        double prec = tol*1.5;

        Handle(ShapeFix_Shape) sfs = new ShapeFix_Shape(copy);
        sfs->SetPrecision(prec);
        sfs->SetMaxTolerance(2*tol);
        sfs->SetMinTolerance(tol);
        auto sfw = sfs->FixWireTool();
        sfw->ModifyGeometryMode() = true;
        sfw->ModifyTopologyMode() = true;
        sfw->FixSmallMode() = true;
        sfw->FixSmall(false, prec);
        sfs->Perform();
        copy = sfs->Shape();

        Handle(ShapeFix_Wireframe) SFWF = new ShapeFix_Wireframe(copy);
        SFWF->SetPrecision(prec);
        SFWF->SetMaxTolerance(2*tol);
        SFWF->SetMinTolerance(tol);
        SFWF->ModeDropSmallEdges() = Standard_True;
        SFWF->FixSmallEdges();
        SFWF->FixWireGaps();
        copy = SFWF->Shape();

        return copy;
    }

    return TopoDS_Shape();
}

}
