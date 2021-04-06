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

Support::Support(bool X, bool Y, bool MZ, CrossSection cross_section):
    X(X), Y(Y), Z(false), MX(false), MY(false), MZ(MZ), S(std::move(cross_section)), fdof(0), mdof(0){
    if(X) ++this->fdof;
    if(Y) ++this->fdof;
    if(MZ) ++this->mdof;
    
}

Support::Support(bool X, bool Y, bool Z, bool MX, bool MY, bool MZ, CrossSection cross_section):
    X(X), Y(Y), Z(Z), MX(MX), MY(MY), MZ(MZ), S(std::move(cross_section)), fdof(0), mdof(0){
    if(X) ++this->fdof;
    if(Y) ++this->fdof;
    if(Z) ++this->fdof;
    if(MX) ++this->mdof;
    if(MY) ++this->mdof;
    if(MZ) ++this->mdof;
}

