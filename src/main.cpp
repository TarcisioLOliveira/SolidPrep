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

#include <iostream>
#include "STEPCAFControl_Reader.hxx"
#include "BRepClass3d_SolidClassifier.hxx"
#include "project_data.hpp"

int main(int argc, char* argv[]){
    // STEPControl_Reader reader;
    // IFSelect_ReturnStatus stat = reader.ReadFile("tmp/square.step");
    // reader.PrintCheckLoad(true, IFSelect_ItemsByEntity);

    // Standard_Integer NbRoots = reader.NbRootsForTransfer();
    // Standard_Integer num = reader.TransferRoots();
    // reader.PrintCheckTransfer(true, IFSelect_ItemsByEntity);
    // TopoDS_Shape result = reader.OneShape();

    // gp_Pnt p(0.0, 0.0, 0.0);
    // BRepClass3d_SolidClassifier insider(result);
    // insider.Perform(p, 0);
    // if(insider.State() == TopAbs_IN || insider.State() == TopAbs_ON){
    //     std::cout << "inside" << std::endl;
    // }
    
    ProjectData proj(argv[0]);

    return 0;
}
