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
#include <STEPCAFControl_Writer.hxx>
#include <STEPControl_StepModelType.hxx>
#include "project_data.hpp"
#include "sizing/beam_sizing.hpp"
#include "utils.hpp"

int main(int argc, char* argv[]){
    // Bnd_Box bounds;
    // BRepBndLib::Add(this->shape, bounds);
    // bounds.SetGap(0.0);
    // Standard_Real fXMin, fYMin, fZMin, fXMax, fYMax, fZMax;
    // bounds.Get(fXMin, fYMin, fZMin, fXMax, fYMax, fZMax);
    // std::cout << fXMax << " " << fXMin << " " << fYMax << " " << fYMin << " " << fZMax << " " << fZMin << std::endl;

    
    ProjectData proj(argv[1]);
    // auto path = proj.pathfinder->find_path(proj.forces[0], proj.supports[0].get_shape());
    // for(auto& p : path){
    //     std::cout << "(" << p.X() << ", " << p.Y() << ")" << std::endl;
    // }
    // proj.sizer.reset(new sizing::BeamSizing(&proj));

    TopoDS_Shape s = proj.sizer->run();
    utils::shape_to_file("test.step", s);

    return 0;
}
