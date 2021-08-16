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
#include "logger.hpp"
#include "project_data.hpp"
#include "sizing/beam_sizing.hpp"
#include "utils.hpp"
#include <gmsh.h>
#include <set>
#include "meshing/gmsh.hpp"
#include "finite_element/direct_solver.hpp"
#include "visualization.hpp"
#include "topology_optimization/minimal_volume.hpp"
#include <thread>
#include "sizing/standard_sizing.hpp"

int main(int argc, char* argv[]){
    // Bnd_Box bounds;
    // BRepBndLib::Add(this->shape, bounds);
    // bounds.SetGap(0.0);
    // Standard_Real fXMin, fYMin, fZMin, fXMax, fYMax, fZMax;
    // bounds.Get(fXMin, fYMin, fZMin, fXMax, fYMax, fZMax);
    // std::cout << fXMax << " " << fXMin << " " << fYMax << " " << fYMin << " " << fZMax << " " << fZMin << std::endl;

    ProjectData proj(argv[1]);
    finite_element::DirectSolver fem_beam;
    sizing::StandardSizing sizer(&proj, &fem_beam);
    auto shape = sizer.run();
    utils::shape_to_file("test.step", shape);
    return 0;

    meshing::Gmsh mesh(4, 1, utils::PROBLEM_TYPE_2D);
    // auto shape2 = proj.sizer->run();
    // utils::shape_to_file("test.step", shape2);
    // auto m = mesh.mesh(shape);
    //
    // auto m = mesh.mesh(proj.ground_structure->shape);
    auto m = mesh.mesh(shape);
    std::vector<MeshElement*> elems;
    std::vector<double> loads;
    mesh.prepare_for_FEM(m, MeshElementFactory::GT9, &proj);
    finite_element::DirectSolver fem;
    fem.calculate_displacements(&proj, &mesh);
    //topology_optimization::MinimalVolume mv(9, proj.material->get_max_Von_Mises_2D(), &proj);

    Visualization v;
    v.start();
    v.load_mesh(&mesh, proj.type);

    v.show();
    //auto f = [&](){
        v.wait();
    // };
    // std::thread t(f);
    // mv.optimize(&v, &fem, &mesh);

    logger::quick_log("Finished.");

    //t.join();

    v.end();


    // TopoDS_Shape s = proj.sizer->run();
    // utils::shape_to_file("test.step", s);

    return 0;
}
