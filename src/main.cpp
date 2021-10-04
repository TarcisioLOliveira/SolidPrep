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
#include <chrono>

int main(int argc, char* argv[]){
    // Bnd_Box bounds;
    // BRepBndLib::Add(this->shape, bounds);
    // bounds.SetGap(0.0);
    // Standard_Real fXMin, fYMin, fZMin, fXMax, fYMax, fZMax;
    // bounds.Get(fXMin, fYMin, fZMin, fXMax, fYMax, fZMax);
    // std::cout << fXMax << " " << fXMin << " " << fYMax << " " << fYMin << " " << fZMax << " " << fZMin << std::endl;

    double size_time = 0;
    ProjectData proj(argv[1]);
    auto start_sizing = std::chrono::high_resolution_clock::now();
    finite_element::DirectSolver fem_beam;
    sizing::StandardSizing* sizer = new sizing::StandardSizing(&proj, &fem_beam);
    auto shape = sizer->run();
    auto stop_sizing = std::chrono::high_resolution_clock::now();
    auto sizing_duration = std::chrono::duration_cast<std::chrono::seconds>(stop_sizing-start_sizing);
    size_time = sizing_duration.count()/60.0;
    delete sizer;
    utils::shape_to_file("test.step", shape);
    //return 0;

    meshing::Gmsh mesh(3, 1, utils::PROBLEM_TYPE_2D);
    // auto shape2 = proj.sizer->run();
    // utils::shape_to_file("test.step", shape2);
    // auto m = mesh.mesh(shape);
    //
    // auto m = mesh.mesh(proj.ground_structure->shape);

    // GroundStructure beam_test("beams.step", 1, utils::PROBLEM_TYPE_2D);
    // auto m = mesh.mesh(beam_test.shape);
    
    auto start_mesh = std::chrono::high_resolution_clock::now();
    //auto m = mesh.mesh(proj.ground_structure->shape);
    auto m = mesh.mesh(shape);
    std::vector<MeshElement*> elems;
    std::vector<double> loads;
    mesh.prepare_for_FEM(m, MeshElementFactory::GT9, &proj);
    finite_element::DirectSolver fem;
    //auto u = fem.calculate_displacements(&proj, &mesh);
    topology_optimization::MinimalVolume mv(6, 10, &proj);//proj.material->get_max_Von_Mises_2D(), &proj);
    auto stop_mesh = std::chrono::high_resolution_clock::now();
    // std::vector<double> stresses;
    // stresses.reserve(mesh.element_list.size());
    // for(auto& e:mesh.element_list){
    //     stresses.push_back(e->get_stress_at(e->get_centroid(), u));
    // }

    Visualization v;
    v.start();
    v.load_mesh(&mesh, proj.type);
    // v.update_stress_view(stresses);

    v.show();
    auto f = [&](){
      v.wait();
    };
    std::thread t(f);
    auto start_to = std::chrono::high_resolution_clock::now();
    TopoDS_Shape result = mv.optimize(&v, &fem, &mesh);
    auto stop_to = std::chrono::high_resolution_clock::now();
    auto to_duration = std::chrono::duration_cast<std::chrono::seconds>(stop_to-start_to);
    double to_time = to_duration.count()/60.0;
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(stop_to-start_sizing);
    double total_time = total_duration.count()/60.0;
    auto mesh_duration = std::chrono::duration_cast<std::chrono::seconds>(stop_mesh-start_mesh);
    double mesh_time = mesh_duration.count()/60.0;

    logger::quick_log("Sizing time: ", size_time, " minutes");
    logger::quick_log("Topology optimization time (including preparations for TO): ", to_time+mesh_time, " minutes");
    logger::quick_log("Total optimization time: ", to_time+size_time+mesh_time, " minutes");
    logger::quick_log("Total runtime (including GUI loading, excluding saving result as STEP): ", total_time, " minutes");

    //utils::shape_to_file("result.step", result);

    logger::quick_log("Finished.");

    t.join();

    v.end();


    // TopoDS_Shape s = proj.sizer->run();
    // utils::shape_to_file("test.step", s);

    return 0;
}
