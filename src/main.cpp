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

#include <iostream>
#include <STEPCAFControl_Writer.hxx>
#include <STEPControl_StepModelType.hxx>
#include "logger.hpp"
#include "project_data.hpp"
#include "sizing/beam_sizing.hpp"
#include "utils.hpp"
#include <gmsh.h>
#include <ratio>
#include <set>
#include "meshing/gmsh.hpp"
#include "finite_element/direct_solver.hpp"
#include "visualization.hpp"
#include "topology_optimization/minimal_volume.hpp"
#include <thread>
#include "sizing/standard_sizing.hpp"
#include <chrono>

int main(int argc, char* argv[]){
    (void)argc;

    double size_time = 0;

    // Load project file
    ProjectData proj(argv[1]);
    auto start_sizing = std::chrono::high_resolution_clock::now();
    TopoDS_Shape shape;

    // Sizing step
    std::vector<ElementShape> m;

    if(proj.analysis == ProjectData::COMPLETE || proj.analysis == ProjectData::BEAMS_ONLY){
        shape = proj.sizer->run();
        auto stop_sizing = std::chrono::high_resolution_clock::now();
        auto sizing_duration = std::chrono::duration_cast<std::chrono::seconds>(stop_sizing-start_sizing);
        size_time = sizing_duration.count()/60.0;
        utils::shape_to_file("sized.step", shape);
        proj.geometries[0]->shape = shape;
    }

    // Meshing
    auto start_mesh = std::chrono::high_resolution_clock::now();
    proj.topopt_mesher->mesh(proj.forces, proj.supports);
    std::vector<MeshElement*> elems;
    std::vector<double> loads;
    auto stop_mesh = std::chrono::high_resolution_clock::now();
    if(proj.analysis == ProjectData::FEA_ONLY || proj.analysis == ProjectData::BEAMS_ONLY){

        // Finite element analysis
        auto start_fea = std::chrono::high_resolution_clock::now();
        auto u = proj.topopt_fea->calculate_displacements(proj.topopt_mesher.get(), proj.topopt_mesher->load_vector);
        auto stop_fea = std::chrono::high_resolution_clock::now();
        double fea_duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop_fea-start_fea).count()/60000.0;

        std::vector<double> stresses;
        std::vector<double> stressesX;
        std::vector<double> stressesY;
        std::vector<double> stressesXY;
        size_t s_size = 0;
        for(auto& g:proj.geometries){
            s_size += g->mesh.size();
        }
        stresses  .reserve(s_size);
        stressesX .reserve(s_size);
        stressesY .reserve(s_size);
        stressesXY.reserve(s_size);
        for(auto& g:proj.geometries){
            const auto D = g->get_D(0);
            for(auto& e:g->mesh){
                stresses.push_back(e->get_stress_at(D, e->get_centroid(), u));
                auto tensor = e->get_stress_tensor(D, e->get_centroid(), u);
                stressesX.push_back(tensor[0]);
                stressesY.push_back(tensor[3]);
                stressesXY.push_back(tensor[1]);
            }
        }

        if(proj.analysis == ProjectData::BEAMS_ONLY){
            logger::quick_log("");
            logger::quick_log("Sizing time: ", size_time, " minutes");
            logger::quick_log("");
        }
        logger::quick_log("");
        logger::quick_log("FEA time: ", fea_duration, " minutes");
        logger::quick_log("");

        // Display results
        Visualization v;
        v.start();
        v.load_mesh(proj.topopt_mesher.get(), proj.type);

        auto stressview_VM = v.add_view("Von Mises Stress",       ViewHandler::ViewType::ELEMENTAL, ViewHandler::DataType::STRESS);
        auto stressview_X  = v.add_view("Normal Stress (X axis)", ViewHandler::ViewType::ELEMENTAL, ViewHandler::DataType::STRESS);
        auto stressview_Y  = v.add_view("Normal Stress (Y axis)", ViewHandler::ViewType::ELEMENTAL, ViewHandler::DataType::STRESS);
        auto stressview_XY = v.add_view("Shear Stress",           ViewHandler::ViewType::ELEMENTAL, ViewHandler::DataType::STRESS);

        stressview_VM->update_view(stresses);
        stressview_X ->update_view(stressesX);
        stressview_Y ->update_view(stressesY);
        stressview_XY->update_view(stressesXY);
        v.show();
        v.wait();
        v.end();
    } else if(proj.analysis == ProjectData::OPTIMIZE_ONLY || proj.analysis == ProjectData::COMPLETE){

        // Display progress
        Visualization v;
        v.start();
        v.load_mesh(proj.topopt_mesher.get(), proj.type);

        v.show();
        auto f = [&](){
          v.wait();
        };
        std::thread t(f);

        auto start_to = std::chrono::high_resolution_clock::now();

        // Optimization
        TopoDS_Shape result = proj.topopt->optimize(&v, proj.topopt_fea.get(), proj.topopt_mesher.get());

        // Display time
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

        logger::quick_log("Finished.");

        //t.join();
        v.hide();
        v.show();
        v.wait();

        v.end();

    }

    return 0;
}
