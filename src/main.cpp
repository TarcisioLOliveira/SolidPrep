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

#include <bits/chrono.h>
#include <iostream>
#include <STEPCAFControl_Writer.hxx>
#include <STEPControl_StepModelType.hxx>
#include "logger.hpp"
#include "project_data.hpp"
#include "sizing/beam_sizing.hpp"
#include "utils.hpp"
#include <gmsh.h>
#include <petsc.h>
#include <ratio>
#include <set>
#include "meshing/gmsh.hpp"
#include "finite_element/direct_solver.hpp"
#include "visualization.hpp"
#include "topology_optimization/minimal_volume.hpp"
#include <string>
#include <thread>
#include "sizing/standard_sizing.hpp"
#include <chrono>
#include <mpich-x86_64/mpi.h>
#include <cblas.h>
#include <time.h>
#include <Eigen/Core>
#include <Message.hxx>
#include <Message_PrinterOStream.hxx>
#include "spview.hpp"

int main(int argc, char* argv[]){
    MPI_Init(NULL, NULL);

    // std::cout << "pid: " << getpid() << std::endl;
    // sleep(10);
    
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    PetscInitialize(&argc, &argv, NULL, NULL);

    if(mpi_id == 0){
        Eigen::initParallel();

        logger::log_assert(argc > 1, logger::ERROR, "missing path to configuration file.");

    }
    std::unique_ptr<Visualization> v;

    double size_time = 0;

    Message::DefaultMessenger()->RemovePrinters(STANDARD_TYPE(Message_PrinterOStream));

    // Load project file
    std::unique_ptr<ProjectData> proj(std::make_unique<ProjectData>(argv[1]));
    TopoDS_Shape shape;

    auto start_sizing = std::chrono::high_resolution_clock::now();
    if(mpi_id == 0){
        v.reset(new Visualization());
        v->start();
        if(proj->analysis == ProjectData::COMPLETE || proj->analysis == ProjectData::BEAMS_ONLY){
            shape = proj->sizer->run();
            auto stop_sizing = std::chrono::high_resolution_clock::now();
            auto sizing_duration = std::chrono::duration_cast<std::chrono::seconds>(stop_sizing-start_sizing);
            size_time = sizing_duration.count()/60.0;
            utils::shape_to_file("sized.step", shape);
            proj->geometries[0]->shape = shape;
        }
    }

    // Meshing
    auto start_mesh = std::chrono::high_resolution_clock::now();
    // May be removed to use Gmsh with MPI
    if(mpi_id == 0){
        proj->topopt_mesher->mesh(proj->forces, proj->supports, proj->springs);
        for(auto& f:proj->fields){
            f->generate();
        }
        proj->topopt_mesher->apply_boundary_conditions(proj->forces, proj->supports, proj->springs, proj->internal_loads, proj->sub_problems);
    }
    std::vector<MeshElement*> elems;
    std::vector<double> loads;
    auto stop_mesh = std::chrono::high_resolution_clock::now();
    if(proj->analysis == ProjectData::FEA_ONLY || proj->analysis == ProjectData::BEAMS_ONLY){

        // Finite element analysis
        auto start_fea = std::chrono::high_resolution_clock::now();
        auto l = proj->topopt_mesher->load_vector;
        std::vector<double> u(proj->topopt_mesher->max_dofs, 0);
        proj->topopt_fea->generate_matrix(proj->topopt_mesher.get());
        proj->topopt_fea->calculate_displacements_global(proj->topopt_mesher.get(), l, u);
        auto stop_fea = std::chrono::high_resolution_clock::now();

        if(mpi_id == 0){
            double fea_duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop_fea-start_fea).count()/60000.0;

            if(proj->type == utils::PROBLEM_TYPE_2D){
                std::vector<double> stresses;
                std::vector<double> stressesX;
                std::vector<double> stressesY;
                std::vector<double> stressesXY;
                std::vector<double> strainX;
                std::vector<double> strainY;
                std::vector<double> strainXY;
                size_t s_size = 0;
                for(auto& g:proj->geometries){
                    s_size += g->mesh.size();
                }
                stresses  .reserve(s_size);
                stressesX .reserve(s_size);
                stressesY .reserve(s_size);
                stressesXY.reserve(s_size);
                strainX .reserve(s_size);
                strainY .reserve(s_size);
                strainXY.reserve(s_size);

                for(auto& g:proj->geometries){
                    for(auto& e:g->mesh){
                        const gp_Pnt c = e->get_centroid();
                        const auto D = g->materials.get_D(e.get(), c);
                        stresses.push_back(e->get_stress_at(D, e->get_centroid(), u));
                        auto tensor = e->get_stress_tensor(D, e->get_centroid(), u);
                        stressesX.push_back(tensor[0]);
                        stressesY.push_back(tensor[3]);
                        stressesXY.push_back(tensor[1]);
                        tensor = e->get_strain_tensor(e->get_centroid(), u);
                        strainX.push_back(tensor[0]);
                        strainY.push_back(tensor[3]);
                        strainXY.push_back(tensor[1]);
                    }
                }


                if(proj->analysis == ProjectData::BEAMS_ONLY){
                    logger::quick_log("");
                    logger::quick_log("Sizing time: ", size_time, " minutes");
                    logger::quick_log("");
                }
                logger::quick_log("");
                logger::quick_log("FEA time: ", fea_duration, " minutes");
                logger::quick_log("");
                logger::quick_log("Compliance: ", cblas_ddot(u.size(), u.data(), 1, proj->topopt_mesher->global_load_vector.data(), 1));
                logger::quick_log("");

                // Display results
                v->load_mesh(proj->topopt_mesher.get(), proj->type);

                auto stressview_VM = v->add_view("Von Mises Stress",       spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto stressview_X  = v->add_view("Normal Stress (X axis)", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto stressview_Y  = v->add_view("Normal Stress (Y axis)", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto stressview_XY = v->add_view("Shear Stress",           spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto strainview_X  = v->add_view("Normal Strain (X axis)", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto strainview_Y  = v->add_view("Normal Strain (Y axis)", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto strainview_XY = v->add_view("Shear Strain",           spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto displ = v->add_view("Displacement",           spview::defs::ViewType::VECTOR, spview::defs::DataType::OTHER);
                auto force = v->add_view("Force",           spview::defs::ViewType::VECTOR, spview::defs::DataType::OTHER);

                for(auto& f:proj->fields){
                    f->initialize_views(v.get());
                }

                stressview_VM->update_view(stresses);
                stressview_X ->update_view(stressesX);
                stressview_Y ->update_view(stressesY);
                stressview_XY->update_view(stressesXY);
                strainview_X ->update_view(strainX);
                strainview_Y ->update_view(strainY);
                strainview_XY->update_view(strainXY);
                displ->update_view(u);
                force->update_view(proj->topopt_mesher->global_load_vector);

                for(auto& f:proj->fields){
                    f->display_views();
                }

                v->wait();
                v->end();
            } else if(proj->type == utils::PROBLEM_TYPE_3D){
                std::vector<double> stresses;
                std::vector<double> stressesX;
                std::vector<double> stressesY;
                std::vector<double> stressesZ;
                std::vector<double> stressesXY;
                std::vector<double> stressesXZ;
                std::vector<double> stressesYZ;
                std::vector<double> strainX;
                std::vector<double> strainY;
                std::vector<double> strainZ;
                std::vector<double> strainXY;
                std::vector<double> strainXZ;
                std::vector<double> strainYZ;
                size_t s_size = 0;
                for(auto& g:proj->geometries){
                    s_size += g->mesh.size();
                }
                stresses  .reserve(s_size);
                stressesX .reserve(s_size);
                stressesY .reserve(s_size);
                stressesZ .reserve(s_size);
                stressesXY.reserve(s_size);
                stressesXZ.reserve(s_size);
                stressesYZ.reserve(s_size);
                strainX .reserve(s_size);
                strainY .reserve(s_size);
                strainZ .reserve(s_size);
                strainXY.reserve(s_size);
                strainXZ.reserve(s_size);
                strainYZ.reserve(s_size);
                double max_stress = 0;
                gp_Pnt max_point(0,0,0);
                for(auto& g:proj->geometries){
                    for(auto& e:g->mesh){
                        const gp_Pnt c = e->get_centroid();
                        const auto D = g->materials.get_D(e.get(), c);
                        stresses.push_back(e->get_stress_at(D, e->get_centroid(), u));
                        if(stresses.back() > max_stress){
                            max_stress = stresses.back();
                            max_point = c;
                        }
                        auto tensor = e->get_stress_tensor(D, e->get_centroid(), u);
                        stressesX.push_back(tensor[0]);
                        stressesY.push_back(tensor[4]);
                        stressesZ.push_back(tensor[8]);
                        stressesXY.push_back(tensor[1]);
                        stressesXZ.push_back(tensor[2]);
                        stressesYZ.push_back(tensor[5]);
                        tensor = e->get_strain_tensor(e->get_centroid(), u);
                        strainX.push_back(tensor[0]);
                        strainY.push_back(tensor[4]);
                        strainZ.push_back(tensor[8]);
                        strainXY.push_back(tensor[1]);
                        strainXZ.push_back(tensor[2]);
                        strainYZ.push_back(tensor[5]);
                    }
                }

                logger::quick_log("Max stress: ", max_stress, " at ", max_point.X(), max_point.Y(), max_point.Z());
                if(proj->analysis == ProjectData::BEAMS_ONLY){
                    logger::quick_log("");
                    logger::quick_log("Sizing time: ", size_time, " minutes");
                    logger::quick_log("");
                }
                logger::quick_log("");
                logger::quick_log("FEA time: ", fea_duration, " minutes");
                logger::quick_log("");
                logger::quick_log("Compliance: ", cblas_ddot(u.size(), u.data(), 1, proj->topopt_mesher->global_load_vector.data(), 1));
                logger::quick_log("");

                // Display results
                v->load_mesh(proj->topopt_mesher.get(), proj->type);

                auto stressview_VM = v->add_view("Von Mises Stress",        spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto stressview_X  = v->add_view("Normal Stress (X axis)",  spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto stressview_Y  = v->add_view("Normal Stress (Y axis)",  spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto stressview_Z  = v->add_view("Normal Stress (Z axis)",  spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto stressview_XY = v->add_view("Shear Stress (XY plane)", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto stressview_XZ = v->add_view("Shear Stress (XZ plane)", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto stressview_YZ = v->add_view("Shear Stress (YZ plane)", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto strainview_X  = v->add_view("Normal Strain (X axis)",  spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto strainview_Y  = v->add_view("Normal Strain (Y axis)",  spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto strainview_Z  = v->add_view("Normal Strain (Z axis)",  spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto strainview_XY = v->add_view("Shear Strain (XY plane)", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto strainview_XZ = v->add_view("Shear Strain (XZ plane)", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto strainview_YZ = v->add_view("Shear Strain (YZ plane)", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
                auto displ = v->add_view("Displacement",           spview::defs::ViewType::VECTOR, spview::defs::DataType::OTHER);
                auto force = v->add_view("Force",           spview::defs::ViewType::VECTOR, spview::defs::DataType::OTHER);

                for(auto& f:proj->fields){
                    f->initialize_views(v.get());
                }

                stressview_VM->update_view(stresses);
                stressview_X ->update_view(stressesX);
                stressview_Y ->update_view(stressesY);
                stressview_Z ->update_view(stressesZ);
                stressview_XY->update_view(stressesXY);
                stressview_XZ->update_view(stressesXZ);
                stressview_YZ->update_view(stressesYZ);
                strainview_X ->update_view(strainX);
                strainview_Y ->update_view(strainY);
                strainview_Z ->update_view(strainZ);
                strainview_XY->update_view(strainXY);
                strainview_XZ->update_view(strainXZ);
                strainview_YZ->update_view(strainYZ);
                displ->update_view(u);
                force->update_view(proj->topopt_mesher->global_load_vector);

                for(auto& f:proj->fields){
                    f->display_views();
                }

                v->wait();
                v->end();
            }
        }
    } else if(proj->analysis == ProjectData::OPTIMIZE_ONLY || proj->analysis == ProjectData::COMPLETE){

        if(mpi_id == 0){
            // Display progress
            v->load_mesh(proj->topopt_mesher.get(), proj->type);

            //proj->topopt->initialize_views(&v);
            proj->optimizer->initialize_views(v.get());
            for(auto& f:proj->fields){
                f->initialize_views(v.get());
            }
            for(auto& f:proj->fields){
                f->display_views();
            }
        }

        auto start_to = std::chrono::high_resolution_clock::now();

        // Optimization
        //TopoDS_Shape result = proj->topopt->optimize(proj->topopt_fea.get(), proj->topopt_mesher.get());
        TopoDS_Shape result = proj->optimizer->optimize(proj->topopt_fea.get(), proj->topopt_mesher.get());

        if(mpi_id == 0){
            // Display time
            auto stop_to = std::chrono::high_resolution_clock::now();
            auto to_duration = std::chrono::duration_cast<std::chrono::seconds>(stop_to-start_to);
            double to_time = to_duration.count()/60.0;
            auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(stop_to-start_sizing);
            double total_time = total_duration.count()/60.0;
            auto mesh_duration = std::chrono::duration_cast<std::chrono::seconds>(stop_mesh-start_mesh);
            double mesh_time = mesh_duration.count()/60.0;

            if(result != TopoDS_Shape()){
                utils::shape_to_file("result.step", result);
            }

            if(proj->analysis == ProjectData::COMPLETE){
                logger::quick_log("Sizing time: ", size_time, " minutes");
            }
            logger::quick_log("Topology optimization time (including preparations for TO): ", to_time+mesh_time, " minutes");
            logger::quick_log("Total optimization time: ", to_time+size_time+mesh_time, " minutes");
            logger::quick_log("Total runtime (including GUI loading, excluding saving result as STEP): ", total_time, " minutes");

            logger::quick_log("Finished.");

            //t.join();
            v->wait();

            v->end();
        }
    }

    // Finalize MPI objects before MPI
    proj.reset(nullptr);

    PetscFinalize();

    MPI_Finalize();

    return 0;
}
