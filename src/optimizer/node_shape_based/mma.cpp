/*
 *   Copyright (C) 2024 Tarcísio Ladeia de Oliveira.
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

#include <mpi.h>
#include "optimizer/node_shape_based/mma.hpp"
#include "logger.hpp"
#include "optimization/MMASolver.hpp"
#include "project_data.hpp"

namespace optimizer::node_shape_based{

MMA::MMA(const projspec::DataMap& data):
    NodeShapeBasedOptimizer(std::move(data.proj->shape_handler)),
    data(data.proj),
    xtol_abs(data.get_double("xtol_abs", 1e-10)),
    ftol_rel(data.get_double("ftol_rel", 1e-10)),
    asyminit(data.get_double("asyminit")),
    asymdec (data.get_double("asymdec")),
    asyminc (data.get_double("asyminc")),
    minfac  (data.get_double("minfac")),
    maxfac  (data.get_double("maxfac")),
    c(data.get_double("c")),
    save_result(data.get_bool("save_result", false)),
    viz(nullptr)
{

}

void MMA::initialize_views(Visualization* viz){
    this->viz = viz;

    this->node_view = viz->add_view("Shape Displacement", spview::defs::ViewType::VECTOR, spview::defs::DataType::DISPLACEMENT);
    this->stress_view = viz->add_view("Von Mises Stress", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);

    for(auto& f:this->objective){
        f->initialize_views(viz);
    }
    for(auto& f:this->constraints){
        f.fun->initialize_views(viz);
    }
}

TopoDS_Shape MMA::optimize(SolverManager* fem, Meshing* mesh){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id == 0){
        logger::quick_log("Preparing for optimization...");
    }

    if(mpi_id == 0){
        this->initialize_optimizer(mesh);
        this->shape_handler.obtain_affected_nodes();
    }
    auto stress_render = this->stresses;

    size_t x_size = this->shape_handler.get_number_of_variables();
    if(mpi_id == 0){
        for(auto& f:this->objective){
            f->initialize(this);
        }
        for(auto& f:this->constraints){
            f.fun->initialize(this);
        }
    }

    auto start_to = std::chrono::high_resolution_clock::now();

    std::vector<double> x(x_size, 0);

    size_t M = 0;
    for(auto& f:this->constraints){
        for(auto& t:f.types){
            if(t == Constraint::Type::GREATER_THAN || t == Constraint::Type::LESS_THAN){
                ++M;
            } else if(t == Constraint::Type::EQUAL){
                M += 2;
            }
        }
    }

    optimization::MMASolver mma(x_size, M, 0, c, 1);
    mma.SetAsymptotes(this->asyminit, this->asymdec, this->asyminc);
    mma.SetBoundFactors(this->minfac, this->maxfac);

    double ff = 1e-7;
    std::vector<double> dftmp(x.size());
    std::vector<double> dgtmp(x.size());
    std::vector<double> df(x.size());
    std::vector<double> dg(x.size()*M);
    std::vector<double> g(M);

    std::vector<double> xmin(x.size(), -1.0);
    std::vector<double> xmax(x.size(),  1.0);

    if(mpi_id == 0){
        logger::quick_log("Done.");
    }

    double fnew = 0;
    std::vector<double> u(mesh->max_dofs, 0);
    auto loads = mesh->load_vector;

    // int max_innerit = 30;
    double ch = 1.0;
    int iter;
	for (iter = 0; (ch > this->xtol_abs && (this->objective_weights.front() == 0 || std::abs((ff-fnew)/ff) > this->ftol_rel)); ++iter){

        fnew = ff;

        fem->generate_matrix(mesh);
        for(size_t i = 0; i < mesh->node_positions.size(); ++i){
            auto& lv = mesh->load_vector[i];
            auto& l = loads[i];
            std::copy(lv.begin(), lv.end(), l.begin());
        }
        fem->calculate_displacements_global(mesh, loads, u);
        if(mpi_id == 0){
            this->get_stresses(mesh->geometries, false, u, fem->D_vec, this->stresses);
            std::copy(stresses.begin(), stresses.end(), stress_render.begin());
            this->stress_view->update_view(stress_render);
        }
        if(this->objective.size() == 1){
            ff = this->objective_weights.front()*this->objective.front()->calculate_with_gradient(this, u, df);
            if(mpi_id == 0){
                if(this->objective_weights.front() != 1.0){
                    for(auto& v:df){
                        v *= this->objective_weights.front();
                    }
                }
            }
        } else {
            ff = 0;
            std::fill(df.begin(), df.end(), 0);
            for(size_t i = 0; i < this->objective.size(); ++i){
                ff += this->objective_weights[i]*this->objective[i]->calculate_with_gradient(this, u, dftmp);
                if(mpi_id == 0){
                    for(size_t j = 0; j < dftmp.size(); ++j){
                        df[j] += this->objective_weights[i]*dftmp[j];
                    }
                }
            }
        }
        this->shape_handler.filter_gradient(df);

        MPI_Bcast(&ff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if(M == 1){
            auto& c = this->constraints[0];
            if(c.types[0] == Constraint::Type::LESS_THAN){
                g[0] = c.fun->calculate_with_gradient(this, u, dg)
                       - c.bounds[0];
            } else if(c.types[0] == Constraint::Type::GREATER_THAN){
                g[0] = - c.fun->calculate_with_gradient(this, u, dg)
                       + c.bounds[0];
            }
            if(mpi_id == 0){
                if(c.types[0] == Constraint::Type::GREATER_THAN){
                    for(auto& d:dg){
                        d *= -1.0;
                    }
                }
            }
            this->shape_handler.filter_gradient(dg);
        } else if(M > 1){
            size_t g_id = 0;
            for(size_t i = 0; i < this->constraints.size(); ++i){
                auto& c = this->constraints[i];
                double val = c.fun->calculate_with_gradient(this, u, dgtmp);

                this->shape_handler.filter_gradient(dgtmp);

                for(size_t k = 0; k < c.types.size(); ++k){
                    if(c.types[k] == Constraint::Type::LESS_THAN){
                        g[g_id] = val - c.bounds[k];
                        if(mpi_id == 0){
                            for(size_t j = 0; j < x_size; ++j){
                                dg[M*j + g_id] = dgtmp[j];
                            }
                        }
                        ++g_id;
                    } else if(c.types[k] == Constraint::Type::GREATER_THAN){
                        g[g_id] = - val + c.bounds[k];
                        if(mpi_id == 0){
                            for(size_t j = 0; j < x_size; ++j){
                                dg[M*j + g_id] = -dgtmp[j];
                            }
                        }
                        ++g_id;
                    } else if(c.types[k] == Constraint::Type::EQUAL){
                        g[g_id] = val - c.bounds[k];
                        g[g_id + 1] = -g[g_id];
                        if(mpi_id == 0){
                            for(size_t j = 0; j < x_size; ++j){
                                dg[M*j + g_id] = dgtmp[j];
                                dg[M*j + g_id + 1] = -dgtmp[j];
                            }
                        }
                        g_id += 2;
                    }
                }
            }
        }

        if(mpi_id == 0){
            for(auto& f:this->objective){
                f->update();
            }
            for(auto& f:this->constraints){
                f.fun->update();
            }

            mma.Update(x.data(), df.data(), g.data(), dg.data(), xmin.data(), xmax.data());

            ch = 0.0;
            for(size_t i = 0; i < x.size(); ++i){
                ch = std::max(ch, std::abs(x[i]));
            }

            this->shape_handler.filter_displacement(x);
            //for (size_t i = 0; i < x.size(); ++i) {
            //    if(std::abs(x[i]) > 1e-4){
            //        logger::quick_log(i, std::abs(x[i]));
            //    }
            //}

            this->shape_handler.update_nodes(x);

            this->node_view->update_view(this->shape_handler.get_shape_displacement());

            // Update volumes after resizing elements
            this->get_volumes(mesh->geometries, mesh->thickness, this->volumes);

            std::fill(x.begin(), x.end(), 0);
        }
        MPI_Bcast(&ch, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(mpi_id == 0){
            logger::quick_log("");
            logger::quick_log("");
            logger::quick_log("Iteration: ", iter);
            logger::quick_log("Results: ", ff);
            logger::quick_log(g);
            logger::quick_log("");
            logger::quick_log("Design var change: ", ch);
            logger::quick_log("fobj change: ", std::abs((ff-fnew)/ff));
            logger::quick_log("");
        }
	}

    if(mpi_id == 0){
        logger::quick_log("");
        auto stop_to = std::chrono::high_resolution_clock::now();
        logger::quick_log("Final result: ", ff);
        auto to_duration = std::chrono::duration_cast<std::chrono::seconds>(stop_to-start_to);
        double to_time = to_duration.count();
        double it_time = to_time/iter;
        logger::quick_log("Time (topology optimization): ", to_time, " seconds");
        logger::quick_log("Time per iteration (topology optimization): ", it_time, " seconds");
        logger::quick_log("Number of iterations (topology optimization): ", iter);

        if(this->save_result){
            return this->make_shape(mesh->geometries, this->data->type);
        }
    }

    return TopoDS_Shape();
}

using namespace projspec;
const bool MMA::reg = Factory<NodeShapeBasedOptimizer>::add(
    [](const DataMap& data){
        return std::make_unique<MMA>(data);
    },
    ObjectRequirements{
        "mma",
        {
            DataEntry{.name = "xtol_abs", .type = TYPE_DOUBLE, .required = false},
            DataEntry{.name = "ftol_rel", .type = TYPE_DOUBLE, .required = false},
            DataEntry{.name = "save_result", .type = TYPE_BOOL, .required = false},
            DataEntry{.name = "asyminit", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "asymdec", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "asyminc", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "minfac", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "maxfac", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "c", .type = TYPE_DOUBLE, .required = true},
        }
    }
);

}
