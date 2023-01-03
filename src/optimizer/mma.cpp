/*
 *   Copyright (C) 2023 Tarcísio Ladeia de Oliveira.
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

#include "optimizer/mma.hpp"
#include "logger.hpp"
#include "optimization/MMASolver.hpp"
#include <chrono>

namespace optimizer{

MMA::MMA(DensityFilter* filter, Projection* projection, ProjectData* data, std::unique_ptr<DensityBasedFunction> objective, std::vector<std::unique_ptr<DensityBasedFunction>> constraints, std::vector<double> constraint_bounds, double pc, double rho_init, double xtol_abs, double ftol_rel, double result_threshold, bool save):
    data(data), rho_init(rho_init), xtol_abs(xtol_abs), ftol_rel(ftol_rel), pc(pc), result_threshold(result_threshold), save_result(save), objective(std::move(objective)), constraints(std::move(constraints)), constraint_bounds(std::move(constraint_bounds)), filter(filter), projection(projection), viz(nullptr)
    {}

void MMA::initialize_views(Visualization* viz){
    this->viz = viz;

    this->stress_view = viz->add_view("Von Mises Stress", ViewHandler::ViewType::ELEMENTAL, ViewHandler::DataType::STRESS);
    this->density_view = viz->add_view("Elemental Density", ViewHandler::ViewType::ELEMENTAL, ViewHandler::DataType::DENSITY);
}

TopoDS_Shape MMA::optimize(FiniteElement* fem, Meshing* mesh){
    logger::quick_log("Preparing for optimization...");

    size_t fem_steps = 1;
    fem_steps += this->objective->additional_steps();
    for(auto& f:this->constraints){
        fem_steps += f->additional_steps();
    }

    fem->set_steps(fem_steps);

    this->objective->initialize();
    for(auto& f:this->constraints){
        f->initialize();
    }

    size_t x_size = 0;
    size_t elem_number = 0;
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            x_size += g->mesh.size();
        }
        elem_number += g->mesh.size();
    }

    auto start_to = std::chrono::high_resolution_clock::now();

    std::vector<double> x(x_size, this->rho_init);
    std::vector<double> new_x(x_size, this->rho_init);

    this->filter->initialize(mesh, x_size);

    size_t M = this->constraints.size();

    optimization::MMASolver mma(x_size, M, 0, 1e6, 1); //1e5
    mma.SetAsymptotes(0.05, 0.7, 1.2);

    double ff;
    std::vector<double> dftmp(x.size());
    std::vector<double> dgtmp1(x.size());
    std::vector<double> dgtmp2(x.size());
    std::vector<double> df(x.size());
    std::vector<double> dg(x.size()*M);
    std::vector<double> g(M);

    // std::vector<double> xnew(x);
    std::vector<double> xold(x);

    std::vector<double> xmin;
    std::vector<double> xmax;

    xmin = std::vector<double>(x.size(), 0.0);
    xmax = std::vector<double>(x.size(), 1.0);

    logger::quick_log("Done.");

    double fnew = 0;

    this->projection->project_densities(new_x);
    auto u = fem->calculate_displacements(mesh, mesh->load_vector, new_x, this->pc);

    ff = this->objective->calculate_with_gradient(u, new_x, dftmp);
    this->projection->project_gradient(dftmp, new_x);
    this->filter->filter_gradient(dftmp, df);

    if(M == 1){
        g[0] = this->constraints[0]->calculate_with_gradient(u,new_x, dgtmp1)
               - this->constraint_bounds[0];
        this->projection->project_gradient(dgtmp1, new_x);
        this->filter->filter_gradient(dgtmp1, dg);
    } else if(M > 1){
        for(size_t i = 0; i < M; ++i){
            g[i] = this->constraints[i]->calculate_with_gradient(u,new_x, dgtmp1)
                   - this->constraint_bounds[i];
            this->projection->project_gradient(dgtmp1, new_x);
            this->filter->filter_gradient(dgtmp1, dgtmp2);
            for(size_t j = 0; j < x_size; ++j){
                dg[M*j + i] = dgtmp2[j];
            }
        }
    }

    // int max_innerit = 30;
    double ch = 1.0;
    int iter;
	for (iter = 0; (ch > this->xtol_abs && std::abs(ff-fnew)/ff > this->ftol_rel); ++iter){

        this->objective->update();
        for(auto& f:this->constraints){
            f->update();
        }

        fnew = ff;

        mma.Update(x.data(), df.data(), g.data(), dg.data(), xmin.data(), xmax.data());
        this->projection->update(iter);

        ch = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            ch = std::max(ch, std::abs(xold[i] - x[i]));
            xold[i] = x[i];
        }

        this->projection->project_densities(new_x);
        u = fem->calculate_displacements(mesh, mesh->load_vector, new_x, this->pc);

        ff = this->objective->calculate_with_gradient(u, new_x, dftmp);
        this->projection->project_gradient(dftmp, new_x);
        this->filter->filter_gradient(dftmp, df);

        if(M == 1){
            g[0] = this->constraints[0]->calculate_with_gradient(u,new_x, dgtmp1)
                   - this->constraint_bounds[0];
            this->projection->project_gradient(dgtmp1, new_x);
            this->filter->filter_gradient(dgtmp1, dg);
        } else if(M > 1){
            for(size_t i = 0; i < M; ++i){
                g[i] = this->constraints[i]->calculate_with_gradient(u,new_x, dgtmp1)
                       - this->constraint_bounds[i];
                this->projection->project_gradient(dgtmp1, new_x);
                this->filter->filter_gradient(dgtmp1, dgtmp2);
                for(size_t j = 0; j < x_size; ++j){
                    dg[M*j + i] = dgtmp2[j];
                }
            }
        }

        logger::quick_log("");
        logger::quick_log("");
        logger::quick_log("Iteration: ", iter);
        logger::quick_log("Results: ", ff);
        logger::quick_log(g);
        logger::quick_log("");
        logger::quick_log("Design var change: ", ch);
        logger::quick_log("fobj change: ", std::abs(ff-fnew)/ff);
	}

    auto stop_to = std::chrono::high_resolution_clock::now();
    logger::quick_log("Final result: ", ff);
    auto to_duration = std::chrono::duration_cast<std::chrono::seconds>(stop_to-start_to);
    double to_time = to_duration.count();
    double it_time = to_time/iter;
    logger::quick_log("Time per iteration (topology optimization): ", it_time, " seconds");
    logger::quick_log("Number of iterations (topology optimization): ", iter);

    if(this->save_result){
        return this->make_shape(new_x, mesh->geometries, this->result_threshold);
    }

    return TopoDS_Shape();
}

}
