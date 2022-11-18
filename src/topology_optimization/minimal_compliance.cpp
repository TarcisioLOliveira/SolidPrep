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

#include "topology_optimization/minimal_compliance.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include <cmath>
#include <cblas.h>
#include <BRepBuilderAPI_Copy.hxx>
#include <chrono>
#include "project_data.hpp"
#include <lapacke.h>
#include <vector>
#include "optimization/MMASolver.hpp"

namespace topology_optimization{


MinimalCompliance::MinimalCompliance(DensityFilter* filter, Projection* projection, ProjectData* data, double Vfinal, double xtol_abs, double ftol_rel, double result_threshold, bool save, int pc):
    data(data), Vfinal(Vfinal), xtol_abs(xtol_abs), ftol_rel(ftol_rel), result_threshold(result_threshold), save_result(save), pc(pc), elem_number(0), filter(filter), projection(projection), viz(nullptr), fem(nullptr), mesh(nullptr), max_V(0), cur_V(0), alpha(1){}

void MinimalCompliance::initialize_views(Visualization* viz){
    this->viz = viz;

    this->density_view = viz->add_view("Elemental Density", ViewHandler::ViewType::ELEMENTAL, ViewHandler::DataType::MATERIAL);
    // this->disp_view = viz->add_view("Displacement", ViewHandler::ViewType::VECTOR, ViewHandler::DataType::DISPLACEMENT);
}

TopoDS_Shape MinimalCompliance::optimize(FiniteElement* fem, Meshing* mesh){

    logger::quick_log("Preparing for optimization...");

    this->fem = fem;
    this->mesh = mesh;

    size_t x_size = 0;
    this->elem_number = 0;
    for(const auto& g:this->mesh->geometries){
        if(g->do_topopt){
            x_size += g->mesh.size();
        }
        this->elem_number += g->mesh.size();
    }

    this->grad_V = std::vector<double>(x_size, 0);
    this->alpha = 1;

    this->fem->set_steps(1);

    auto v_it = this->grad_V.begin();
    for(const auto& g:this->mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                *v_it = e->get_volume(data->thickness);
                this->max_V += *v_it;

                ++v_it;
            }
        }
    }

    this->cur_V = this->max_V;

    std::vector<double> x(x_size, this->Vfinal);
    std::vector<double> new_x(x_size, this->Vfinal);

    this->filter->initialize(this->mesh, x_size);

    auto start_to = std::chrono::high_resolution_clock::now();

    optimization::MMASolver mma(x.size(), 1, 0, 1e4, 1); //1e5
    mma.SetAsymptotes(0.5, 0.7, 1.2);

    double ff;
    std::vector<double> df(x.size());
    std::vector<double> dg(x.size());
    std::vector<double> dftmp(x.size());
    std::vector<double> dgtmp(x.size());
    std::vector<double> g(1);

    std::vector<double> xold(x);

    std::vector<double> xmin;
    std::vector<double> xmax;

    xmin = std::vector<double>(x.size(), 0.0);
    xmax = std::vector<double>(x.size(), 1.0);

    logger::quick_log("Done.");

    double fnew = this->max_V;
    std::vector<double> gnew(1);
    g[0] = 1e3;

    this->projection->project_densities(new_x);
    ff = this->fobj_grad(new_x, dftmp);
    this->projection->project_gradient(dftmp, new_x);
    this->filter->filter_gradient(dftmp, df);
    g[0] = this->fc_norm_grad(new_x, dgtmp);
    this->projection->project_gradient(dgtmp, new_x);
    this->filter->filter_gradient(dgtmp, dg);

    // int max_innerit = 30;
    double ch = 1.0;
    int iter;
	for (iter = 0; (ch > this->xtol_abs && std::abs(ff-fnew)/fnew > this->ftol_rel); ++iter){

        fnew = ff;
        gnew = g;

        mma.Update(x.data(), df.data(), g.data(), dg.data(), xmin.data(), xmax.data());
        this->projection->update(iter);

        ch = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            ch = std::max(ch, std::abs(xold[i] - x[i]));
            xold[i] = x[i];
        }

        this->filter->filter_densities(x, new_x);
        this->projection->project_densities(new_x);
        ff = this->fobj_grad(new_x, dftmp);
        this->projection->project_gradient(dftmp, new_x);
        this->filter->filter_gradient(dftmp, df);
        g[0] = this->fc_norm_grad(new_x, dgtmp);
        this->projection->project_gradient(dgtmp, new_x);
        this->filter->filter_gradient(dgtmp, dg);

        logger::quick_log("");
        logger::quick_log("");
        logger::quick_log("Iteration: ", iter);
        logger::quick_log("Results: ", ff, g[0]);//+this->Smax);
        logger::quick_log("");
        logger::quick_log("Design var change: ", ch);
        logger::quick_log("Compliance change: ", std::abs(ff-fnew)/fnew);
        logger::quick_log("Volume change: ", std::abs(g[0]-gnew[0])/this->max_V);
	}
    
    auto stop_to = std::chrono::high_resolution_clock::now();
    logger::quick_log("Final volume: ", this->cur_V);
    auto to_duration = std::chrono::duration_cast<std::chrono::seconds>(stop_to-start_to);
    double to_time = to_duration.count();
    double it_time = to_time/iter;
    logger::quick_log("Time per iteration (topology optimization): ", it_time, " seconds");
    logger::quick_log("Number of iterations (topology optimization): ", iter);
   
    logger::quick_log(" "); 

    TopoDS_Shape result;
    if(this->save_result){
        logger::quick_log("Saving resulting topology...");
        std::cout << "\r" << 0 << "%         ";
        result = BRepBuilderAPI_Copy(this->data->geometries[0]->shape);
        for(size_t i = 0; i < x.size(); ++i){
            if(x[i] >= this->result_threshold){
                result = utils::cut_shape(result, this->data->geometries[0]->mesh[i]->get_shape());
            }
            double pc = i/(double)(x.size()-1);
            std::cout << "\r" << pc*100 << "%         ";
        }
        result = utils::cut_shape(this->data->geometries[0]->shape, result);
        std::cout << "\r" << 100 << "%         ";
        logger::quick_log(" "); 
    }

    // Validate results
    bool validate = true;
    for(const auto& g:this->mesh->geometries){
        if(g->number_of_materials() > 1){
            validate = false;
        }
    }
    if(validate){
        this->mesh->prune(this->data->forces, this->data->supports, new_x, this->result_threshold);
        std::vector<double> u = this->fem->calculate_displacements(this->mesh, this->mesh->load_vector);
        double c = cblas_ddot(u.size(), this->mesh->load_vector.data(), 1, u.data(), 1);
        logger::quick_log(" ");
        logger::quick_log("Compliance error: ", 100*(ff/c-1), "%");
        logger::quick_log(" ");
    }

    return result;
}

double MinimalCompliance::fobj(const std::vector<double>& x){
    double pc = this->pc;
    double c = 0;

    std::vector<double> u = this->fem->calculate_displacements(this->mesh, this->mesh->load_vector, x, pc);

    c = cblas_ddot(u.size(), this->mesh->load_vector.data(), 1, u.data(), 1);

    return c;

}
double MinimalCompliance::fobj_grad(const std::vector<double>& x, std::vector<double>& grad){
    double pc = this->pc;
    double c = 0;

    std::vector<double> u = this->fem->calculate_displacements(this->mesh, this->mesh->load_vector, x, pc);

    this->density_view->update_view(x);
    //this->disp_view->update_view(u);
    this->viz->redraw();

    size_t i = 0;
    for(const auto& g:this->mesh->geometries){
        const size_t num_den = g->number_of_densities_needed();
        const size_t num_mat = g->number_of_materials();
        if(g->do_topopt){
            if(num_den == 1){
                const auto D = g->get_D(0);
                size_t j = i;
                for(const auto& e:g->mesh){
                    double uKu = e->get_compliance(D, this->mesh->thickness, u);
                    grad[i] = -pc*std::pow(x[i], pc-1)*uKu;

                    ++i;
                }
                if(num_mat == 2){
                    const auto D = g->get_D(1);
                    for(const auto& e:g->mesh){
                        double uKu = e->get_compliance(D, this->mesh->thickness, u);
                        grad[j] += pc*std::pow(1 - x[j], pc-1)*uKu;

                        ++j;
                    }
                }
            } else {
                logger::log_assert(num_den == 1, logger::ERROR, "minimal compliance problems that require more than 1 design variables are currently not supported.");
            }
        }
    }
    c = cblas_ddot(u.size(), this->mesh->load_vector.data(), 1, u.data(), 1);

    return c;
}

double MinimalCompliance::fc_norm(const std::vector<double>& x){
    double V = 0;

    for(size_t i = 0; i < x.size(); ++i){
        V += x[i]*this->grad_V[i];
    }
    this->cur_V = V;

    return V - this->Vfinal*this->max_V;

}
double MinimalCompliance::fc_norm_grad(const std::vector<double>& x, std::vector<double>& grad){
    double V = 0;

    for(size_t i = 0; i < x.size(); ++i){
        V += x[i]*this->grad_V[i];
        grad[i] = this->grad_V[i];
    }
    this->cur_V = V;

    return V - this->Vfinal*this->max_V;
}

}
