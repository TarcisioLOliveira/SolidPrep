/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
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

#include <chrono>
#include <cblas.h>
#include <BRepBuilderAPI_Copy.hxx>
#include "project_data.hpp"
#include "optimization/MMASolver.hpp"
#include "topology_optimization/compliance_constraint_simple.hpp"

namespace topology_optimization {

ComplianceConstraintSimple::ComplianceConstraintSimple(DensityFilter* filter, Projection* projection, double c_max, ProjectData* data, double rho_init, double xtol_abs, double result_threshold, bool save, int P, int pc):
    c_max(c_max), data(data), rho_init(rho_init), xtol_abs(xtol_abs), result_threshold(result_threshold), save_result(save), P(P), pc(pc), filter(filter), projection(projection), viz(nullptr), fem(nullptr), mesh(nullptr), elem_number(0){}

void ComplianceConstraintSimple::initialize_views(Visualization* viz){
    this->viz = viz;

    this->stress_view = viz->add_view("Von Mises Stress", ViewHandler::ViewType::ELEMENTAL, ViewHandler::DataType::STRESS);
    this->density_view = viz->add_view("Density", ViewHandler::ViewType::ELEMENTAL, ViewHandler::DataType::DENSITY);
}

TopoDS_Shape ComplianceConstraintSimple::optimize(FiniteElement* fem, Meshing* mesh){

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

    this->fem->set_steps(2);

    auto V_it = this->grad_V.begin();
    for(const auto& g:this->mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                *V_it = e->get_volume(data->thickness);
                ++V_it;
            }
        }
    }

    std::vector<double> x(x_size, this->rho_init);
    std::vector<double> new_x(x_size, this->rho_init);

    this->filter->initialize(mesh, x_size);

    auto start_to = std::chrono::high_resolution_clock::now();
    
    //optimization::GCMMASolver gcmma(x.size(), 1, 0, 1e6, 1); //1e5
    optimization::MMASolver mma(x.size(), 1, 0, 1e6, 1); //1e5
    mma.SetAsymptotes(0.05, 0.7, 1.2);

    double ff;
    std::vector<double> dftmp(x.size());
    std::vector<double> dgtmp(x.size());
    std::vector<double> df(x.size());
    std::vector<double> dg(x.size());
    std::vector<double> g(1);

    // std::vector<double> xnew(x);
    std::vector<double> xold(x);

    std::vector<double> xmin;
    std::vector<double> xmax;

    xmin = std::vector<double>(x.size(), 0.0);
    xmax = std::vector<double>(x.size(), 1.0);

    logger::quick_log("Done.");

    double fnew = 0;
    std::vector<double> gnew(1);
    g[0] = 1e3;

    this->projection->project_densities(new_x);
    ff = this->fobj_grad(new_x, dftmp);
    this->projection->project_gradient(dftmp, new_x);
    this->filter->filter_gradient(dftmp, df);
    g[0] = this->fc_grad(new_x, dgtmp);
    this->projection->project_gradient(dgtmp, new_x);
    this->filter->filter_gradient(dgtmp, dg);

    // int max_innerit = 30;
    double ch = 1.0;
    int iter;
	for (iter = 0; ch > this->xtol_abs; ++iter){

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
        g[0] = this->fc_grad(new_x, dgtmp);
        this->projection->project_gradient(dgtmp, new_x);
        this->filter->filter_gradient(dgtmp, dg);

        logger::quick_log("");
        logger::quick_log("");
        logger::quick_log("Iteration: ", iter);
        logger::quick_log("Results: ", ff, g[0]);//+this->Smax);
        logger::quick_log("");
        logger::quick_log("Design var change: ", ch);
	}
    
    auto stop_to = std::chrono::high_resolution_clock::now();
    logger::quick_log("Final compliance: ", ff);
    auto to_duration = std::chrono::duration_cast<std::chrono::seconds>(stop_to-start_to);
    double to_time = to_duration.count();
    double it_time = to_time/iter;
    logger::quick_log("Time per iteration (topology optimization): ", it_time, " seconds");
    logger::quick_log("Number of iterations (topology optimization): ", iter);
   
    logger::quick_log(" "); 
    if(this->save_result){
        logger::quick_log("Saving resulting topology...");
        std::cout << "\r" << 0 << "%         ";
        TopoDS_Shape result = BRepBuilderAPI_Copy(this->data->geometries[0]->shape);
        for(size_t i = 0; i < new_x.size(); ++i){
            if(new_x[i] >= this->result_threshold){
                result = utils::cut_shape(result, this->data->geometries[0]->mesh[i]->get_shape());
            }
            double pc = i/(double)(new_x.size()-1);
            std::cout << "\r" << pc*100 << "%         ";
        }
        result = utils::cut_shape(this->data->geometries[0]->shape, result);
        return result;
        std::cout << "\r" << 100 << "%         ";
        logger::quick_log(" "); 
    }

    return TopoDS_Shape();
}

double ComplianceConstraintSimple::fobj_grad(const std::vector<double>& x, std::vector<double>& grad){

    double pc = this->pc;
    double pt = 1.0/2;

    this->u = this->fem->calculate_displacements(this->mesh, this->mesh->load_vector, x, pc);

    // Calculating stresses
    std::vector<double> fl(u.size(), 0);

    // Calculating global stress
    double Spn = 0;
    double Smax = 0;

    std::vector<double> stress_list(elem_number);
    auto x_it = x.cbegin();
    auto v_it = this->grad_V.cbegin();
    auto stress_it = stress_list.begin();
    for(const auto& g:this->mesh->geometries){
        const size_t num_den = g->number_of_densities_needed();
        const size_t num_mat = g->number_of_materials();
        if(g->do_topopt){
            if(num_den == 1){
                if(num_mat == 1){
                    const auto D = g->get_D(0);
                    for(const auto& e:g->mesh){
                        double S = e->get_stress_at(D, e->get_centroid(), u);
                        double Se = std::pow(*x_it, pt)*S;
                        *stress_it = *x_it*S;//std::pow(this->new_x[i],pc)*S;//Se;
                        if(Se > Smax){
                            Smax = Se;
                        }
                        double v = *v_it;
                        Spn += v*std::pow(Se, P);
                        e->get_virtual_load(D, v*std::pow(*x_it, pt*P)*std::pow(S, P-2), e->get_centroid(), u, fl);

                        ++x_it;
                        ++v_it;
                        ++stress_it;
                    }
                } else {
                    for(const auto& e:g->mesh){
                        const auto D = g->get_D_topopt(*x_it, 1);
                        double S = e->get_stress_at(D, e->get_centroid(), u);
                        // No relaxation when there is no void?
                        // Will do for now
                        double Se = S;
                        *stress_it = S;
                        if(Se > Smax){
                            Smax = Se;
                        }
                        double v = *v_it;
                        Spn += v*std::pow(Se, P);
                        e->get_virtual_load(D, v*std::pow(*x_it, pt*P)*std::pow(S, P-2), e->get_centroid(), u, fl);

                        ++x_it;
                        ++v_it;
                        ++stress_it;
                    }
                }
            } else {
                logger::log_assert(num_den == 1, logger::ERROR, "minimal volume problems that require more than 1 design variables are currently not supported.");
            }
        } else {
            const auto D = g->get_D(0);
            for(const auto& e:g->mesh){
                double S = e->get_stress_at(D, e->get_centroid(), u);
                *stress_it = S;
                ++stress_it;
            }
        }
    }
    this->stress_view->update_view(stress_list);
    this->density_view->update_view(x);
    this->viz->redraw();

    Spn = std::pow(Spn, 1.0/P);

    double result = Spn;

    double Sg = std::pow(Spn, 1 - P);

    logger::quick_log("Calculating adjoint problem...{");
    auto l = this->fem->calculate_displacements(this->mesh, fl, x, pc);
    logger::quick_log("} Done.");

    logger::quick_log("Calculating stress gradient...");

    x_it = x.cbegin();
    v_it = this->grad_V.cbegin();
    auto grad_it = grad.begin();
    for(const auto& g:this->mesh->geometries){
        const size_t num_den = g->number_of_densities_needed();
        const size_t num_mat = g->number_of_materials();
        if(g->do_topopt){
            if(num_den == 1){
                const auto D = g->get_D(0);
                auto x_it2    = x_it;
                auto v_it2    = v_it;
                auto grad_it2 = grad_it;
                for(const auto& e:g->mesh){
                    double lKu = pc*std::pow(*x_it, pc-1)*e->get_compliance(D, this->data->thickness, u, l);
                    double v = *v_it;
                    double S = e->get_stress_at(D, e->get_centroid(), u);
                    double Se = pt*v*std::pow(*x_it, pt*P-1)*std::pow(S, P);

                    *grad_it = Sg*(Se - lKu);

                    ++x_it;
                    ++v_it;
                    ++grad_it;
                }
                if(num_mat == 2){
                    const auto D = g->get_D(1);
                    for(const auto& e:g->mesh){
                        double lKu = -pc*std::pow(1 - *x_it2, pc-1)*e->get_compliance(D, this->data->thickness, u, l);
                        double v = *v_it2;
                        double S = e->get_stress_at(D, e->get_centroid(), u);
                        double Se = -pt*v*std::pow(1 - *x_it2, pt*P-1)*std::pow(S, P);

                        *grad_it2 += Sg*(Se - lKu);

                        ++x_it2;
                        ++v_it2;
                        ++grad_it2;
                    }
                }
            } else {
                logger::log_assert(num_den == 1, logger::ERROR, "minimal volume problems that require more than 1 design variables are currently not supported.");
            }
        }
    }
    logger::quick_log("Done.");

    return result;
}

double ComplianceConstraintSimple::fc_grad(const std::vector<double>& x, std::vector<double>& grad){
    double pc = this->pc;
    double c = 0;

    size_t i = 0;
    for(const auto& g:this->mesh->geometries){
        const size_t num_den = g->number_of_densities_needed();
        const size_t num_mat = g->number_of_materials();
        if(g->do_topopt){
            if(num_den == 1){
                const auto D = g->get_D(0);
                size_t j = i;
                for(const auto& e:g->mesh){
                    double uKu = e->get_compliance(D, this->mesh->thickness, this->u);
                    grad[i] = -pc*std::pow(x[i], pc-1)*uKu;

                    ++i;
                }
                if(num_mat == 2){
                    const auto D = g->get_D(1);
                    for(const auto& e:g->mesh){
                        double uKu = e->get_compliance(D, this->mesh->thickness, this->u);
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

    return c - this->c_max;
}

}
