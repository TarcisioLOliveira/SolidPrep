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

#include "topology_optimization/minimal_volume.hpp"
#include "element/TRI3.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include <cmath>
#include <cblas.h>
#include <BRepBuilderAPI_Copy.hxx>
#include <chrono>
#include "project_data.hpp"
#include <lapacke.h>
#include <vector>
#include "optimization/GCMMASolver.hpp"
#include "optimization/MMASolver.hpp"

namespace topology_optimization{


MinimalVolume::MinimalVolume(DensityFilter* filter, Projection* projection, double Smax, ProjectData* data, double rho_init, double xtol_abs, double Vfrac_abs, double result_threshold, bool save, int P, int pc):
    Smax(Smax), data(data), rho_init(rho_init), xtol_abs(xtol_abs), Vfrac_abs(Vfrac_abs), result_threshold(result_threshold), save_result(save), P(P), pc(pc), filter(filter), projection(projection), viz(nullptr), fem(nullptr), mesh(nullptr), c(1), max_V(0), cur_V(0), alpha(1), Spn(1), Sm(1), elem_number(0){}

void MinimalVolume::initialize_views(Visualization* viz){
    this->viz = viz;

    this->stress_view = viz->add_view("Von Mises Stress", ViewHandler::ViewType::ELEMENTAL, ViewHandler::DataType::STRESS);
    this->density_view = viz->add_view("Elemental Density", ViewHandler::ViewType::ELEMENTAL, ViewHandler::DataType::DENSITY);
}

TopoDS_Shape MinimalVolume::optimize(FiniteElement* fem, Meshing* mesh){

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

    this->c  = 1;
    this->grad_V = std::vector<double>(x_size, 0);
    this->alpha = 1;

    this->fem->set_steps(2);

    auto V_it = this->grad_V.begin();
    for(const auto& g:this->mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                *V_it = e->get_volume(data->thickness);
                this->max_V += *V_it;
                ++V_it;
            }
        }
    }

    this->cur_V = this->rho_init*this->max_V;

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
	for (iter = 0; (ch > this->xtol_abs && std::abs(ff-fnew)/this->max_V > this->Vfrac_abs) || std::abs(this->Sm/(this->c*this->Spn) - 1) > 1e-1; ++iter){

        update_c();

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
        logger::quick_log("Volume change: ", std::abs(ff-fnew)/this->max_V);
        logger::quick_log("Stress difference: ", std::abs(this->Sm/(this->c*this->Spn) - 1));
	}
    /* GCMMA
	for (int iter = 0; (ch > this->xtol_abs && std::abs(ff-fnew)/this->max_V > this->Vfrac_abs) || std::abs(this->Sm/(this->c*this->Spn) - 1) > 1e-4; ++iter){

        update_c();

        ff = this->fobj_grad(x, df);
        g[0] = this->fc_norm_grad(x, dg);

        // GCMMA version
        gcmma.OuterUpdate(xnew.data(), x.data(), ff, df.data(),
            g.data(), dg.data(), xmin.data(), xmax.data());

        // Check conservativity
        fnew = this->fobj(xnew);
        gnew[0] = this->fc_norm(xnew);
        logger::quick_log("Candidate: ", fnew, gnew[0]);//+this->Smax);
        logger::quick_log("");

        bool conserv = gcmma.ConCheck(fnew, gnew.data());

        int inneriter = 0;
        while(!conserv && inneriter < max_innerit){// && ch > this->xtol_abs 
            // Inner iteration update
            gcmma.InnerUpdate(xnew.data(), fnew, gnew.data(), x.data(), ff,
                df.data(), g.data(), dg.data(), xmin.data(), xmax.data());

            // Check conservativity
            fnew = this->fobj(xnew);
            gnew[0] = this->fc_norm(xnew);

            logger::quick_log("Candidate: ", fnew, gnew[0]);//+this->Smax);
            logger::quick_log("");

            conserv = gcmma.ConCheck(fnew, gnew.data());
            ++inneriter;
        }
        if(inneriter < max_innerit){// || ch > this->xtol_abs){
            logger::quick_log("Candidate is conservative.");
        } else {
            logger::quick_log("Maximum inner iterations.");
        }
        ch = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            ch = std::max(ch, std::abs(xnew[i] - x[i]));
            x[i] = xnew[i];
        }

        logger::quick_log("Inner converged to: ", fnew, gnew[0]);//+this->Smax);

        logger::quick_log("Design var change: ", ch);
        logger::quick_log("Volume change: ", std::abs(ff-fnew)/this->max_V);
        logger::quick_log("Stress difference: ", std::abs(this->Sm/(this->c*this->Spn) - 1));
	}
    */
    
    auto stop_to = std::chrono::high_resolution_clock::now();
    logger::quick_log("Final volume: ", this->cur_V);
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


double MinimalVolume::fobj(const std::vector<double>& x){
    double V = 0;

    #pragma omp parallel for reduction(+:V)
    for(size_t i = 0; i < x.size(); ++i){
        V += x[i]*this->grad_V[i];
    }

    this->cur_V = V;
    return V;
}
double MinimalVolume::fobj_grad(const std::vector<double>& x, std::vector<double>& grad){
    double V = 0;

    #pragma omp parallel for reduction(+:V)
    for(size_t i = 0; i < x.size(); ++i){
        grad[i] = this->grad_V[i];
        V += x[i]*this->grad_V[i];
    }

    this->cur_V = V;
    return V;

}

void MinimalVolume::update_c(){
    double new_c = this->Sm/this->Spn;
    if(this->c == 0){
        this->c = new_c;
        return;
    }
    double old_c = this->c;

    this->c = this->alpha*new_c + (1 - this->alpha)*this->c;
    double T = 1;
    this->alpha = std::min({std::pow(old_c/new_c,T), std::pow(new_c/old_c, T)});
}

double MinimalVolume::fc_norm(const std::vector<double>& x){
    (void)x;
    double pc = this->pc;
    double pt = 1.0/2;

    std::vector<double> u(this->fem->calculate_displacements(this->mesh, this->mesh->load_vector, x, pc));

    // Calculating global stress
    int P = this->P;
    double Spn = 0;
    double Smax = 0;

    auto x_it = x.begin();
    auto v_it = this->grad_V.begin();
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
                        if(Se > Smax){
                            Smax = Se;
                        }
                        double v = *v_it;
                        Spn += v*std::pow(Se, P);
                        ++x_it;
                        ++v_it;
                    }
                } else {
                    for(const auto& e:g->mesh){
                        const auto D = g->get_D_topopt(*x_it, 1);
                        double S = e->get_stress_at(D, e->get_centroid(), u);
                        // No relaxation when there is no void?
                        // Will do for now
                        double Se = S;
                        if(Se > Smax){
                            Smax = Se;
                        }
                        double v = *v_it;
                        Spn += v*std::pow(Se, P);
                        ++x_it;
                        ++v_it;
                    }
                }
            } else {
                logger::log_assert(num_den == 1, logger::ERROR, "minimal volume problems that require more than 1 design variables are currently not supported.");
            }
        }
    }

    Spn = std::pow(Spn, 1.0/P);

    double result = this->c*Spn;

    return result - this->Smax;
}
double MinimalVolume::fc_norm_grad(const std::vector<double>& x, std::vector<double>& grad){

    double pc = this->pc;
    double pt = 1.0/2;

    std::vector<double> u(this->fem->calculate_displacements(this->mesh, this->mesh->load_vector, x, pc));

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

    double result = this->c*Spn;
    // if(result < Smax){
    //     this->c = Smax/Spn;
    //     result = Smax;
    //     this->alpha = 0.5;
    // }

    this->Spn = Spn;
    this->Sm = Smax;

    double Sg = this->c*std::pow(Spn, 1 - P);

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

    logger::quick_log(result, this->c, Spn, Smax, this->Smax, this->alpha, this->cur_V);

    return result - this->Smax;
}

}
