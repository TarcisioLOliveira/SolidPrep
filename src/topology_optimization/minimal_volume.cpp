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


MinimalVolume::MinimalVolume(double r_o, double Smax, ProjectData* data, double rho_init, double xtol_abs, double Vfrac_abs, double result_threshold, bool save, int P, int pc):
    r_o(r_o), Smax(Smax), data(data), rho_init(rho_init), xtol_abs(xtol_abs), Vfrac_abs(Vfrac_abs), result_threshold(result_threshold), save_result(save), P(P), pc(pc), viz(nullptr), fem(nullptr), mesh(nullptr), c(1), new_x(), max_V(0), cur_V(0), alpha(1), neighbors(), p(), w(), Spn(1), Sm(1){}


TopoDS_Shape MinimalVolume::optimize(Visualization* viz, FiniteElement* fem, Meshing* mesh){

    logger::quick_log("Preparing for optimization...");

    this->viz = viz;
    this->fem = fem;
    this->mesh = mesh;
    this->c  = 1;
    this->new_x = std::vector<double>(mesh->element_list.size(), this->rho_init);
    this->grad_V = std::vector<double>(mesh->element_list.size(), 0);
    this->alpha = 1;
    this->neighbors = std::vector<std::vector<size_t>>(mesh->element_list.size()), std::vector<double>();
    this->p = std::vector<double>(mesh->element_list.size()*3);
    this->w = std::vector<double>(mesh->element_list.size());

    
    for(size_t i = 0; i < mesh->element_list.size(); ++i){
        gp_Pnt c = mesh->element_list[i]->get_centroid();
        this->p[3*i] = c.X();
        this->p[3*i+1] = c.Y();
        this->p[3*i+2] = c.Z();
    }

    // Uses more memory but is much faster
    for(size_t i = 0; i < mesh->element_list.size(); ++i){
        this->grad_V[i] = mesh->element_list[i]->get_volume(data->thickness);
        this->max_V += this->grad_V[i];
        for(size_t j = i; j < mesh->element_list.size(); ++j){
            // double dist = data.mesh->element_list[i]->get_centroid().Distance(data.mesh->element_list[j]->get_centroid());
            double dist = std::sqrt(std::pow(this->p[3*i] - this->p[3*j], 2) + std::pow(this->p[3*i+1] - this->p[3*j+1], 2) + std::pow(this->p[3*i+2] - this->p[3*j+2], 2));
            if(dist <= this->r_o){
                this->neighbors[i].push_back(j);
                this->neighbors[j].push_back(i);
                double wj = 1 - dist/this->r_o;
                this->w[i] += wj;
                this->w[j] += wj;
            }
        }
    }

    this->cur_V = this->max_V;

    std::vector<double> x = this->new_x;

    logger::quick_log("Done.");
    auto start_to = std::chrono::high_resolution_clock::now();
    
    //optimization::GCMMASolver gcmma(x.size(), 1, 0, 1e6, 1); //1e5
    optimization::MMASolver mma(x.size(), 1, 0, 1e4, 1); //1e5
    mma.SetAsymptotes(0.05, 0.2, 1.01);

    double ff;
    std::vector<double> df(x.size());
    std::vector<double> dg(x.size());
    std::vector<double> g(1);

    // std::vector<double> xnew(x);
    std::vector<double> xold(x);

    std::vector<double> xmin;
    std::vector<double> xmax;

    xmin = std::vector<double>(x.size(), 0.0);
    xmax = std::vector<double>(x.size(), 1.0);

    double fnew = this->max_V;
    std::vector<double> gnew(1);
    g[0] = 1e3;

    ff = this->fobj_grad(x, df);
    g[0] = this->fc_norm_grad(x, dg);

    // int max_innerit = 30;
    double ch = 1.0;
    int iter;
	for (iter = 0; (ch > this->xtol_abs && std::abs(ff-fnew)/this->max_V > this->Vfrac_abs) || std::abs(this->Sm/(this->c*this->Spn) - 1) > 1e-1; ++iter){

        update_c();

        fnew = ff;
        gnew = g;

        mma.Update(x.data(), df.data(), g.data(), dg.data(), xmin.data(), xmax.data());

        ch = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            ch = std::max(ch, std::abs(xold[i] - x[i]));
            xold[i] = x[i];
        }

        ff = this->fobj_grad(x, df);
        g[0] = this->fc_norm_grad(x, dg);

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
        for(size_t i = 0; i < this->new_x.size(); ++i){
            if(this->new_x[i] >= this->result_threshold){
                result = utils::cut_shape(result, mesh->element_list[i]->get_shape());
            }
            double pc = i/(double)(this->new_x.size()-1);
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

    // Density filtering
    for(size_t i = 0; i < x.size(); ++i){
        this->new_x[i] = 0;
        for(const auto& j:this->neighbors[i]){
             double dist = std::sqrt(std::pow(this->p[3*i] - this->p[3*j], 2) + std::pow(this->p[3*i+1] - this->p[3*j+1], 2) + std::pow(this->p[3*i+2] - this->p[3*j+2], 2));
             double wj = 1 - dist/this->r_o;
             this->new_x[i] += wj*x[j];
        }
        this->new_x[i] /= this->w[i];
        V += this->new_x[i]*this->grad_V[i];
    }

    this->cur_V = V;
    return V;
}
double MinimalVolume::fobj_grad(const std::vector<double>& x, std::vector<double>& grad){
    double V = 0;

    // Density filtering
    for(size_t i = 0; i < x.size(); ++i){
        grad[i] = 0;
        this->new_x[i] = 0;
        for(const auto& j:this->neighbors[i]){
             double dist = std::sqrt(std::pow(this->p[3*i] - this->p[3*j], 2) + std::pow(this->p[3*i+1] - this->p[3*j+1], 2) + std::pow(this->p[3*i+2] - this->p[3*j+2], 2));
             double wj = 1 - dist/this->r_o;
             this->new_x[i] += wj*x[j];
             // Yes, it is w[j] instead of w[i] here
             grad[i] += wj*this->grad_V[j]/this->w[j];
        }
        this->new_x[i] /= this->w[i];
        V += this->new_x[i]*this->grad_V[i];
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
    double pc = this->pc;
    double pt = 1.0/2;

    std::vector<double> u(this->fem->calculate_displacements(this->data, this->mesh, this->new_x, pc));

    // Calculating global stress
    int P = this->P;
    double Spn = 0;
    double Smax = 0;

    const auto D = this->data->geometries[0]->get_D();
    std::vector<double> stress_list(this->mesh->element_list.size());
    for(size_t i = 0; i < this->mesh->element_list.size(); ++i){
        auto& e = this->mesh->element_list[i];
        double S = e->get_stress_at(D, e->get_centroid(), u);
        double Se = std::pow(this->new_x[i], pt)*S;
        if(Se > Smax){
            Smax = Se;
        }
        double v = this->grad_V[i];
        Spn += v*std::pow(Se, P);
    }

    Spn = std::pow(Spn, 1.0/P);

    double result = this->c*Spn;

    return result - this->Smax;
}
double MinimalVolume::fc_norm_grad(const std::vector<double>& x, std::vector<double>& grad){

    double pc = this->pc;
    double pt = 1.0/2;

    std::vector<double> u(this->fem->calculate_displacements(this->data, this->mesh, this->new_x, pc));

    // Calculating stresses
    std::vector<double> fl(u.size(), 0);

    // Calculating global stress
    int P = this->P;
    double Spn = 0;
    double Smax = 0;

    const auto D = this->data->geometries[0]->get_D();
    std::vector<double> stress_list(this->mesh->element_list.size());
    for(size_t i = 0; i < this->mesh->element_list.size(); ++i){
        auto& e = this->mesh->element_list[i];
        double S = e->get_stress_at(D, e->get_centroid(), u);
        double Se = std::pow(this->new_x[i], pt)*S;
        stress_list[i] = this->new_x[i]*S;//std::pow(this->new_x[i],pc)*S;//Se;
        if(Se > Smax){
            Smax = Se;
        }
        double v = this->grad_V[i];
        Spn += v*std::pow(Se, P);

        e->get_virtual_load(D, v*std::pow(this->new_x[i], pt*P)*std::pow(S, P-2), e->get_centroid(), u, fl);
    }
    this->viz->update_stress_view(stress_list);
    this->viz->update_density_view(this->new_x);

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
    auto l = this->fem->calculate_displacements(this->data, this->mesh, this->new_x, pc, true, fl);
    logger::quick_log("} Done.");

    logger::quick_log("Calculating stress gradient...");
    std::vector<double> grad_tmp(grad);
    for(size_t i = 0; i < this->mesh->element_list.size(); ++i){
        auto& e = this->mesh->element_list[i];
        double lKu = pc*std::pow(this->new_x[i], pc-1)*e->get_compliance(D, this->data->thickness, u, l);
        double v = this->grad_V[i];
        double S = e->get_stress_at(D, e->get_centroid(), u);
        double Se = pt*v*std::pow(this->new_x[i], pt*P-1)*std::pow(S, P);

        grad_tmp[i] = Sg*(Se - lKu);
    }
    // Sensitivity filtering (implied by the density filtering)
    for(size_t i = 0; i < x.size(); ++i){
         grad[i] = 0;
         for(const auto& j:this->neighbors[i]){
             double dist = std::sqrt(std::pow(this->p[3*i] - this->p[3*j], 2) + std::pow(this->p[3*i+1] - this->p[3*j+1], 2) + std::pow(this->p[3*i+2] - this->p[3*j+2], 2));
             double wj = 1 - dist/this->r_o;
             // Yes, it is w[j] instead of w[i] here
             grad[i] += wj*grad_tmp[j]/this->w[j];
         }
     }
    logger::quick_log("Done.");

    logger::quick_log(result, this->c, Spn, Smax, this->Smax, this->alpha, this->cur_V);

    return result - this->Smax;
}

}
