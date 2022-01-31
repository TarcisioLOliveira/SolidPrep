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

#include "topology_optimization/minimal_volume.hpp"
#include "element/TRI3.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include <cmath>
#include <nlopt.hpp>
#include <cblas.h>
#include <BRepBuilderAPI_Copy.hxx>
#include <chrono>
#include "project_data.hpp"
#include <lapacke.h>
#include <vector>
#include "optimization/GCMMASolver.hpp"

namespace topology_optimization{


MinimalVolume::MinimalVolume(double r_o, double Smax, ProjectData* data, double rho_init, double xtol_abs, double Vfrac_abs, double result_threshold, bool save, int P, int pc):
    r_o(r_o), Smax(Smax), data(data), rho_init(rho_init), xtol_abs(xtol_abs), Vfrac_abs(Vfrac_abs), result_threshold(result_threshold), save_result(save), P(P), pc(pc){}


TopoDS_Shape MinimalVolume::optimize(Visualization* viz, FiniteElement* fem, Meshing* mesh){

    logger::quick_log("Preparing for optimization...");

    struct Data{
        Visualization* viz;
        FiniteElement* fem;
        Meshing* mesh;
        MinimalVolume* mv;
        double c;
        std::vector<double> new_x;
        std::vector<double> d;
        double max_V;
        double cur_V;
        std::vector<double> grad_V;
        double alpha;
        std::vector<std::vector<size_t>> neighbors;
        std::vector<double> u;
        int it_num;
        double xtol_abs;
        std::vector<double> p;
        std::vector<double> w;
        double Spn;
        double Sm;
        // std::vector<double> grad_obj;
        // std::vector<double> grad_S;
    };

    Data data{viz, fem, mesh, this, 1, std::vector<double>(mesh->element_list.size(), this->rho_init), std::vector<double>(mesh->element_list.size(), this->rho_init), 0, 0, std::vector<double>(mesh->element_list.size(), 0), 1, std::vector<std::vector<size_t>>(mesh->element_list.size()), std::vector<double>(), 0, this->xtol_abs, std::vector<double>(mesh->element_list.size()*3), std::vector<double>(mesh->element_list.size()), 1, 1};

    
    for(size_t i = 0; i < mesh->element_list.size(); ++i){
        gp_Pnt c = mesh->element_list[i]->get_centroid();
        data.p[3*i] = c.X();
        data.p[3*i+1] = c.Y();
        data.p[3*i+2] = c.Z();
    }

    // Uses more memory but is much faster
    for(size_t i = 0; i < mesh->element_list.size(); ++i){
        data.grad_V[i] = mesh->element_list[i]->get_volume();
        data.max_V += data.grad_V[i];
        for(size_t j = i; j < mesh->element_list.size(); ++j){
            // double dist = data.mesh->element_list[i]->get_centroid().Distance(data.mesh->element_list[j]->get_centroid());
            double dist = std::sqrt(std::pow(data.p[3*i] - data.p[3*j], 2) + std::pow(data.p[3*i+1] - data.p[3*j+1], 2) + std::pow(data.p[3*i+2] - data.p[3*j+2], 2));
            if(dist <= data.mv->r_o){
                data.neighbors[i].push_back(j);
                data.neighbors[j].push_back(i);
                double wj = 1 - dist/data.mv->r_o;
                data.w[i] += wj;
                data.w[j] += wj;
            }
        }
    }

    data.cur_V = data.max_V;

    // (Le et al. 2010)
    auto f = [](const std::vector<double>& x, std::vector<double>& grad, void* f_data)->double{
        // Getting the data
        Data* data = static_cast<Data*>(f_data);

        double V = 0;
        // double change = 0;

        // Density filtering
        for(size_t i = 0; i < x.size(); ++i){
            // change = std::max(change, std::abs(data->d[i] - x[i]));
            grad[i] = 0;
            data->new_x[i] = 0;
            for(const auto& j:data->neighbors[i]){
                 double dist = std::sqrt(std::pow(data->p[3*i] - data->p[3*j], 2) + std::pow(data->p[3*i+1] - data->p[3*j+1], 2) + std::pow(data->p[3*i+2] - data->p[3*j+2], 2));
                 double wj = 1 - dist/data->mv->r_o;
                 data->new_x[i] += wj*x[j];
                 // Yes, it is w[j] instead of w[i] here
                 grad[i] += /*(100/data->max_V)* */wj*data->grad_V[j]/data->w[j];
            }
            data->new_x[i] /= data->w[i];
            V += data->new_x[i]*data->grad_V[i];
        }
        // data->d = x;

        // logger::quick_log("");
        // logger::quick_log("Change: ", change);
        // logger::quick_log("");
        // if(data->it_num > 1 && change < data->xtol_abs){
        //     data->cur_V = V;
        //     throw nlopt::forced_stop();
        // }
        // if((data->it_num > 0) && V <= 0.0011*data->max_V){
        //     data->cur_V = V;
        //     throw nlopt::forced_stop();
        // }

        ++data->it_num;

        data->cur_V = V;
        return V;//100*V/data->max_V;
    };
    auto fc = [](const std::vector<double>& x, std::vector<double>& grad, void* f_data)->double{
        // Getting the data
        Data* data = static_cast<Data*>(f_data);

        // double pc = std::min(1+data->it_num*0.1, data->mv->pc*1.0);//data->mv->pc;
        // double pt = 1.0/std::max(1.0, 3-data->it_num*0.1);
        double pc = data->mv->pc;
        double pt = 1.0/2;

        data->u = data->fem->calculate_displacements(data->mv->data, data->mesh, data->new_x, pc);

        // Calculating stresses
        std::vector<double> fl(data->u.size(), 0);

        // Calculating global stress
        int P = data->mv->P;
        double Spn = 0;
        double Smax = 0;

        std::vector<double> stress_list(data->mesh->element_list.size());
        for(size_t i = 0; i < data->mesh->element_list.size(); ++i){
            auto& e = data->mesh->element_list[i];
            double S = e->get_stress_at(e->get_centroid(), data->u);
            double Se = std::pow(data->new_x[i], pt)*S;
            stress_list[i] = data->new_x[i]*S;//std::pow(data->new_x[i],pc)*S;//Se;
            if(Se > Smax){
                Smax = Se;
            }
            double v = data->grad_V[i];
            Spn += v*std::pow(Se, P);

            e->get_virtual_load(v*std::pow(data->new_x[i], pt*P)*std::pow(S, P-2), e->get_centroid(), data->u, fl);
        }
        data->viz->update_stress_view(stress_list);
        //data->viz->update_density_view(data->new_x);

        Spn = std::pow(Spn, 1.0/P);

        double result = data->c*Spn;
        if(result < Smax){
            data->c = Smax/Spn;
            result = Smax;
            data->alpha = 0.5;
        }
        // } else {
        //     double new_c = Smax/Spn;
        //     if(data->c == 0){
        //         data->c = new_c;
        //     }
        //     double T = 1;
        //     double alpha_i = std::min({std::pow(data->c/new_c,T), std::pow(new_c/data->c, T)});
        //     data->alpha = alpha_i;
        //     data->c = data->alpha*new_c + (1 - data->alpha)*data->c;
        // }

        data->Spn = Spn;
        data->Sm = Smax;

        double Sg = data->c*std::pow(Spn, 1 - P);

        logger::quick_log("Calculating adjoint problem...{");
        auto l = data->fem->calculate_displacements(data->mv->data, data->mesh, data->new_x, pc, fl);
        logger::quick_log("} Done.");

        logger::quick_log("Calculating stress gradient...");
        std::vector<double> grad_tmp(grad);
        for(size_t i = 0; i < data->mesh->element_list.size(); ++i){
            auto& e = data->mesh->element_list[i];
            double lKu = pc*std::pow(data->new_x[i], pc-1)*e->get_compliance(data->u, l);
            double v = data->grad_V[i];
            double S = e->get_stress_at(e->get_centroid(), data->u);
            double Se = pt*v*std::pow(data->new_x[i], pt*P-1)*std::pow(S, P);

            grad_tmp[i] = Sg*(Se - lKu);
            // grad[i] = Sg*(Se - lKu);
        }
        // Sensitivity filtering (implied by the density filtering)
        for(size_t i = 0; i < x.size(); ++i){
             grad[i] = 0;
             for(const auto& j:data->neighbors[i]){
                 double dist = std::sqrt(std::pow(data->p[3*i] - data->p[3*j], 2) + std::pow(data->p[3*i+1] - data->p[3*j+1], 2) + std::pow(data->p[3*i+2] - data->p[3*j+2], 2));
                 double wj = 1 - dist/data->mv->r_o;
                 // Yes, it is w[j] instead of w[i] here
                 grad[i] += wj*grad_tmp[j]/data->w[j];
             }
         }
        logger::quick_log("Done.");

        logger::quick_log(result, data->c, Spn, Smax, data->mv->Smax, data->alpha, data->cur_V);

        return result - data->mv->Smax;
    
    };

    auto fobj = [](const std::vector<double>& x, void* f_data)->double{
        // Getting the data
        Data* data = static_cast<Data*>(f_data);

        double V = 0;

        // Density filtering
        for(size_t i = 0; i < x.size(); ++i){
            data->new_x[i] = 0;
            for(const auto& j:data->neighbors[i]){
                 double dist = std::sqrt(std::pow(data->p[3*i] - data->p[3*j], 2) + std::pow(data->p[3*i+1] - data->p[3*j+1], 2) + std::pow(data->p[3*i+2] - data->p[3*j+2], 2));
                 double wj = 1 - dist/data->mv->r_o;
                 data->new_x[i] += wj*x[j];
            }
            data->new_x[i] /= data->w[i];
            V += data->new_x[i]*data->grad_V[i];
        }
        // Check connectivity
        // bool changed = true;
        // while(changed){
        //     changed = false;
        //     for(size_t i = 0; i < x.size(); ++i){
        //         if(data->new_x[i] < 1e-3){
        //             continue;
        //         }
        //         bool connected = false;
        //         auto& e1 = data->mesh->element_list[i];
        //         for(const auto& j:data->neighbors[i]){
        //             bool shares_node = false;
        //             auto& e2 = data->mesh->element_list[j];
        //             for(auto& n1:e1->nodes){
        //                 for(auto& n2:e2->nodes){
        //                     if(n1->id == n2->id){
        //                         shares_node = true;
        //                         break;
        //                     }
        //                 }
        //                 if(shares_node){
        //                     break;
        //                 }
        //             }
        //             if(shares_node){
        //                 if(data->new_x[j] >= 1e-3){
        //                     connected = true;
        //                     break;
        //                 }
        //             }
        //         }
        //         if(!connected){
        //             data->new_x[i] = 0;
        //             changed = true;
        //         }
        //     }
        // }


        data->cur_V = V;
        return V;//(100/data->max_V)*V;// - data->mv->rho_init*data->max_V;
    };
    auto fcobj = [](const std::vector<double>& x, void* f_data)->double{
        // Getting the data
        Data* data = static_cast<Data*>(f_data);

        // double pc = std::min(1+data->it_num*0.1, data->mv->pc*1.0);//data->mv->pc;
        // double pt = 1.0/std::max(1.0, 3-data->it_num*0.1);
        double pc = data->mv->pc;
        double pt = 1.0/2;

        data->u = data->fem->calculate_displacements(data->mv->data, data->mesh, data->new_x, pc);

        // Calculating global stress
        int P = data->mv->P;
        double Spn = 0;
        double Smax = 0;

        std::vector<double> stress_list(data->mesh->element_list.size());
        for(size_t i = 0; i < data->mesh->element_list.size(); ++i){
            auto& e = data->mesh->element_list[i];
            double S = e->get_stress_at(e->get_centroid(), data->u);
            double Se = std::pow(data->new_x[i], pt)*S;
            if(Se > Smax){
                Smax = Se;
            }
            double v = data->grad_V[i];
            Spn += v*std::pow(Se, P);
        }

        Spn = std::pow(Spn, 1.0/P);

        double result = data->c*Spn;

        return result - data->mv->Smax;
    
    };

    auto update_c = [&](){
        double new_c = data.Sm/data.Spn;
        if(data.c == 0){
            data.c = new_c;
            return;
        }
        double old_c = data.c;

        data.c = data.alpha*new_c + (1 - data.alpha)*data.c;
        // if(data.alpha < 0.95 || data.it_num <= 1){
        double T = 1;
        data.alpha = std::min({std::pow(old_c/new_c,T), std::pow(new_c/old_c, T)});
        // } else {
        //     data.alpha = 1;
        // }
        
    };

    // nlopt::opt MMA(nlopt::LD_MMA, mesh->element_list.size());
    // MMA.set_min_objective(f, &data);
    // MMA.add_inequality_constraint(fc, &data, 1e-10);
    // MMA.set_param("verbosity", 5);
    // MMA.set_xtol_abs(this->xtol_abs);
    // MMA.set_param("inner_maxeval", 15);
    // MMA.set_param("dual_algorithm", nlopt::LN_SBPLX);

    std::vector<double> x = data.new_x;

    // if(data.move_limit >= 1.0){

    // MMA.set_lower_bounds(0.001);
    // MMA.set_upper_bounds(1.0);

    // } else {

    //     std::vector<double> xmin(x.size(), std::max(this->rho_init - data.move_limit, 0.001));
    //     std::vector<double> xmax(x.size(), std::min(this->rho_init + data.move_limit, 1.0));

    //     MMA.set_lower_bounds(xmin);
    //     MMA.set_upper_bounds(xmax);

    // }

    // data.opt = &MMA;

    logger::quick_log("Done.");
    auto start_to = std::chrono::high_resolution_clock::now();

    nlopt::result r = nlopt::FORCED_STOP;
    // double opt_f = 0;
    // try{
    //     r = MMA.optimize(x, opt_f);
    // } catch(nlopt::forced_stop& f){
    //     r = nlopt::FORCED_STOP;
    // }
    
    optimization::GCMMASolver gcmma(x.size(), 1, 0, 1e6, 1); //1e5

    double ff;
    std::vector<double> df(x.size());
    std::vector<double> dg(x.size());
    std::vector<double> g(1);

    std::vector<double> xnew(x);

    std::vector<double> xmin;
    std::vector<double> xmax;

    xmin = std::vector<double>(x.size(), 0.001);
    xmax = std::vector<double>(x.size(), 1.0);

    double fnew = data.max_V;
    std::vector<double> gnew(1);
    g[0] = 1e3;

    // fnew = fobj(xnew, &data);
    // gnew[0] = fcobj(xnew, &data);

    int max_innerit = 30;
    double ch = 1.0;
	for (int iter = 0; (ch > this->xtol_abs && std::abs(ff-fnew)/data.max_V > this->Vfrac_abs) || std::abs(data.Sm/(data.c*data.Spn) - 1) > 1e-4; ++iter){

        update_c();

        ff = f(x, df, &data);
        g[0] = fc(x, dg, &data);

        // GCMMA version
        gcmma.OuterUpdate(xnew.data(), x.data(), ff, df.data(),
            g.data(), dg.data(), xmin.data(), xmax.data());

        // Check conservativity
        fnew = fobj(xnew, &data);
        gnew[0] = fcobj(xnew, &data);
        logger::quick_log("Candidate: ", fnew, gnew[0]);//+this->Smax);
        logger::quick_log("");

        bool conserv = gcmma.ConCheck(fnew, gnew.data());

        int inneriter = 0;
        while(!conserv && inneriter < max_innerit){// && ch > this->xtol_abs 
            // Inner iteration update
            gcmma.InnerUpdate(xnew.data(), fnew, gnew.data(), x.data(), ff,
                df.data(), g.data(), dg.data(), xmin.data(), xmax.data());

            // Check conservativity
            fnew = fobj(xnew, &data);
            gnew[0] = fcobj(xnew, &data);

            logger::quick_log("Candidate: ", fnew, gnew[0]);//+this->Smax);
            logger::quick_log("");

            conserv = gcmma.ConCheck(fnew, gnew.data());
            ++inneriter;

            // ch = 0.0;
            // for (size_t i = 0; i < x.size(); ++i) {
            //     ch = std::max(ch, std::abs(xnew[i] - x[i]));
            // }
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
        logger::quick_log("Volume change: ", std::abs(ff-fnew)/data.max_V);
        logger::quick_log("Stress difference: ", std::abs(data.Sm/(data.c*data.Spn) - 1));
	}
    
    auto stop_to = std::chrono::high_resolution_clock::now();
    logger::quick_log("Final volume: ", data.cur_V);
    auto to_duration = std::chrono::duration_cast<std::chrono::seconds>(stop_to-start_to);
    double to_time = to_duration.count();
    double it_time = to_time/data.it_num;
    logger::quick_log("Time per iteration (topology optimization): ", it_time, " seconds");
    logger::quick_log("Number of iterations (topology optimization): ", data.it_num);
   
    logger::quick_log(" "); 
    if(this->save_result && (r > 0 || r == nlopt::FORCED_STOP)){
        logger::quick_log("Saving resulting topology...");
        std::cout << "\r" << 0 << "%         ";
        TopoDS_Shape result = BRepBuilderAPI_Copy(this->data->ground_structure->shape);
        for(size_t i = 0; i < data.new_x.size(); ++i){
            if(data.new_x[i] >= this->result_threshold){
                result = utils::cut_shape(result, mesh->element_list[i]->get_shape());
            }
            double pc = i/(double)(data.new_x.size()-1);
            std::cout << "\r" << pc*100 << "%         ";
        }
        result = utils::cut_shape(this->data->ground_structure->shape, result);
        return result;
        std::cout << "\r" << 100 << "%         ";
        logger::quick_log(" "); 
    }

    return TopoDS_Shape();
}


}
