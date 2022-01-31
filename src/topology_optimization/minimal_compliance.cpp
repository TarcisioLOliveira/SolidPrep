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

#include "topology_optimization/minimal_compliance.hpp"
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

namespace topology_optimization{


MinimalCompliance::MinimalCompliance(double r_o, ProjectData* data, double Vfinal, double xtol_abs, double result_threshold, bool save, int pc):
    r_o(r_o), data(data), Vfinal(Vfinal), xtol_abs(xtol_abs), result_threshold(result_threshold), save_result(save), pc(pc){}


TopoDS_Shape MinimalCompliance::optimize(Visualization* viz, FiniteElement* fem, Meshing* mesh){

    logger::quick_log("Preparing for optimization...");

    struct Data{
        Visualization* viz;
        FiniteElement* fem;
        Meshing* mesh;
        MinimalCompliance* mv;
        double c;
        std::vector<double> new_x;
        std::vector<double> d;
        double max_V;
        double cur_V;
        std::vector<double> grad_V;
        double alpha;
        std::vector<std::vector<size_t>> neighbors;
        int it_num;
        double xtol_abs;
        std::vector<double> p;
        std::vector<double> w;
    };

    Data data{viz, fem, mesh, this, 1, std::vector<double>(mesh->element_list.size(), this->Vfinal), std::vector<double>(mesh->element_list.size(), Vfinal), 0, 0, std::vector<double>(mesh->element_list.size(), 0), 1, std::vector<std::vector<size_t>>(mesh->element_list.size()), 0, this->xtol_abs, std::vector<double>(mesh->element_list.size()*3), std::vector<double>(mesh->element_list.size())};

    
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
    // for(auto& v:data.grad_V){
    //     v /= data.max_V;
    // }
    // data.max_V = 1;
    data.cur_V = data.max_V;

    // (Le et al. 2010)
    auto f = [](const std::vector<double>& x, std::vector<double>& grad, void* f_data)->double{
        // Getting the data
        Data* data = static_cast<Data*>(f_data);

        double change = 0;
        double pc = data->mv->pc;
        double c = 0;

        std::vector<double> grad_tmp(x.size());

        // Density filtering
        for(size_t i = 0; i < x.size(); ++i){
            change = std::max(change, std::abs(data->d[i] - x[i]));
            data->new_x[i] = 0;
            for(const auto& j:data->neighbors[i]){
                // double dist = data->mesh->element_list[i]->get_centroid().Distance(data->mesh->element_list[j]->get_centroid());
                double dist = std::sqrt(std::pow(data->p[3*i] - data->p[3*j], 2) + std::pow(data->p[3*i+1] - data->p[3*j+1], 2) + std::pow(data->p[3*i+2] - data->p[3*j+2], 2));
                double wj = 1 - dist/data->mv->r_o;
                data->new_x[i] += wj*x[j];
            }
            data->new_x[i] /= data->w[i];
        }
        data->viz->update_density_view(data->new_x);

        // logger::quick_log("");
        logger::quick_log("Change: ", change);
        // logger::quick_log("");
        if(data->it_num > 1 && change < data->xtol_abs){
            throw nlopt::forced_stop();
        }

        std::vector<double> u = data->fem->calculate_displacements(data->mv->data, data->mesh, data->new_x, pc);

        for(size_t i = 0; i < data->new_x.size(); ++i){
            auto& e = data->mesh->element_list[i];
            double uKu = e->get_compliance(u);
            grad_tmp[i] = -pc*std::pow(data->new_x[i], pc-1)*uKu;
            c += std::pow(data->new_x[i], pc)*uKu;
        }
        for(size_t i = 0; i < x.size(); ++i){
             grad[i] = 0;
             for(const auto& j:data->neighbors[i]){
                 // double dist = data->mesh->element_list[i]->get_centroid().Distance(data->mesh->element_list[j]->get_centroid());
                 double dist = std::sqrt(std::pow(data->p[3*i] - data->p[3*j], 2) + std::pow(data->p[3*i+1] - data->p[3*j+1], 2) + std::pow(data->p[3*i+2] - data->p[3*j+2], 2));
                 double wj = 1 - dist/data->mv->r_o;
                 // Yes, it is w[j] instead of w[i] here
                 grad[i] += wj*grad_tmp[j]/data->w[j];
             }
         }

        ++data->it_num;

        data->d = x;

        return c;
    };
    auto fc = [](const std::vector<double>& x, std::vector<double>& grad, void* f_data)->double{
        // Getting the data
        Data* data = static_cast<Data*>(f_data);

        double V = 0;

        for(size_t i = 0; i < x.size(); ++i){
            grad[i] = 0;
            V += data->new_x[i]*data->grad_V[i];
            for(const auto& j:data->neighbors[i]){
                // double dist = data->mesh->element_list[i]->get_centroid().Distance(data->mesh->element_list[j]->get_centroid());
                double dist = std::sqrt(std::pow(data->p[3*i] - data->p[3*j], 2) + std::pow(data->p[3*i+1] - data->p[3*j+1], 2) + std::pow(data->p[3*i+2] - data->p[3*j+2], 2));
                double wj = 1 - dist/data->mv->r_o;
                // Yes, it is w[j] instead of w[i] here
                grad[i] += wj*data->grad_V[j]/data->w[j];
            }
        }
        data->cur_V = V;

        if(data->it_num > 1 && V < 0.0011*data->max_V){
            data->cur_V = V;
            throw nlopt::forced_stop();
        }

        return V - data->mv->Vfinal*data->max_V;
    };

    nlopt::opt MMA(nlopt::LD_MMA, mesh->element_list.size());
    MMA.set_min_objective(f, &data);
    MMA.set_lower_bounds(0.001);
    MMA.set_upper_bounds(1);
    MMA.add_inequality_constraint(fc, &data, 1e-10);
    MMA.set_param("verbosity", 5);
    MMA.set_xtol_abs(this->xtol_abs);
    //MMA.set_param("inner_maxeval", 5);

    double opt_f = 0;

    std::vector<double> x = data.new_x;

    logger::quick_log("Done.");
    auto start_to = std::chrono::high_resolution_clock::now();

    nlopt::result r = nlopt::FORCED_STOP;
    try{
        r = MMA.optimize(x, opt_f);
    } catch(nlopt::forced_stop& f){
        r = nlopt::FORCED_STOP;
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
