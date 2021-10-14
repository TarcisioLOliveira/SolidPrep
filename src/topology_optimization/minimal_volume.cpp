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
#include <nlopt.hpp>
#include <cblas.h>
#include <BRepBuilderAPI_Copy.hxx>
#include <chrono>
#include "project_data.hpp"

namespace topology_optimization{

MinimalVolume::MinimalVolume(double r_o, double Smax, ProjectData* data, double rho_init, double ftol_rel, double result_threshold, bool save):
    r_o(r_o), Smax(Smax), data(data), rho_init(rho_init), ftol_rel(ftol_rel), result_threshold(result_threshold), save_result(save){}


TopoDS_Shape MinimalVolume::optimize(Visualization* viz, FiniteElement* fem, Meshing* mesh){

    logger::quick_log("Preparing for optimization...");

    struct Data{
        Visualization* viz;
        FiniteElement* fem;
        Meshing* mesh;
        MinimalVolume* mv;
        double c;
        std::vector<double> new_x;
        double max_V;
        double cur_V;
        std::vector<double> grad_V;
        double alpha;
        std::vector<std::vector<size_t>> neighbors;
        std::vector<double> u;
        int it_num;
        double ftol_rel;
    };

    Data data{viz, fem, mesh, this, 1, std::vector<double>(mesh->element_list.size(), this->rho_init), 0, 0, std::vector<double>(mesh->element_list.size(), 0), 1, std::vector<std::vector<size_t>>(mesh->element_list.size()), std::vector<double>(), 0, this->ftol_rel};

    // Uses more memory but is much faster
    for(size_t i = 0; i < mesh->element_list.size(); ++i){
        data.grad_V[i] = mesh->element_list[i]->get_volume();
        data.max_V += data.grad_V[i];
        for(size_t j = i; j < mesh->element_list.size(); ++j){
            double dist = data.mesh->element_list[i]->get_centroid().Distance(data.mesh->element_list[j]->get_centroid());
            if(dist <= data.mv->r_o){
                data.neighbors[i].push_back(j);
                data.neighbors[j].push_back(i);
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

        double V = 0;
        double pc = 3;

        // Density filtering
        for(size_t i = 0; i < x.size(); ++i){
            double w = 0;
            data->new_x[i] = 0;
            for(const auto& j:data->neighbors[i]){
                double dist = data->mesh->element_list[i]->get_centroid().Distance(data->mesh->element_list[j]->get_centroid());
                double wj = 1 - dist/data->mv->r_o;
                data->new_x[i] += wj*x[j];
                w += wj;
            }
            data->new_x[i] /= w;
            V += data->new_x[i]*data->grad_V[i];
        }
        data->u = data->fem->calculate_displacements(data->mv->data, data->mesh, data->new_x, pc);
        grad = data->grad_V;
        //for(size_t i = 0; i < data->new_x.size(); ++i){
        //    auto& e = data->mesh->element_list[i];
        //    grad[i] = -pc*std::pow(data->new_x[i], pc-1)*e->get_compliance(data->u);
        //}

        ++data->it_num;

        if(data->it_num > 1 && std::abs(V - data->cur_V)/V < data->ftol_rel){
            data->cur_V = V;
            throw nlopt::forced_stop();
        }

        data->cur_V = V;
        return V;
        //return cblas_ddot(data->u.size(), data->u.data(), 1, data->mesh->load_vector.data(), 1);
    };
    auto fc = [](const std::vector<double>& x, std::vector<double>& grad, void* f_data)->double{
        // Getting the data
        Data* data = static_cast<Data*>(f_data);
        double pc = 3;
        double pt = 1.0/2;


        // Calculating stresses
        std::vector<double> fl(data->u.size(), 0);

        // Calculating global stress
        int P = 20;
        double Spn = 0;
        double Smax = 0;

        std::vector<double> stress_list(data->mesh->element_list.size());
        for(size_t i = 0; i < data->mesh->element_list.size(); ++i){
            auto& e = data->mesh->element_list[i];
            double S = e->get_stress_at(e->get_centroid(), data->u);
            double Se = std::pow(data->new_x[i], pt)*S;
            stress_list[i] = Se;
            if(Se > Smax){
                Smax = Se;
            }
            double v = data->new_x[i]*data->grad_V[i];
            Spn += v*std::pow(Se, P);

            e->get_virtual_load(P*v*std::pow(data->new_x[i], pt*P)*std::pow(S, P-2), e->get_centroid(), data->u, fl);
            //e->get_virtual_load(v*std::pow(data->new_x[i], pt)/(S), e->get_centroid(), u, fl);
        }
        data->viz->update_stress_view(stress_list);
        //data->viz->update_density_view(data->new_x);

        Spn = std::pow(Spn, 1.0/P);
        double new_c = Smax/Spn;
        if(data->c == 0){
            data->c = new_c;
        }
        double Sg = data->c*std::pow(Spn, 1 - P)/P;

        double result = data->c*Spn;

        logger::quick_log("Calculating adjoint problem...{");
        auto l = data->fem->calculate_displacements(data->mv->data, data->mesh, data->new_x, pc, fl);
        logger::quick_log("} Done.");

        logger::quick_log("Calculating stress gradient...");
        // std::vector<double> grad_tmp(grad);
        for(size_t i = 0; i < data->mesh->element_list.size(); ++i){
            auto& e = data->mesh->element_list[i];
            double lKu = pc*std::pow(data->new_x[i], pc-1)*e->get_compliance(data->u, l);
            double v = data->new_x[i]*data->grad_V[i];
            double S = e->get_stress_at(e->get_centroid(), data->u);
            double Se = (pt*P+1)*v*std::pow(data->new_x[i], pt*P-1)*std::pow(S, P);
            //double Se = (pt+1)*v*std::pow(data->new_x[i], pt-1)*S;

            //grad_tmp[i] = Sg*(Se - lKu);
            grad[i] = Sg*(Se - lKu);
            // grad[i] = Se - lKu;
        }
        // Sensitivity filtering
        // for(size_t i = 0; i < x.size(); ++i){
        //     double w = 0;
        //     grad[i] = 0;
        //     for(const auto& j:data->neighbors[i]){
        //         double dist = data->mesh->element_list[i]->get_centroid().Distance(data->mesh->element_list[j]->get_centroid());
        //         double wj = 1 - dist/data->mv->r_o;
        //         grad[i] += wj*data->new_x[j]*grad_tmp[j];
        //         w += wj;
        //     }
        //     grad[i] /= w*data->new_x[i];
        // }
        logger::quick_log("Done.");

        double T = 1;
        double alpha_i = std::min({1.0, std::pow(data->c/new_c,T), std::pow(new_c/data->c, T)});
        data->alpha = alpha_i;

        data->c = data->alpha*new_c + (1 - data->alpha)*data->c;

        logger::quick_log(result, data->c, Spn, Smax, data->mv->Smax, data->alpha, data->cur_V);

        return result;
    
    };
    nlopt::opt MMA(nlopt::LD_MMA, mesh->element_list.size());
    MMA.set_min_objective(f, &data);
    MMA.set_lower_bounds(0.001);
    MMA.set_upper_bounds(1);
    MMA.add_inequality_constraint(fc, &data, this->Smax);
    MMA.set_param("verbosity", 5);
    MMA.set_ftol_rel(data.ftol_rel);

    double opt_f = 0;

    std::vector<double> density = data.new_x;

    logger::quick_log("Done.");
    auto start_to = std::chrono::high_resolution_clock::now();
    nlopt::result r = nlopt::FORCED_STOP;
    try{
        r = MMA.optimize(density, opt_f);
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
    logger::quick_log("Saving resulting topology...");
    std::cout << "\r" << 0 << "%         ";
    if(this->save_result && (r > 0 || r == nlopt::FORCED_STOP)){
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
    }
    std::cout << "\r" << 100 << "%         ";
    logger::quick_log(" "); 

    return TopoDS_Shape();
}

}
