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

namespace topology_optimization{

MinimalVolume::MinimalVolume(double r_o, double Smax, ProjectData* data):
    r_o(r_o), Smax(Smax), data(data){}


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
    };

    Data data{viz, fem, mesh, this, 1, std::vector<double>(mesh->element_list.size(), 0.3), 0, 0, std::vector<double>(mesh->element_list.size(), 0), 1, std::vector<std::vector<size_t>>(mesh->element_list.size()), std::vector<double>()};

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
        data->cur_V = V;
        data->u = data->fem->calculate_displacements(data->mv->data, data->mesh, data->new_x, pc);
        for(size_t i = 0; i < data->new_x.size(); ++i){
            auto& e = data->mesh->element_list[i];
            grad[i] = -pc*std::pow(data->new_x[i], pc-1)*e->get_compliance(data->u);
        }

        // return V;
        return cblas_ddot(data->u.size(), data->u.data(), 1, data->mesh->load_vector.data(), 1);
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
        // data->viz->update_density_view(data->new_x);

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

            // grad_tmp[i] = Sg*(Se - lKu);
            grad[i] = Sg*(Se - lKu);
            // grad[i] = Se - lKu;
        }
        // // Sensitivity filtering
        // for(size_t i = 0; i < x.size(); ++i){
        //     double w = 0;
        //     grad[i] = 0;
        //     for(const auto& j:data->neighbors[i]){
        //         double dist = data->mesh->element_list[i]->get_centroid().Distance(data->mesh->element_list[j]->get_centroid());
        //         double wj = 1 - dist/data->mv->r_o;
        //         grad[i] += wj*grad_tmp[j];
        //         w += wj;
        //     }
        //     grad[i] /= w;
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
    MMA.set_xtol_rel(1e-6);

    double opt_f = 0;

    std::vector<double> density = data.new_x;

    logger::quick_log("Done.");
    nlopt::result r = MMA.optimize(density, opt_f);
    logger::quick_log("Final volume: ", data.cur_V);

    if(r > 0){

    }

    return TopoDS_Shape();
}

}
