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
#include "logger.hpp"
#include "utils.hpp"
#include <nlopt.hpp>

namespace topology_optimization{

MinimalVolume::MinimalVolume(double r_o, double Smax, ProjectData* data):
    r_o(r_o), Smax(Smax), data(data){}


TopoDS_Shape MinimalVolume::optimize(Visualization* viz, FiniteElement* fem, Meshing* mesh){
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
    };

    Data data{viz, fem, mesh, this, 1, std::vector<double>(), 0, 0, std::vector<double>(mesh->element_list.size(), 0), 1};

    std::vector<double> density(mesh->element_list.size(), 1);

    if(this->data->type == utils::PROBLEM_TYPE_2D){
        for(size_t i = 0; i < mesh->element_list.size(); ++i){
            data.grad_V[i] = mesh->element_list[i]->get_volume()*this->data->thickness*1e-3;
            data.max_V += data.grad_V[i];
        }
    } else if(this->data->type == utils::PROBLEM_TYPE_3D){
        for(size_t i = 0; i < mesh->element_list.size(); ++i){
            data.grad_V[i] = mesh->element_list[i]->get_volume();
            data.max_V += data.grad_V[i];
        }
    }
    data.cur_V = data.max_V;

    // (Le et al. 2010)
    auto f = [](const std::vector<double>& x, std::vector<double>& grad, void* f_data)->double{
        // Getting the data
        Data* data = static_cast<Data*>(f_data);

        double V = 0;
        grad = data->grad_V;

        // Density filtering
        std::vector<double> new_x(x.size(), 0);
        for(size_t i = 0; i < x.size(); ++i){
            double w = 0;
            for(size_t j = 0; j < x.size(); ++j){
                double dist = data->mesh->element_list[i]->get_centroid().Distance(data->mesh->element_list[j]->get_centroid());
                if(dist <= data->mv->r_o){
                    double wj = 1 - dist/data->mv->r_o;
                    new_x[i] += wj*x[j];
                    w += wj;
                }
            }
            new_x[i] /= w;
            V += new_x[i]*data->grad_V[i];
        }

        data->new_x = std::move(new_x);
        data->cur_V = V;

        return V;
    };
    auto fc = [](const std::vector<double>& x, std::vector<double>& grad, void* f_data)->double{
        // Getting the data
        Data* data = static_cast<Data*>(f_data);
        double pc = 3;
        double pt = 1.0/2;

        std::vector<float> u;

        std::vector<double> x_u(data->new_x);
        for(auto& d:x_u){
            d = std::pow(d, pc);
        }

        // Calculating stresses
        u = data->fem->calculate_displacements(data->mv->data, data->mesh, x_u);
        data->fem->calculate_stresses(data->mesh, u, x_u);
        data->viz->update_view();

        std::vector<float> fl(u.size(), 0);

        // Calculating global stress
        int P = 8;
        double Spn = 0;
        double Smax = 0;
        for(size_t i = 0; i < data->mesh->element_list.size(); ++i){
            auto& e = data->mesh->element_list[i];
            double S = e->get_stress_at(e->get_centroid(), u);
            double Se = S*std::pow(data->new_x[i], pt);
            if(S > Smax){
                Smax = Se;
            }
            double v = data->new_x[i]*data->grad_V[i];
            Spn += v*std::pow(Se, P);

            e->get_virtual_load(P, e->get_centroid(), u, fl);
        }

        logger::quick_log("Calculating adjoint problem...{");
        auto l = data->fem->calculate_displacements(data->mv->data, data->mesh, x_u, fl);
        logger::quick_log("} Done.");

        Spn = std::pow(Spn, 1.0/P);
        double Sg = data->c*std::pow(Spn, 1 - P)/P;

        double result = data->c*Spn;

        logger::quick_log("Calculating stress gradient...");
        for(size_t i = 0; i < data->mesh->element_list.size(); ++i){
            double uKl = 0;
            auto& e = data->mesh->element_list[i];
            if(data->mv->data->type == utils::PROBLEM_TYPE_2D){
                uKl = std::pow(data->new_x[i], pc)*e->get_compliance(u, l)*data->mv->data->thickness*1e-3;
            } else if(data->mv->data->type == utils::PROBLEM_TYPE_3D){
                uKl = std::pow(data->new_x[i], pc)*e->get_compliance(u, l);
            }
            double S = e->get_stress_at(e->get_centroid(), u);
            double Se = std::pow(S*std::pow(data->new_x[i], pt), P);
            
            grad[i] = Sg*(Se + uKl);
        }
        logger::quick_log("Done.");

        double new_c = Smax/Spn;
        double alpha_i = std::max(std::min(1.0, data->c/new_c), 0.0001);

        data->c = data->alpha*new_c + (1 - data->alpha)*data->c;

        logger::quick_log(result, data->c, Spn, Smax, data->mv->Smax, data->alpha);
        data->alpha = alpha_i;

        return result - data->mv->Smax;
    
    };
    nlopt::opt MMA(nlopt::LD_MMA, mesh->element_list.size());
    MMA.set_min_objective(f, &data);
    MMA.set_lower_bounds(0.001);
    MMA.set_upper_bounds(1);
    MMA.add_inequality_constraint(fc, &data);
    MMA.set_param("verbosity", 5);
    MMA.set_xtol_abs(0.1);

    double opt_f = 0;
    nlopt::result r = MMA.optimize(density, opt_f);

    if(r > 0){

    }

}

}
