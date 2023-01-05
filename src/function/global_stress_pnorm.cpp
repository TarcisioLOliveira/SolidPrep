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

#include "function/global_stress_pnorm.hpp"
#include <algorithm>

namespace function{

GlobalStressPnorm::GlobalStressPnorm(const Meshing* const mesh, FiniteElement* fem, double pc, double P, double pt):
    mesh(mesh), fem(fem), pc(pc), P(P), pt(pt){}

double GlobalStressPnorm::calculate(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x){
    (void)x;
    auto grad_V = op->get_volumes();
    auto stresses = op->get_stresses();

    // Calculating global stress
    double Spn = 0;
    double Smax = 0;

    auto stress_it = stresses.cbegin();
    auto x_it = x.cbegin();
    auto v_it = grad_V.cbegin();
    for(const auto& g:this->mesh->geometries){
        const size_t num_den = g->number_of_densities_needed();
        const size_t num_mat = g->number_of_materials();
        if(g->do_topopt){
            if(num_den == 1){
                if(num_mat == 1){
                    for(auto si = stress_it; si < stress_it+g->mesh.size(); ++si, ++x_it, ++v_it){
                        const double S = *stress_it;
                        const double Se = std::pow(*x_it, pt)*S;
                        if(Se > Smax){
                            Smax = Se;
                        }
                        const double v = *v_it;
                        Spn += v*std::pow(Se, P);
                    }
                    stress_it += g->mesh.size();
                } else {
                    logger::log_assert(num_mat == 1, logger::ERROR, "problems that require more than 1 materials are currently not supported.");
                }
            } else {
                logger::log_assert(num_den == 1, logger::ERROR, "problems that require more than 1 design variables are currently not supported.");
            }
        } else {
            for(auto si = stress_it; si < stress_it+g->mesh.size(); ++si, ++v_it){
                const double S = *stress_it;
                const double Se = S;
                if(Se > Smax){
                    Smax = Se;
                }
                const double v = *v_it;
                Spn += v*std::pow(Se, P);
            }
            stress_it += g->mesh.size();
        }
    }

    Spn = std::pow(Spn, 1.0/P);

    return Spn;
}

double GlobalStressPnorm::calculate_with_gradient(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
    // Calculating stresses
    std::vector<double> fl(u.size(), 0);

    auto grad_V = op->get_volumes();
    auto stresses = op->get_stresses();

    // Calculating global stress
    double Spn = 0;
    double Smax = 0;

    auto stress_it = stresses.cbegin();
    auto x_it = x.cbegin();
    auto v_it = grad_V.cbegin();
    for(const auto& g:this->mesh->geometries){
        const size_t num_den = g->number_of_densities_needed();
        const size_t num_mat = g->number_of_materials();
        if(g->do_topopt){
            if(num_den == 1){
                if(num_mat == 1){
                    for(auto si = stress_it; si < stress_it+g->mesh.size(); ++si, ++x_it, ++v_it){
                        const double S = *stress_it;
                        const double Se = std::pow(*x_it, pt)*S;
                        if(Se > Smax){
                            Smax = Se;
                        }
                        const double v = *v_it;
                        Spn += v*std::pow(Se, P);
                    }
                    stress_it += g->mesh.size();
                } else {
                    logger::log_assert(num_mat == 1, logger::ERROR, "problems that require more than 1 materials are currently not supported.");
                }
            } else {
                logger::log_assert(num_den == 1, logger::ERROR, "problems that require more than 1 design variables are currently not supported.");
            }
        } else {
            for(auto si = stress_it; si < stress_it+g->mesh.size(); ++si, ++v_it){
                const double S = *stress_it;
                const double Se = S;
                if(Se > Smax){
                    Smax = Se;
                }
                const double v = *v_it;
                Spn += v*std::pow(Se, P);
            }
            stress_it += g->mesh.size();
        }
    }

    Spn = std::pow(Spn, 1.0/P);

    double result = Spn;

    double Sg = std::pow(Spn, 1 - P);

    logger::quick_log("Calculating adjoint problem...{");
    auto l = this->fem->calculate_displacements(this->mesh, fl, x, pc);
    logger::quick_log("} Done.");

    logger::quick_log("Calculating stress gradient...");

    x_it = x.cbegin();
    v_it = grad_V.cbegin();
    stress_it = stresses.cbegin();
    auto grad_it = grad.begin();
    for(const auto& g:this->mesh->geometries){
        const size_t num_den = g->number_of_densities_needed();
        const size_t num_mat = g->number_of_materials();
        if(g->do_topopt){
            if(num_den == 1){
                const auto D = g->get_D(0);
                for(const auto& e:g->mesh){
                    double lKu = pc*std::pow(*x_it, pc-1)*e->get_compliance(D, this->mesh->thickness, u, l);
                    double v = *v_it;
                    double S = *stress_it;
                    double Se = pt*v*std::pow(*x_it, pt*P-1)*std::pow(S, P);

                    *grad_it = Sg*(Se - lKu);

                    ++x_it;
                    ++v_it;
                    ++grad_it;
                    ++stress_it;
                }
                if(num_mat == 2){
                    logger::log_assert(num_mat == 1, logger::ERROR, "problems that require more than 1 materials are currently not supported.");
                }
            } else {
                logger::log_assert(num_den == 1, logger::ERROR, "minimal volume problems that require more than 1 design variables are currently not supported.");
            }
        }
    }
    logger::quick_log("Done.");

    logger::quick_log(result, Spn, Smax);

    return result;
}

}

