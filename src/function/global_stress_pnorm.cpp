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

void GlobalStressPnorm::initialize(){
    size_t x_size = 0;
    for(const auto& g:this->mesh->geometries){
        if(g->do_topopt){
            x_size += g->mesh.size();
        }
        this->elem_number += g->mesh.size();
    }

    this->grad_V.resize(x_size, 0);

    auto V_it = this->grad_V.begin();
    for(const auto& g:this->mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                *V_it = e->get_volume(mesh->thickness);
                ++V_it;
            }
        }
    }
}

double GlobalStressPnorm::calculate(const std::vector<double>& u, const std::vector<double>& x){
    (void)x;

    // Calculating global stress
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

    return Spn;
}
double GlobalStressPnorm::calculate_with_gradient(const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
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
                    double lKu = pc*std::pow(*x_it, pc-1)*e->get_compliance(D, this->mesh->thickness, u, l);
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
                        double lKu = -pc*std::pow(1 - *x_it2, pc-1)*e->get_compliance(D, this->mesh->thickness, u, l);
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

    logger::quick_log(result, Spn, Smax);

    return result;
}

}

