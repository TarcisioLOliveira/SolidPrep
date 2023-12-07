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

#include "function/global_stress_pnorm_normalized.hpp"
#include <algorithm>
#include <mpich-x86_64/mpi.h>

namespace function{

GlobalStressPnormNormalized::GlobalStressPnormNormalized(const Meshing* const mesh, FiniteElement* fem, double pc, double P, double pt, double psiK, double psiS):
    mesh(mesh), fem(fem), pc(pc), P(P), pt(pt), psiK(psiK), psiS(psiS){}

void GlobalStressPnormNormalized::update(){
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

double GlobalStressPnormNormalized::calculate(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x){
    auto grad_V = op->get_volumes();
    //auto stresses = op->get_stresses();
    std::vector<double> stresses(grad_V.size());
    auto s_it = stresses.begin();
    auto x_it0 = x.begin();
    for(const auto& g:this->mesh->geometries){
        g->get_stresses(u, pt, this->psiS, x_it0, s_it);
    }

    // Calculating global stress
    double Spn = 0;
    double Smax = 0;

    auto stress_it = stresses.cbegin();
    auto x_it = x.cbegin();
    auto v_it = grad_V.cbegin();
    for(auto si = stress_it; si < stresses.end(); ++si, ++x_it, ++v_it){
        const double Se = *stress_it;
        if(Se > Smax){
            Smax = Se;
        }
        const double v = *v_it;
        Spn += v*std::pow(Se, P);
    }

    Spn = std::pow(Spn, 1.0/P);

    double result = this->c*Spn;

    return result;
}
double GlobalStressPnormNormalized::calculate_with_gradient(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    // Calculating stresses
    std::vector<double> fl(u.size(), 0);

    auto grad_V = op->get_volumes();
    //auto stresses = op->get_stresses();
    std::vector<double> stresses(grad_V.size());
    auto s_it = stresses.begin();
    auto x_it0 = x.begin();
    for(const auto& g:this->mesh->geometries){
        g->get_stresses(u, 0.0, this->psiS, x_it0, s_it);
    }

    // Calculating global stress
    double Spn = 0;
    double Smax = 0;
    double result = 0;
    double Sg = 0;
    const size_t s_size = this->mesh->elem_info->get_D_dimension();

    auto stress_it = stresses.cbegin();
    auto x_it = x.cbegin();
    auto v_it = grad_V.cbegin();

    if(mpi_id == 0){
        for(const auto& g:this->mesh->geometries){
            if(g->do_topopt){
                const size_t num_den = g->number_of_densities_needed();
                auto D = std::vector<double>(s_size*s_size, 0);
                if(g->with_void){
                    for(const auto& e:g->mesh){
                        const auto c = e->get_centroid();
                        g->materials.get_D(x_it, this->psiS, e.get(), c, D);
                        const double S = *stress_it;
                        const double Se = std::pow(*x_it, pt)*S;
                        if(Se > Smax){
                            Smax = Se;
                        }
                        const double v = *v_it;
                        e->get_virtual_load(D, v*std::pow(*x_it, pt*P)*std::pow(S, P-2), c, u, fl);
                        Spn += v*std::pow(Se, P);

                        x_it += num_den;
                        ++v_it;
                        ++stress_it;
                    }
                } else {
                    for(const auto& e:g->mesh){
                        const auto c = e->get_centroid();
                        g->materials.get_D(x_it, this->psiS, e.get(), c, D);
                        const double Se = *stress_it;
                        if(Se > Smax){
                            Smax = Se;
                        }
                        const double v = *v_it;
                        e->get_virtual_load(D, v*std::pow(Se, P-2), c, u, fl);
                        Spn += v*std::pow(Se, P);

                        x_it += num_den;
                        ++v_it;
                        ++stress_it;
                    }
                }
            } else {
                for(const auto& e:g->mesh){
                    const auto c = e->get_centroid();
                    const auto D = g->materials.get_D(e.get(), c);
                    const double S = *stress_it;
                    const double Se = S;
                    if(Se > Smax){
                        Smax = Se;
                    }
                    const double v = *v_it;
                    e->get_virtual_load(D, v*std::pow(Se, P-2), c, u, fl);
                    Spn += v*std::pow(Se, P);

                    ++v_it;
                    ++stress_it;
                }
            }
        }


        Spn = std::pow(Spn, 1.0/P);

        result = this->c*Spn;

        this->Spn = Spn;
        this->Sm = Smax;

        Sg = this->c*std::pow(Spn, 1 - P);

        logger::quick_log("Calculating adjoint problem...{");
    }
    auto l = this->fem->calculate_displacements(this->mesh, fl, x, pc);
    if(mpi_id == 0){
        logger::quick_log("} Done.");

        logger::quick_log("Calculating stress gradient...");

        x_it = x.cbegin();
        v_it = grad_V.cbegin();
        stress_it = stresses.cbegin();
        auto grad_it = grad.begin();
        for(const auto& g:this->mesh->geometries){
            const size_t num_den = g->number_of_densities_needed();
            if(g->do_topopt){
                std::vector<std::vector<double>> gradD_K(num_den, std::vector<double>(s_size*s_size, 0));
                std::vector<std::vector<double>> gradD_S(num_den, std::vector<double>(s_size*s_size, 0));
                std::vector<double> D_S(std::vector<double>(s_size*s_size, 0));
                if(g->with_void){
                    for(const auto& e:g->mesh){
                        const auto c = e->get_centroid();
                        g->materials.get_gradD(x_it, psiK, e.get(), c, gradD_K);
                        g->materials.get_gradD(x_it, psiS, e.get(), c, gradD_S);
                        g->materials.get_D(x_it, psiS, e.get(), c, D_S);

                        const double lKu = pc*std::pow(*x_it, pc-1)*e->get_compliance(gradD_K[0], this->mesh->thickness, u, l);
                        const double v = *v_it;
                        const double S = *stress_it;
                        const double Se = pt*v*std::pow(*x_it, pt*P-1)*std::pow(S, P);

                        *grad_it = Sg*(Se - lKu);

                        ++x_it;
                        ++grad_it;

                        for(size_t i = 1; i < num_den; ++i){
                            const double lKu = std::pow(*x_it, pc)*e->get_compliance(gradD_K[i], this->mesh->thickness, u, l);
                            const double Se = e->von_Mises_derivative(D_S, gradD_S[i], v*std::pow(*x_it, pt*P)*std::pow(S, P-2), c, u);

                            *grad_it = Sg*(Se - lKu);

                            ++x_it;
                            ++grad_it;
                        }

                        ++v_it;
                        ++stress_it;
                    }
                } else {
                    for(const auto& e:g->mesh){
                        const auto c = e->get_centroid();
                        g->materials.get_gradD(x_it, psiK, e.get(), c, gradD_K);
                        g->materials.get_gradD(x_it, psiS, e.get(), c, gradD_S);
                        g->materials.get_D(x_it, psiS, e.get(), c, D_S);

                        const double v = *v_it;
                        const double S = *stress_it;
                        for(size_t i = 0; i < num_den; ++i){
                            const double lKu = e->get_compliance(gradD_K[i], this->mesh->thickness, u, l);
                            const double Se = e->von_Mises_derivative(D_S, gradD_S[i], v*std::pow(S, P-2), c, u);

                            *grad_it = Sg*(Se - lKu);

                            ++x_it;
                            ++grad_it;
                        }
                        ++v_it;
                        ++stress_it;
                    }
                }
            } else {
                v_it += g->mesh.size();
                stress_it += g->mesh.size();
            }
        }
        logger::quick_log("Done.");

        logger::quick_log(result, this->c, Spn, Smax, this->alpha);
    }

    return result;
}

}

