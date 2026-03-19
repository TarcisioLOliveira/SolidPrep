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

#ifndef FINITE_ELEMENT_HPP
#define FINITE_ELEMENT_HPP

#include <vector>
#include "math/matrix.hpp"
#include "meshing.hpp"
#include "geometry.hpp"

class GlobalStiffnessMatrix;

class FiniteElement{
    public:
    static std::string get_name(){
        return "finite_element";
    }
    const double K_MIN = 1e-6;

    enum ContactType{
        RIGID,
        FRICTIONLESS_PENALTY,
        FRICTIONLESS_DISPL_LOG,
        FRICTIONLESS_DISPL_SIMPLE
    };

    struct ContactData{
        FiniteElement::ContactType contact_type = FiniteElement::ContactType::RIGID;
        double rtol_abs = 0;
        double max_step = 1;
        double step_tol = 0;
        double EPS_DISPL_SIMPLE = 0;
        double HK = 0;
        double HMU = 0;
    };

    FiniteElement(const ContactData& data);

    virtual ~FiniteElement() = default;

    void generate_matrix(const Meshing* const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda);

    void calculate_displacements(const Meshing* const mesh, std::vector<double>& load, const std::vector<double>& u0, std::vector<double>& lambda);

    void calculate_adjoint(const Meshing* const mesh, std::vector<double>& load, const std::vector<double>& u0, std::vector<double>& lambda);

    virtual std::vector<double> calculate_forces(const Meshing* const mesh, const std::vector<double>& displacements) const;

    inline const GlobalStiffnessMatrix* view_matrix() const{
        return this->matrix;
    }

    protected:
    inline void set_global_matrix(GlobalStiffnessMatrix* m){
        this->matrix = m;
    }

    const ContactData contact_data;
    GlobalStiffnessMatrix* matrix = nullptr;
    size_t u_size, l_num;
    bool first_adjoint = true;

    virtual void generate_matrix_base(const Meshing* const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const ContactType type) = 0;

    virtual void solve(std::vector<double>& load) = 0;
    virtual void reset_hessian() = 0;

    void solve_rigid(std::vector<double>& load);
    void solve_rigid_nonlinear(const Meshing* const mesh, std::vector<double>& load, const std::vector<double>& u0);
    void solve_frictionless_displ_simple(const Meshing* const mesh, std::vector<double>& load, std::vector<double>& lambda, const std::vector<double>& u0);
    void solve_frictionless_displ_log(const Meshing* const mesh, std::vector<double>& load, std::vector<double>& lambda, const std::vector<double>& u0);
    void solve_frictionless_penalty(const Meshing* const mesh, std::vector<double>& load, const std::vector<double>& u0);

    void adjoint_frictionless_displ_log(const Meshing* const mesh, std::vector<double>& load, std::vector<double>& lambda, const std::vector<double>& u_ext);

    private:

    template<typename F1, typename F2>
    double regula_falsi(F1 get_g, F2 get_rnorm_est, const double rnorm){
        double g1 = 0, g2 = 1;
        double a1 = 0, a2 = this->contact_data.max_step;
        g1 = get_g(a1);
        g2 = get_g(a2);
        logger::quick_log("g1", g1, "a1", a1);
        logger::quick_log("g2", g2, "a2", a2);
        if(!(g1 < 0 && g2 < 0) && !(std::abs(g1/g2) > 1e2 && g2 < 0.1)){
            if((g1 > 0 && g2 > 0) || std::abs(g2/g1) > 1e1 || std::abs(g1/g2) > 1e1){
                logger::quick_log("Improving starting points...");
                if(g1 > 0 && g2 > 0){
                    double rnorm_tmp = rnorm;
                    if(g1 < g2){
                        a1 = -rnorm_tmp/g1;
                        g1 = get_g(a1);
                        rnorm_tmp = get_rnorm_est(a1);
                        while(g1 < 0 && a1 < -1){
                            a1 *= 0.1;
                            g1 = get_g(a1);
                        }
                        if(g1 > 0){
                            a1 *= 10;
                            g1 = get_g(a1);
                        }
                    } else {
                        a2 = -rnorm_tmp/g2;
                        g2 = get_g(a2);
                        rnorm_tmp = get_rnorm_est(a2);
                    }
                }
                if(std::abs(g2/g1) > 1e1 && g2*g1 < 0){
                    double mtp = 0.8;
                    while(true){
                        while(g2*g1 < 0 && std::abs(g2/g1) > 1e1){
                            a2 *= mtp;
                            g2 = get_g(a2);
                        }
                        if(g2*g1 > 0){
                            a2 /= mtp;
                            g2 = get_g(a2);
                        }
                        if(std::abs(g2/g1) > 1e1){
                            mtp = (mtp - 1.0)*0.8 + 1.0;
                        } else {
                            break;
                        }
                    }
                } else if(std::abs(g1/g2) > 1e1 && g2*g1 < 0){
                    double mtp = 1.5;
                    a1 = 0.1;//-0.8*g1/(g2-g1);
                    g1 = get_g(a1);
                    while(g1 > 0){
                        a1 *= 0.1;
                        g1 = get_g(a1);
                    }
                    if(g2*g1 < 0){
                        while(true){
                            while(g2*g1 < 0 && std::abs(g1/g2) > 1e1 && a1 < 1.0){
                                a1 *= mtp;
                                g1 = get_g(a1);
                            }
                            if(g2*g1 > 0){
                                a1 /= mtp;
                                g1 = get_g(a1);
                            }
                            if(std::abs(g1/g2) > 1e1){
                                mtp = (mtp - 1.0)*0.8 + 1.0;
                            } else {
                                break;
                            }
                        }
                    }
                }
                logger::quick_log("g1", g1, "a1", a1);
                logger::quick_log("g2", g2, "a2", a2);
            }
            const double DIFF = 1e-1;
            const double stop = DIFF*std::min(std::abs(g1), std::abs(g2));
            if(g2*g1 < 0){
                logger::quick_log("Regula falsi...");
                if(g1 > g2){
                    std::swap(g1, g2);
                    std::swap(a1, a2);
                }
                // Regula falsi Illinois
                while(std::min(std::abs(g2), std::abs(g1)) > stop){
                    double diff = 0.5*g2 - g1;
                    if(std::abs(diff) < 1e-13){
                        break;
                    }
                    double c = (0.5*g2*a1 - g1*a2)/diff;
                    double gc = get_g(c);
                    if(gc > 0){
                        a2 = c;
                        g2 = gc;
                    } else {
                        a1 = c;
                        g1 = gc;
                    }
                }
                if(std::abs(g1) < std::abs(g2)){
                    a2 = a1;
                    logger::quick_log(a1, g1);
                } else {
                    logger::quick_log(a2, g2);
                }
            } else if(g2 > 0 && g1 > 0){
                a2 = 0;
            }
        }

        return a2;
    }
};

#endif
