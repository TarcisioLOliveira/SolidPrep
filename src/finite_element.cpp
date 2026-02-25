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

#include <algorithm>
#include <cblas.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include "finite_element.hpp"
#include "general_solver/mumps_general.hpp"
#include "logger.hpp"
#include "global_stiffness_matrix.hpp"
#include "math/matrix.hpp"
#include "utils.hpp"

FiniteElement::FiniteElement(const ContactData& data)
    : contact_data(data)
{

}

void FiniteElement::generate_matrix(const Meshing* const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda){
    this->u_size = u_size;
    this->l_num = l_num;
    this->generate_matrix_base(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, lambda, this->contact_data.contact_type);

}

void FiniteElement::calculate_displacements(const Meshing* const mesh, std::vector<double>& load, const std::vector<double>& u0, std::vector<double>& lambda){
    switch(this->contact_data.contact_type){
        case RIGID:
            this->solve_rigid(load);
            break;
        case FRICTIONLESS_PENALTY:
            this->solve_frictionless_penalty(mesh, load, u0);
            break;
        case FRICTIONLESS_DISPL_LOG:
            this->solve_frictionless_displ_log(mesh, load, lambda, u0);
            break;
        case FRICTIONLESS_DISPL_SIMPLE:
            this->solve_frictionless_displ_simple(mesh, load, lambda, u0);
            break;
    }
}

void FiniteElement::calculate_adjoint(const Meshing* const mesh, std::vector<double>& load, const std::vector<double>& u0, std::vector<double>& lambda){
    (void) lambda;
    switch(this->contact_data.contact_type){
        case RIGID:
            this->solve_rigid(load);
            break;
        case FRICTIONLESS_PENALTY:
            logger::log_assert(false, logger::ERROR, "adjoint equation for penalty contact not implemented");
            //this->solve_frictionless_penalty(mesh, load, u0);
            break;
        case FRICTIONLESS_DISPL_LOG:
            this->adjoint_frictionless_displ_log(mesh, load, lambda, u0);
            break;
        case FRICTIONLESS_DISPL_SIMPLE:
            logger::log_assert(false, logger::ERROR, "adjoint equation for displ_simple contact not implemented");
            //this->solve_frictionless_displ_simple(mesh, load, lambda, u0);
            break;
    }
}

void FiniteElement::solve_rigid(std::vector<double>& load){
    this->solve(load);
}

void FiniteElement::adjoint_frictionless_displ_log(const Meshing* const mesh, std::vector<double>& load, std::vector<double>& lambda, const std::vector<double>& u_ext){
    const size_t vec_size = this->u_size;

    std::vector<double> r(vec_size);
    std::copy(load.begin(), load.end(), r.begin());

    if(this->first_adjoint){
        this->reset_hessian();
        this->matrix->add_frictionless_log(mesh, mesh->node_positions[0], lambda, u_ext);
        this->first_adjoint = false;
    }
        
    this->solve(r);

    std::copy(r.begin(), r.begin() + this->u_size, load.begin());
}
void FiniteElement::solve_frictionless_displ_log(const Meshing* const mesh, std::vector<double>& load, std::vector<double>& lambda, const std::vector<double>& u0){
    const size_t vec_size = this->u_size;
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();
    const size_t l_num = mesh->lag_node_map.size();

    double E = 0;
    size_t it = 0;

    std::vector<double> f(vec_size, 0);
    std::vector<double> u(vec_size, 0);
    mesh->de_extend_vector(0, u0, u);
    std::vector<double> r(vec_size);
    std::vector<double> dr(vec_size);
    std::vector<double> u1(vec_size);

    std::vector<double> Ku(vec_size, 0);
    std::vector<double> Kd1(vec_size, 0);
    std::vector<double> u_ext(u0);

    std::copy(load.begin(), load.end(), f.begin());

    std::vector<double> lambda_source(l_num);

    double rnorm = 0;

    double step = 1.0;

    this->first_adjoint = true;

    for(auto& l:lambda){
        l *= 0.9;
    }

    this->reset_hessian();
    this->matrix->dot_vector(u, Ku);
    this->matrix->append_Ku_frictionless_log(mesh, u_ext, lambda, Ku);
    for(size_t i = 0; i < vec_size; ++i){
        r[i] = -(Ku[i] - f[i]);
    }
    this->matrix->add_frictionless_log(mesh, mesh->node_positions[0], u_ext, lambda);
    rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
    step = this->contact_data.max_step;

    do{
        // Solve linear subproblem
        logger::quick_log("r min max", *std::min_element(r.begin(), r.begin()+u_size), *std::max_element(r.begin(), r.begin()+u_size));
        logger::quick_log("u min max", *std::min_element(u.begin(), u.begin()+u_size), *std::max_element(u.begin(), u.begin()+u_size));
        this->solve(r);

        this->reset_hessian();
        std::copy(r.begin(), r.end(), dr.begin());

        double drnorm = std::sqrt(cblas_ddot(dr.size(), dr.data(), 1, dr.data(), 1));

        logger::quick_log("drnorm", drnorm);
        logger::quick_log("dr min max", *std::min_element(dr.begin(), dr.end()), *std::max_element(dr.begin(), dr.end()));


        const auto get_g = [&](double alpha)->double{
            for(size_t i = 0; i < vec_size; ++i){
                u1[i] = u[i] + alpha*dr[i];
            }
            mesh->extend_vector(0, u1, u_ext);
            std::fill(Ku.begin(), Ku.end(), 0);
            std::fill(Kd1.begin(), Kd1.end(), 0);
            this->matrix->dot_vector(u1, Ku);
            this->matrix->dot_vector(dr, Kd1);
            this->matrix->append_Ku_frictionless_log(mesh, u_ext, lambda, Ku);
            this->matrix->append_dKu_frictionless_log(mesh, u1, dr, lambda, Kd1);
            for(size_t i = 0; i < vec_size; ++i){
                r[i] = (Ku[i] - f[i]);
            }

            rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));

            const double g = cblas_ddot(vec_size, Kd1.data(), 1, r.data(), 1)/rnorm;
            return g;
        };

        const auto get_rnorm_est = [&](double alpha)->double{
            for(size_t i = 0; i < vec_size; ++i){
                u1[i] = u[i] + alpha*dr[i];
            }
            mesh->extend_vector(0, u1, u_ext);
            std::fill(Ku.begin(), Ku.end(), 0);
            this->matrix->dot_vector(u1, Ku);
            this->matrix->append_Ku_frictionless_log(mesh, u_ext, lambda, Ku);
            for(size_t i = 0; i < vec_size; ++i){
                r[i] = (Ku[i] - f[i]);
            }

            rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));

            return rnorm;
        };

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

        logger::quick_log("Applying step...");
        step = this->contact_data.max_step;
        if(a2 != this->contact_data.max_step){
            step = 1.0*a2;
        }


        for(size_t i = 0; i < vec_size; ++i){
            u[i] += step*dr[i];
        }
        mesh->extend_vector(0, u, u_ext);

        // Generate lambda update
        {
            const double MU = this->matrix->HMU;
            const double K = this->matrix->HK;
            for(const auto& e:mesh->paired_boundary){
                const auto f = e.elem->lambda_source_log(u_ext, -e.b1->normal, K);

                for(size_t i = 0; i < bnum; ++i){
                    const auto pos = mesh->lag_node_map.at(e.b1->nodes[i]->id);
                    lambda_source[pos] += MU*f[i];
                }
            }
            logger::quick_log("Min delta lambda:", *std::min_element(lambda_source.begin(), lambda_source.end()));
            logger::quick_log("Max delta lambda:", *std::max_element(lambda_source.begin(), lambda_source.end()));
            for(size_t i = 0; i < l_num; ++i){
                lambda[i] += step*lambda_source[i];
            }
            std::fill(lambda_source.begin(), lambda_source.end(), 0.0);
        }

        std::fill(Ku.begin(), Ku.end(), 0);
        this->matrix->dot_vector(u, Ku);
        E = this->matrix->constraint_frictionless_log(mesh, u_ext, lambda);
        for(size_t i = 0; i < vec_size; ++i){
            E += u[i]*(Ku[i]/2.0 - f[i]);
        }

        this->matrix->append_Ku_frictionless_log(mesh, u_ext, lambda, Ku);
        for(size_t i = 0; i < vec_size; ++i){
            r[i] = -(Ku[i] - f[i]);
        }

        rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
       
        this->matrix->add_frictionless_log(mesh, mesh->node_positions[0], u_ext, lambda);

        logger::quick_log("Iteration:", it);
        logger::quick_log("Potential energy:", E);
        logger::quick_log("||r||", rnorm);
        logger::quick_log("step:", step);
        logger::quick_log("Min lambda:", *std::min_element(lambda.begin(), lambda.end()));
        logger::quick_log("Max lambda:", *std::max_element(lambda.begin(), lambda.end()));
        logger::quick_log("");
        double max_gp = -1e100;
        double min_gp = 1e100;
        for(const auto& e:mesh->paired_boundary){
            for(size_t i = 0; i < bnum; ++i){
                const auto n1 = e.b1->nodes[i];
                const auto n2 = e.b2->nodes[i];
                const auto normal = e.b1->normal;
                double gp = 0;
                for(size_t j = 0; j < dof; ++j){
                    auto ni1 = mesh->node_positions[0][n1->u_pos[j]];
                    auto ni2 = mesh->node_positions[0][n2->u_pos[j]];
                    double u1 = 0;
                    double u2 = 0;
                    if(ni1 > -1){
                        u1 = u[ni1];
                    }
                    if(ni2 > -1){
                        u2 = u[ni2];
                    }
                    gp += (u2 - u1)*normal.Coord(1+j);
                }
                max_gp = std::max(gp, max_gp);
                min_gp = std::min(gp, min_gp);
            }
        }
        logger::quick_log("max_gp ", max_gp);
        if(min_gp < 0){
            logger::quick_log("min_gp", min_gp);
        } else {
            logger::quick_log("min_gp ", min_gp);
        }
        logger::quick_log("Gap ratio", std::abs(min_gp/max_gp)*100.0, "%");
        logger::quick_log("");

        logger::log_assert(rnorm < 1e30, logger::ERROR, "Newton process diverged, increase rtol_abs and/or decrease mult");
        ++it;
        if(step == 0){
            break;
        }
    } while((rnorm > this->contact_data.rtol_abs && std::abs(step) > this->contact_data.step_tol));// || it < 2);

    std::copy(u.begin(), u.begin() + u_size, load.begin());
}

void FiniteElement::solve_frictionless_displ_simple(const Meshing* const mesh, std::vector<double>& load, std::vector<double>& lambda, const std::vector<double>& u0){
    const size_t vec_size = this->u_size + this->l_num;
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();
    const size_t kw = mesh->elem_info->get_k_dimension();
    (void) u0;

    double E = 0;
    size_t it = 0;

    std::vector<double> f(vec_size, 0);
    std::vector<double> u(vec_size, 0);
    std::vector<double> r(vec_size);
    std::vector<double> dr(vec_size);
    std::vector<double> u1(vec_size);

    mesh->de_extend_vector(0, u0, u);
    std::copy(lambda.begin(), lambda.end(), u.begin() + this->u_size);


    std::vector<double> Ku(vec_size, 0);
    std::vector<double> Ku2(vec_size, 0);
    std::vector<double> Ku_tmp(vec_size, 0);
    std::vector<double> Kd1(vec_size, 0);
    std::vector<double> Kdd1(vec_size, 0);
    std::vector<double> Kd2(vec_size, 0);
    std::vector<double> u_ext(mesh->global_load_vector.size(), 0);
    std::vector<double> f_ext(mesh->global_load_vector.size(), 0);
    std::vector<double> Klag(l_num, 0);

    std::copy(load.begin(), load.end(), f.begin());

    mesh->extend_vector(0, load, f_ext);
    mesh->extend_vector(0, u, u_ext);

    std::copy(lambda.begin(), lambda.end(), u.begin() + u_size);
    double rnorm = 0;

    double step = this->contact_data.max_step;
    const double LAG = this->matrix->get_lag_displ_simple();

    std::vector<gp_Pnt> points(bnum);
    std::vector<long> u_pos(2*kw);
    std::vector<long> l_pos(bnum);
    std::vector<double> u_e1(kw), u_e2(kw);
    std::vector<double> l_e(bnum);
    std::vector<double> du_e1(kw), du_e2(kw);
    std::vector<double> dl_e(bnum);

    this->matrix->dot_vector(u, Ku);
    this->matrix->append_Ku_frictionless_simple(mesh, u, Ku);
    for(size_t i = 0; i < vec_size; ++i){
        r[i] = -(Ku[i] - f[i]);
    }
    this->matrix->add_frictionless_simple(mesh, mesh->node_positions[0], u_ext, lambda);

    rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
    
    logger::quick_log("Initial ||r||", rnorm);
    double u_var = 0;
    do{
        u_var = 0;

        // Solve linear subproblem
        logger::quick_log("r min max", *std::min_element(r.begin(), r.end()), *std::max_element(r.begin(), r.end()));
        logger::quick_log("r min max constr", *std::min_element(r.begin()+u_size, r.end()), *std::max_element(r.begin()+u_size, r.end()));
        logger::quick_log("u min max", *std::min_element(u.begin(), u.begin()+u_size), *std::max_element(u.begin(), u.begin()+u_size));
        logger::quick_log("l min max", *std::min_element(lambda.begin(), lambda.end()), *std::max_element(lambda.begin(), lambda.end()));
        std::copy(r.begin(), r.end(), dr.begin());
        this->solve(dr);
        logger::quick_log("dr min max constr", *std::min_element(dr.begin()+u_size, dr.end()), *std::max_element(dr.begin()+u_size, dr.end()));


        step = this->contact_data.max_step;
        
        this->reset_hessian();

        double drnorm = std::sqrt(cblas_ddot(dr.size(), dr.data(), 1, dr.data(), 1));
        logger::quick_log("drnorm", drnorm);
        logger::quick_log("dr min max", *std::min_element(dr.begin(), dr.end()), *std::max_element(dr.begin(), dr.end()));

        
        const auto get_g = [&](double alpha)->double{
            for(size_t i = 0; i < vec_size; ++i){
                u1[i] = u[i] + alpha*dr[i];
            }

            std::fill(Ku.begin(), Ku.end(), 0);
            std::fill(Kd1.begin(), Kd1.end(), 0);
            this->matrix->dot_vector(u1, Ku);
            this->matrix->dot_vector(dr, Kd1);
            this->matrix->append_Ku_frictionless_simple(mesh, u1, Ku);
            this->matrix->append_dKu_frictionless_simple(mesh, u1, dr, Kd1);
            for(size_t i = 0; i < vec_size; ++i){
                r[i] = (Ku[i] - f[i]);
            }

            rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
            
            const double g = cblas_ddot(vec_size, Kd1.data(), 1, r.data(), 1)/rnorm;
            return g;
        };
        const auto get_rnorm_est = [&](double alpha)->double{
            for(size_t i = 0; i < vec_size; ++i){
                u1[i] = u[i] + alpha*dr[i];
            }
            mesh->extend_vector(0, u1, u_ext);
            std::fill(Ku.begin(), Ku.end(), 0);
            this->matrix->dot_vector(u1, Ku);
            this->matrix->append_Ku_frictionless_simple(mesh, u_ext, Ku);
            for(size_t i = 0; i < vec_size; ++i){
                r[i] = (Ku[i] - f[i]);
            }

            rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));

            return rnorm;
        };

        double g1 = 0, g2 = 1;
        double a1 = 0, a2 = this->contact_data.max_step;
        const double search_step = 0.01;
        g1 = get_g(a1);
        g2 = get_g(a2);
        if(g1 < 0){
            for(double a = a1; a < this->contact_data.max_step + 1e-7; a += search_step){
                const double g = get_g(a);
                if(g > 0){
                    a1 = a - search_step;
                    a2 = a;
                    break;
                }
            }
        } else {
            for(double a = a1; a > -(this->contact_data.max_step + 1e-7); a -= search_step){
                const double g = get_g(a);
                if(g < 0){
                    a1 = a;
                    a2 = a + search_step;
                    break;
                }
            }
        }
        logger::quick_log("g1", g1, "a1", a1);
        logger::quick_log("g2", g2, "a2", a2);
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
                logger::quick_log("g1", g1, "a1", a1);
            } else {
                a2 = -rnorm_tmp/g2;
                g2 = get_g(a2);
                rnorm_tmp = get_rnorm_est(a2);
                logger::quick_log("g2", g2, "a2", a2);
            }
        }
        if(std::abs(g2/g1) > 1e3 && g2*g1 < 0){
            while(g2*g1 < 0){
                a2 /= 2;
                g2 = get_g(a2);
            }
            a2 *= 2;
            logger::quick_log(a2, g2);
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
        logger::quick_log("Applying step...");
        step = this->contact_data.max_step;
        if(a2 != this->contact_data.max_step){
            step = 1.0*a2;
        }


        for(size_t i = 0; i < vec_size; ++i){
            u[i] += step*dr[i];
        }

        std::fill(Ku.begin(), Ku.end(), 0);
        this->matrix->dot_vector(u, Ku);
        E = 0;
        // This is definitely wrong, but works as an estimate
        for(size_t i = 0; i < u.size(); ++i){
            E += u[i]*(Ku[i]/2 - f[i]);
        }
        this->matrix->append_Ku_frictionless_simple(mesh, u, Ku);
        for(size_t i = 0; i < vec_size; ++i){
            r[i] = -(Ku[i] - f[i]);
        }

        rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));

        std::copy(u.begin() + u_size, u.end(), lambda.begin());

        mesh->extend_vector(0, u, u_ext);
        this->matrix->add_frictionless_simple(mesh, mesh->node_positions[0], u_ext, lambda);

        logger::quick_log("Iteration:", it);
        logger::quick_log("Potential energy:", E);
        logger::quick_log("L", LAG);
        logger::quick_log("||r||", rnorm);
        logger::quick_log("du", u_var);
        const double Kd1norm = std::sqrt(cblas_ddot(Kd1.size(), Kd1.data(), 1, Kd1.data(), 1));
        logger::quick_log("||Kd1||:", Kd1norm);
        logger::quick_log("Step:", step);
        logger::quick_log("");
        double max_gp = -1e100;
        double min_gp = 1e100;
        for(const auto& e:mesh->paired_boundary){
            for(size_t i = 0; i < bnum; ++i){
                const auto n1 = e.b1->nodes[i];
                const auto n2 = e.b2->nodes[i];
                const auto normal = e.b1->normal;
                double gp = 0;
                for(size_t j = 0; j < dof; ++j){
                    auto ni1 = mesh->node_positions[0][n1->u_pos[j]];
                    auto ni2 = mesh->node_positions[0][n2->u_pos[j]];
                    double u1 = 0;
                    double u2 = 0;
                    if(ni1 > -1){
                        u1 = u[ni1];
                    }
                    if(ni2 > -1){
                        u2 = u[ni2];
                    }
                    gp += (u2 - u1)*normal.Coord(1+j);
                }
                max_gp = std::max(gp, max_gp);
                min_gp = std::min(gp, min_gp);
            }
        }
        logger::quick_log("max_gp ", max_gp);
        if(min_gp < 0){
            logger::quick_log("min_gp", min_gp);
        } else {
            logger::quick_log("min_gp ", min_gp);
        }

        logger::log_assert(rnorm < 1e30, logger::ERROR, "Newton process diverged, increase rtol_abs and/or decrease mult");
        ++it;
        if(step == 0){
            break;
        }
    } while((rnorm > this->contact_data.rtol_abs && std::abs(step) > this->contact_data.step_tol) || it < 1);

    std::copy(u.begin(), u.begin() + u_size, load.begin());
}
void FiniteElement::solve_frictionless_penalty(const Meshing* const mesh, std::vector<double>& load, const std::vector<double>& u0){
    const size_t vec_size = this->u_size;
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();

    double E = 0;
    size_t it = 0;

    std::vector<double> f(vec_size, 0);
    std::vector<double> u(vec_size, 0);
    std::vector<double> u1(vec_size, 0);
    std::vector<double> uold1(vec_size, 0);
    std::vector<double> uold2(vec_size, 0);
    std::vector<double> r(vec_size);
    std::vector<double> dr(vec_size);
    std::vector<double> lambda(mesh->lag_node_map.size(), 0);

    mesh->de_extend_vector(0, u0, u);

    // Initiate "exterior" point
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < bnum; ++i){
            const auto n2 = e.b2->nodes[i];
            const auto normal = e.b1->normal;
            for(size_t j = 0; j < dof; ++j){
                auto ni2 = mesh->node_positions[0][n2->u_pos[j]];
                if(u[ni2] == 0){
                    if(std::abs(normal.Coord(1+j)) > 1e-10){
                        const double sign = (normal.Coord(1+j) > 0) ? 1 : -1;
                        u[ni2] = -sign*1e-6;
                    }
                }
            }

            // If I ever start with non zero `u`, edit this
            /*
            for(size_t j = 0; j < dof; ++j){
                auto ni1 = mesh->node_positions[0][n1->u_pos[j]];
                auto ni2 = mesh->node_positions[0][n2->u_pos[j]];
                if(ni1 > -1){
                    du1 = dr[ni1]*normal.Coord(1+j);
                }
                if(ni2 > -1){
                    du2 = dr[ni2]*normal.Coord(1+j);
                }
            }
            gp = (du2 - du1);
            if(gp < -1e-10){
                if(std::abs(du2) < std::abs(du1)){
                    double a = du2/du1;
                    for(size_t j = 0; j < dof; ++j){
                        auto ni1 = mesh->node_positions[0][n1->u_pos[j]];
                        if(ni1 > -1 && std::abs(normal.Coord(1+j)) > 1e-10){
                            dr[ni1] *= a;
                        }
                    }
                } else {
                    double a = du1/du2;
                    for(size_t j = 0; j < dof; ++j){
                        auto ni2 = mesh->node_positions[0][n2->u_pos[j]];
                        if(ni2 > -1 && std::abs(normal.Coord(1+j)) > 1e-10){
                            dr[ni2] *= a;
                        }
                    }
                }
            }
            */
        }
    }

    const auto apply_multiplier = [&](const std::vector<double>& u, std::vector<double>& l){
        for(const auto& e:mesh->paired_boundary){
            for(size_t i = 0; i < bnum; ++i){
                const auto n1 = e.b1->nodes[i];
                const auto n2 = e.b2->nodes[i];
                const auto normal = -e.b1->normal;
                double gp = 0;
                for(size_t j = 0; j < dof; ++j){
                    auto ni1 = mesh->node_positions[0][n1->u_pos[j]];
                    auto ni2 = mesh->node_positions[0][n2->u_pos[j]];
                    double u1 = 0;
                    double u2 = 0;
                    if(ni1 > -1){
                        u1 = u[ni1];
                    }
                    if(ni2 > -1){
                        u2 = u[ni2];
                    }
                    gp += (u2 - u1)*normal.Coord(1+j);
                }
                if(gp > 0){
                    l[mesh->lag_node_map.at(n1->id)] += 1e-3*gp*this->contact_data.EPS_DISPL_SIMPLE;
                }
            }
        }

    };

    std::vector<double> Ku(vec_size, 0);
    std::vector<double> Kd(vec_size, 0);
    std::vector<double> u_ext(mesh->global_load_vector.size(), 0);
    std::vector<double> f_ext(mesh->global_load_vector.size(), 0);
    std::vector<double> Klag(l_num, 0);

    std::vector<bool> is_contact_node(vec_size, false);

    std::copy(load.begin(), load.end(), f.begin());

    mesh->extend_vector(0, load, f_ext);
    mesh->extend_vector(0, u, u_ext);
    double rnorm = 0;

    double step = 1.0;
    std::fill(Ku.begin(), Ku.end(), 0);
    this->matrix->add_contacts(mesh, mesh->node_positions[0], u_ext);
    this->matrix->dot_vector(u, Ku);
    for(size_t i = 0; i < u_size; ++i){
        r[i] = -(Ku[i] - f[i]);
    }
    rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
    logger::quick_log("Initial ||r||:", rnorm);

    step = this->contact_data.max_step;
    
    do{

        // Solve linear subproblem
        std::copy(r.begin(), r.end(), dr.begin());
        this->solve(dr);

        this->reset_hessian();


        // Damped Newton-Rhapson
        /*
            const auto get_g = [&](double alpha)->double{
                for(size_t i = 0; i < u.size(); ++i){
                    u1[i] = u[i] + alpha*dr[i];
                }

                std::fill(Ku.begin(), Ku.end(), 0);
                std::fill(Kd.begin(), Kd.end(), 0);
                mesh->extend_vector(0, u1, u_ext);
                this->matrix->dot_vector(u1, Ku);
                this->matrix->dot_vector(dr, Kd);
                this->matrix->append_Ku_penalty(mesh, mesh->node_positions[0], u_ext, lambda, Ku);
                this->matrix->append_dKu_penalty(mesh, mesh->node_positions[0], u_ext, lambda, dr, Kd);
                for(size_t i = 0; i < vec_size; ++i){
                    r[i] = (Ku[i] - f[i]);
                }

                rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
                return cblas_ddot(vec_size, Kd.data(), 1, r.data(), 1)/rnorm;
            };

            //get_LAG();
            double g1 = 0, g2 = 1;
            double a1 = 0, a2 = this->max_step;
            //if(drnorm > 1e3){
            //    a2 = this->max_step/std::pow(drnorm, 3.0/2.0);
            //}
            g1 = get_g(a1);
            double g0 = g1;
            g2 = get_g(a2);
            logger::quick_log("g0", g0);
            logger::quick_log("g2", g2);
            if(g0 > 0){
                logger::quick_log("Improving starting points...");
                a1 = -1e-20;
                g1 = get_g(a1);
                g0 = g1;
                while(g0 > 0){
                    a1 *= 2;
                    g1 = get_g(a1);
                    g0 = g1;
                    if(g0 < 0) break;
                    g1 = get_g(std::abs(a1));
                    g0 = g1;
                    if(g0 < 0){
                        a1 = -a1;
                    }
                    if(std::abs(a1) >= 1){
                        break;
                    }
                }
            }
            if(std::abs(g2/g1) > 1e4 && g2*g1 < 0){
                a2 = std::abs(g1/g2);
                g2 = get_g(a2);
                while(g2*g1 > 0){
                    a2 *= 2;
                    g2 = get_g(a2);
                }
                logger::quick_log(a2, g2);
            }
            logger::quick_log("g1", g1);
            logger::quick_log("g2", g2);
            if(a1 > -1 && g0 < 0 && a2 > 1e-5){
                const double DIFF = 1e-2;
                const double stop = DIFF*std::min(std::abs(g0), std::abs(g2));
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
                logger::quick_log("Applying step...");
                step = this->max_step;
                if(a2 != this->max_step){
                    step = 0.9*a2;
                    //step = 1.0*a2;
                }
            } else {
                step = 0;
            }
        */
            for(size_t i = 0; i < u.size(); ++i){
                u[i] += step*dr[i];
            }
            apply_multiplier(u, lambda);
            

        mesh->extend_vector(0, u, u_ext);
        std::fill(Ku.begin(), Ku.end(), 0);
        this->matrix->dot_vector(u, Ku);
        E = 0;
        for(size_t i = 0; i < u_size; ++i){
            E += u[i]*(Ku[i]/2 - f[i]);
        }
        this->matrix->append_Ku_penalty(mesh, mesh->node_positions[0], u_ext, lambda, Ku);
        for(size_t i = 0; i < u_size; ++i){
            r[i] = -(Ku[i] - f[i]);
        }

        this->matrix->add_contacts(mesh, mesh->node_positions[0], u_ext);

        logger::quick_log("Iteration:", it);
        logger::quick_log("Potential energy:", E);
        rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
        logger::quick_log("||r||:", rnorm);
        logger::quick_log("Step:", step);
        logger::quick_log("");
        double max_gp = -1e100;
        double min_gp = 1e100;
        for(const auto& e:mesh->paired_boundary){
            for(size_t i = 0; i < bnum; ++i){
                const auto n1 = e.b1->nodes[i];
                const auto n2 = e.b2->nodes[i];
                const auto normal = e.b1->normal;
                double gp = 0;
                for(size_t j = 0; j < dof; ++j){
                    auto ni1 = mesh->node_positions[0][n1->u_pos[j]];
                    auto ni2 = mesh->node_positions[0][n2->u_pos[j]];
                    double u1 = 0;
                    double u2 = 0;
                    if(ni1 > -1){
                        u1 = u[ni1];
                    }
                    if(ni2 > -1){
                        u2 = u[ni2];
                    }
                    gp += (u2 - u1)*normal.Coord(1+j);
                }
                max_gp = std::max(gp, max_gp);
                min_gp = std::min(gp, min_gp);
            }
        }
        logger::quick_log("max_gp ", max_gp);
        logger::quick_log("min_gp", min_gp);

        logger::log_assert(rnorm < 1e30, logger::ERROR, "Newton process diverged, increase rtol_abs and/or decrease mult");
        ++it;
    } while((rnorm > this->contact_data.rtol_abs || it < 1) && step > 0 && it < 50);

    std::copy(u.begin(), u.begin() + u_size, load.begin());
}

std::vector<double> FiniteElement::calculate_forces(const Meshing* const mesh, const std::vector<double>& displacements) const{
    logger::quick_log("Calculating forces...");
    const size_t dof       = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    std::vector<double> results(mesh->node_list.size()*dof, 0);
  
   for(auto& g:mesh->geometries){
       const size_t num_mat = g->number_of_materials();
       if(num_mat == 1){
            for(auto& e:g->mesh){
                const gp_Pnt c = e->get_centroid();
                const auto D = g->materials.get_D(e.get(), c);
                auto f = e->get_internal_loads(D, mesh->thickness, displacements);
                for(size_t n = 0; n < node_num; ++n){
                    auto& node = e->nodes[n];
                    for(size_t i = 0; i < dof; ++i){
                        results[node->id*dof + i] += f[n*dof+i];
                    }
                }
            }
       } else {
           logger::quick_log(num_mat == 1, logger::ERROR, "calculate_forces() not implemented for geometries with more than one material");
       }
    }
    logger::quick_log("Done.");

    return results;
}
