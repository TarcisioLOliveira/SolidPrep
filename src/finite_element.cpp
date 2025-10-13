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
#include <set>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <unordered_map>
#include "finite_element.hpp"
#include "logger.hpp"
#include "global_stiffness_matrix.hpp"
#include "math/matrix.hpp"
#include "project_data.hpp"
#include "utils.hpp"
#include "optimization/MMASolver.hpp"

FiniteElement::FiniteElement(const ContactData& data)
    :contact_type(data.contact_type),
     rtol_abs(data.rtol_abs),
     step_tol(data.step_tol),
     start_lag_simple(data.EPS_DISPL_SIMPLE),
     max_step(data.max_step)
{

}

void FiniteElement::generate_matrix(const Meshing* const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda){
    this->u_size = u_size;
    this->l_num = l_num;
    this->generate_matrix_base(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, lambda, this->contact_type);

    if(this->contact_type == FRICTIONLESS_DISPL_SIMPLE){
        this->matrix->set_lag_displ_simple(this->start_lag_simple);
    } else if(this->contact_type == FRICTIONLESS_DISPL_CONSTR){
        const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
        const size_t num_bound_nodes = mesh->elem_info->get_boundary_nodes_per_element();
        const size_t dof = mesh->elem_info->get_dof_per_node();
        const double t = mesh->thickness;

        std::set<Node*> ref_elem_nodes;
        std::set<Node*> ref_bound_nodes_e1;
        std::map<Node*, Node*> ref_bound_nodes_e2;
        for(auto& b:mesh->paired_boundary){
            auto& e = b.elem->e2;
            auto& be = b.elem;
            for(size_t i = 0; i < num_nodes; ++i){
                ref_elem_nodes.insert(e->nodes[i]);
            }
            for(size_t i = 0; i < num_bound_nodes; ++i){
                for(size_t j = 0; j < num_nodes; ++j){
                    ref_bound_nodes_e1.insert(be->nodes[i]);
                    if(be->nodes[i]->point.IsEqual(e->nodes[j]->point, Precision::Confusion())){
                        ref_bound_nodes_e2[be->nodes[i]] = e->nodes[j];
                    }
                }
            }
        }
        for(auto n:ref_elem_nodes){
            this->constr_id_map[n] = -1;
        }
        long id = 0;
        for(auto n:ref_bound_nodes_e1){
            this->constr_id_map[n] = id;
            this->constr_id_map[ref_bound_nodes_e2[n]] = id;
            ++id;
        }
        ref_elem_nodes.clear();
        ref_bound_nodes_e1.clear();
        ref_bound_nodes_e2.clear();

        const size_t matrix_size = dof*id;
        this->constr_solver.initialize_matrix(true, matrix_size);
        this->constr_force.resize(matrix_size);

        const math::Matrix K(math::Matrix::identity(3));
        std::vector<gp_Pnt> bound(num_bound_nodes);
        std::vector<long> pos(num_nodes*dof);
        math::Matrix M;
        for(auto& b:mesh->paired_boundary){
            auto& c = b.elem;
            auto& e = b.elem->e2;
            for(size_t i = 0; i < num_bound_nodes; ++i){
                bound[i] = c->nodes[i]->point;
            }
            for(size_t i = 0; i < num_nodes; ++i){
                const long id = this->constr_id_map[e->nodes[i]];
                if(id >= 0){
                    for(size_t j = 0; j < dof; ++j){
                        pos[i*dof + j] = id*dof + j;
                    }
                } else {
                    for(size_t j = 0; j < dof; ++j){
                        pos[i*dof + j] = -1;
                    }
                }
            }
            M = e->get_R(K, t, bound);

            this->constr_solver.add_element(M, pos);
        }

        /*
        const size_t matrix_size = id;
        this->constr_solver.initialize_matrix(true, matrix_size);
        this->constr_force.resize(matrix_size);

        std::vector<gp_Pnt> bound(num_bound_nodes);
        std::vector<long> pos(num_nodes);
        math::Matrix M;
        for(auto& b:mesh->paired_boundary){
            auto& c = b.elem;
            auto& e = b.elem->e2;
            for(size_t i = 0; i < num_bound_nodes; ++i){
                bound[i] = c->nodes[i]->point;
            }
            for(size_t i = 0; i < num_nodes; ++i){
                const long id = this->constr_id_map[e->nodes[i]];
                pos[i] = id;
            }
            M = e->robin_1dof(t, bound);

            this->constr_solver.add_element(M, pos);
        }
        */

        this->constr_solver.compute();

        std::set<Node*> unique_nodes;
        for(const auto& e:mesh->paired_boundary){
            for(size_t i = 0; i < num_nodes; ++i){
                //const auto n1 = e.elem->e1->nodes[i];
                const auto n2 = e.elem->e2->nodes[i];
                const long id = this->constr_id_map.at(n2);
                if(id >= 0){
                    unique_nodes.insert(n2);
                }
            }
        }
        this->boundary_reference_nodes.insert(this->boundary_reference_nodes.begin(), unique_nodes.begin(), unique_nodes.end());
    }    
}

void FiniteElement::calculate_displacements(const Meshing* const mesh, std::vector<double>& load, const std::vector<double>& u0, std::vector<double>& lambda, const bool topopt, const std::vector<math::Matrix>& D_cache){
    switch(this->contact_type){
        case RIGID:
            this->solve_rigid(load);
            break;
        case FRICTIONLESS_PENALTY:
            this->solve_frictionless_penalty(mesh, load, u0);
            break;
        case FRICTIONLESS_DISPL_LOG:
            this->solve_frictionless_displ_log(mesh, load, u0);
            break;
        case FRICTIONLESS_DISPL_SIMPLE:
            this->solve_frictionless_displ_simple(mesh, load, lambda, u0);
            break;
        case FRICTIONLESS_DISPL_CONSTR:
            this->solve_frictionless_displ(mesh, load, lambda, topopt, D_cache, u0);
            break;
    }
}

void FiniteElement::solve_rigid(std::vector<double>& load){
    this->solve(load);
}
void FiniteElement::solve_frictionless_displ_log(const Meshing* const mesh, std::vector<double>& load, const std::vector<double>& u0){
    const size_t vec_size = this->u_size + 1;
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t kw = mesh->elem_info->get_k_dimension();
    const size_t dof = mesh->elem_info->get_dof_per_node();
    //const size_t node_num = mesh->elem_info->get_nodes_per_element();
    //auto& node_positions = mesh->node_positions[0];
    (void) u0;

    double E = 0;
    size_t it = 0;

    std::vector<double> f(vec_size, 0);
    std::vector<double> u(vec_size, 0);
    std::vector<double> r(vec_size);
    std::vector<double> G(vec_size);
    std::vector<double> HG(vec_size);
    std::vector<double> r2(vec_size);
    std::vector<double> dr(vec_size);
    std::vector<double> u_test(vec_size);
    std::vector<double> u1(vec_size);
    std::vector<double> u2(vec_size);
    std::vector<double> uold1(vec_size);
    std::vector<double> uold2(vec_size);

    //std::copy(u0.begin(), u0.end(), u.begin());

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

    double& LAG = this->matrix->LAG_DISPL_LOG;
    LAG = 1e6;
    u.back() = LAG;
    double rnorm = 0;
    //double old_rnorm = 0;

    double step = 1.0;
    //double gamma = 0;

    std::vector<gp_Pnt> points(bnum);
    std::vector<long> u_pos(2*kw);
    std::vector<long> l_pos(bnum);
    std::vector<double> u_e1(kw), u_e2(kw);
    std::vector<double> l_e(bnum);
    std::vector<double> du_e1(kw), du_e2(kw);
    std::vector<double> dl_e(bnum);

    // // Initiate interior point
    // for(const auto& e:mesh->paired_boundary){
    //     for(size_t i = 0; i < bnum; ++i){
    //         const auto n2 = e.b2->nodes[i];
    //         const auto normal = e.b1->normal;
    //         for(size_t j = 0; j < dof; ++j){
    //             auto ni2 = mesh->node_positions[0][n2->u_pos[j]];
    //             if(std::abs(normal.Coord(1+j)) > 1e-10){
    //                 const double sign = (normal.Coord(1+j) > 0) ? 1 : -1;
    //                 u[ni2] = -sign*(1e-5);
    //             }
    //         }

    //         // If I ever start with non zero `u`, edit this
    //         /*
    //         for(size_t j = 0; j < dof; ++j){
    //             auto ni1 = mesh->node_positions[0][n1->u_pos[j]];
    //             auto ni2 = mesh->node_positions[0][n2->u_pos[j]];
    //             if(ni1 > -1){
    //                 du1 = dr[ni1]*normal.Coord(1+j);
    //             }
    //             if(ni2 > -1){
    //                 du2 = dr[ni2]*normal.Coord(1+j);
    //             }
    //         }
    //         gp = (du2 - du1);
    //         if(gp < -1e-10){
    //             if(std::abs(du2) < std::abs(du1)){
    //                 double a = du2/du1;
    //                 for(size_t j = 0; j < dof; ++j){
    //                     auto ni1 = mesh->node_positions[0][n1->u_pos[j]];
    //                     if(ni1 > -1 && std::abs(normal.Coord(1+j)) > 1e-10){
    //                         dr[ni1] *= a;
    //                     }
    //                 }
    //             } else {
    //                 double a = du1/du2;
    //                 for(size_t j = 0; j < dof; ++j){
    //                     auto ni2 = mesh->node_positions[0][n2->u_pos[j]];
    //                     if(ni2 > -1 && std::abs(normal.Coord(1+j)) > 1e-10){
    //                         dr[ni2] *= a;
    //                     }
    //                 }
    //             }
    //         }
    //         */
    //     }
    // }

    //const auto norm = [](const std::vector<double>& r)->double{
    //    return std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
    //};

    //this->matrix->dot_vector(u, Ku);
    //std::fill(Ku.begin() + u_size, Ku.end(), 0);
    this->matrix->append_Ku_frictionless_log(mesh, u_ext, Ku);
    for(size_t i = 0; i < vec_size; ++i){
        r[i] = -(Ku[i] - f[i]);
    }
    this->matrix->add_frictionless_log(mesh, mesh->node_positions[0], u_ext);
    const double u_max =  0.1;
    const double u_min = -0.1;
    const double mma_step = 0.001;
    std::vector<double> S(vec_size, mma_step);
    S.back() = 1;
    step = 1;
    rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
    step = this->max_step;

    double u_var = 0;
    do{
        u_var = 0;
        // Solve linear subproblem
        logger::quick_log("r min max", *std::min_element(r.begin(), r.begin()+u_size), *std::max_element(r.begin(), r.begin()+u_size));
        logger::quick_log("u min max", *std::min_element(u.begin(), u.begin()+u_size), *std::max_element(u.begin(), u.begin()+u_size));
        this->solve(r);

        this->reset_hessian();
        std::copy(r.begin(), r.end(), dr.begin());

        double drnorm = std::sqrt(cblas_ddot(dr.size(), dr.data(), 1, dr.data(), 1));

        logger::quick_log("drnorm", drnorm);
        logger::quick_log("dr min max", *std::min_element(dr.begin(), dr.end()), *std::max_element(dr.begin(), dr.end()));

        const auto sign = [](const double x)->double{
            if(x > 0){
                return 1;
            } else if(x < 0){
                return -1;
            }
            return 0;
        };

        const auto update_vars = [&](const size_t i, const double min, const double max){
            //const double s = std::min(std::abs(dr[i]), 1.0)*S[i];
            const double s = S[i];
            const double U = std::min(u[i] + s, max);
            const double UX = U - u[i];
            const double L = std::max(u[i] - s, min);
            const double LX = u[i] - L;
            const double P =  (UX*UX*UX)/(UX + LX) + 1e-11/(max - min);
            const double Q = -(LX*LX*LX)/(UX + LX) - 1e-11/(max - min);
            const double V = dr[i];
            //const double V = -std::abs(dr[i])*sign(r[i]);
            const bool P_nonzero = std::abs(P) > 1e-15;
            const bool Q_nonzero = std::abs(Q) > 1e-15;
            logger::log_assert(P_nonzero && Q_nonzero, logger::ERROR, "P or Q became zero somehow");
            double R = 0;
            if(UX > 0){
                R -= P/UX;
            }
            if(LX > 0){
                R -= Q/LX;  
            }

            const double a = (V - R);
            const double b = -((V - R)*(U + L) - (P - Q));
            const double c = (V - R)*U*L + Q*U - P*L;

            const double delta = b*b - 4*a*c;
            if(delta >= 0 && std::abs(a) > 1e-7){
                const double xt1 = -2*c/(b + std::sqrt(delta));
                const double xt2 = -2*c/(b - std::sqrt(delta));
                //const bool xt1_within = (L < xt1 || std::abs(xt1 - L) < 1e-5) && (xt1 < U || std::abs(U - xt1) < 1e-5);
                //const bool xt2_within = (L < xt2 || std::abs(xt2 - L) < 1e-5) && (xt2 < U || std::abs(U - xt2) < 1e-5);
                if((V-R) > 0){
                    u[i] = std::max(xt1, xt2);
                    u[i] = std::min(u[i], U);
                } else {
                    u[i] = std::min(xt1, xt2);
                    u[i] = std::max(u[i], L);
                }
                //if(std::abs(V) > 0.001){
                //if(!xt1_within && !xt2_within){
                //    if(std::abs(xt1 - L) < 1e-10 || std::abs(xt2 - L) < 1e-10){
                //        u[i] = L;
                //    } else if(std::abs(xt1 - U) < 1e-10 || std::abs(xt2 - U) < 1e-10) {
                //        u[i] = U;
                //    }
                //} else {
                //    logger::log_assert(xt1_within != xt2_within || std::abs(xt1 - xt2) < 1e-5, logger::ERROR, "incorrect assumption about test points placement, L: {} xt1: {} xt2: {} U: {} xi: {} P: {} Q: {} R: {}, V: {}", L, xt1, xt2, U, u[i], P, Q, R, V);
                //    if(xt1_within){
                //        u[i] = xt1;
                //    } else {
                //        u[i] = xt2;
                //    }
                //}
            } else if(delta >= 0){
                const double xt1 = -2*c/(b + std::sqrt(delta));
                u[i] = xt1;
            }
            //logger::log_assert(L - 1e-15 < u[i] && u[i] < U + 1e-15, logger::ERROR, "Update violated asymptotes, L: {} U: {} xi: {} P: {} Q: {} R: {} V: {} P_nonzero: {} Q_nonzero: {}", L, U, u[i], P, Q, R, V, P_nonzero, Q_nonzero);
            //logger::log_assert(min - 1e-15 < u[i] && u[i] < max + 1e-15, logger::ERROR, "Update violated boundaries, L: {} U: {} xi: {} P: {} Q: {} R: {} V: {} P_nonzero: {} Q_nonzero: {}", L, U, u[i], P, Q, R, V, P_nonzero, Q_nonzero);
            // Otherwise, x == x0
            if(it > 1){
                if((uold2[i] - uold1[i])*(uold1[i] - u[i]) < 0){
                    S[i] *= 0.9;
                    S[i] = std::max(S[i], 1e-14);
                    //S[i] = std::abs(uold1[i] - u[i]);
                //} else {
                //    S[i] = std::min(1.05*S[i], mma_step);
                }
            }
            uold2[i] = uold1[i];
            uold1[i] = u[i];
        };
        //u.back() = LAG;
        //#pragma omp parallel
        //{
        //    #pragma omp for
        //    for(size_t i = 0; i < u_size; ++i){
        //        update_vars(i, u_min, u_max);
        //    }
        //    update_vars(u_size, -100, 100);
        //}
        //LAG = u.back();

            /*
            const auto get_g = [&](double alpha)->double{
                for(size_t i = 0; i < u.size(); ++i){
                    u1[i] = u[i] + alpha*dr[i];
                }
                LAG = u1.back();
                mesh->extend_vector(0, u1, u_ext);
                std::fill(Ku.begin(), Ku.end(), 0);
                std::fill(Kd1.begin(), Kd1.end(), 0);
                this->matrix->dot_vector(u1, Ku);
                this->matrix->dot_vector(dr, Kd1);
                this->matrix->append_Ku_frictionless_log(mesh, u_ext, Ku);
                this->matrix->append_dKu_frictionless_log(mesh, u_ext, u1, dr, alpha, Kd1);
                for(size_t i = 0; i < vec_size; ++i){
                    r[i] = (Ku[i] - f[i]);
                }

                rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));

                return cblas_ddot(vec_size, Kd1.data(), 1, r.data(), 1)/rnorm;
            };

            double g1 = 0, g2 = 1;
            double a1 = 0, a2 = this->max_step;
            g1 = get_g(a1);
            double g0 = g1;
            g2 = get_g(a2);
            logger::quick_log("g1", g0, "a1", a1);
            logger::quick_log("g2", g2, "a2", a2);
            logger::quick_log("Improving starting points...");
            if(g0 > 0){
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
            if(a1 > -1){
                if(std::abs(g2/g1) > 1e4 && g2*g1 < 0){
                    a2 = 0.1*std::abs(g1/g2);
                    g2 = get_g(a2);
                    while(g2*g1 > 0){
                        a2 *= 2;
                        g2 = get_g(a2);
                    }
                    logger::quick_log(a2, g2);
                }
                const double DIFF = 1e-4;
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
                    //step = 0.8*a2;
                    step = 1.0*a2;
                }
            } else {
                step = 0;
            }
            */

            for(size_t i = 0; i < vec_size; ++i){
                u[i] += step*dr[i];
            }
            LAG = u.back();

            std::fill(Ku.begin(), Ku.end(), 0);
            this->matrix->dot_vector(u, Ku);
            mesh->extend_vector(0, u, u_ext);
            this->matrix->append_Ku_frictionless_log(mesh, u_ext, Ku);
            for(size_t i = 0; i < vec_size; ++i){
                r[i] = -(Ku[i] - f[i]);
            }
            logger::quick_log("Constraint:", -r.back());
            logger::quick_log("Multiplier:", LAG);


        for(size_t i = 0; i < vec_size; ++i){
            u_var = std::max(std::abs(u[i] - uold2[i]), u_var);
        }

        const double rnorm2 = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
        //if(rnorm2 > rnorm){
        //    step *= 0.7;
        //}
        rnorm = rnorm2;
       
        this->matrix->add_frictionless_log(mesh, mesh->node_positions[0], u_ext);
        E = 0;
        // This is definitely wrong, but works as an estimate
        for(size_t i = 0; i < vec_size; ++i){
            E += u[i]*(Ku[i]/2 - f[i]);
        }

        logger::quick_log("Iteration:", it);
        logger::quick_log("Potential energy:", E);
        logger::quick_log("||r||", rnorm);
        const double Kd1norm = std::sqrt(cblas_ddot(Kd1.size(), Kd1.data(), 1, Kd1.data(), 1));
        logger::quick_log("||Kd1||:", Kd1norm);
        logger::quick_log("step:", step);
        logger::quick_log("du:", u_var);
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
        //if(b > 1e-6 || step <= 0){
        //    break;
        //}
        if(step == 0){
            break;
        }
    } while((rnorm > this->rtol_abs && std::abs(step) > this->step_tol) || it < 2);

    std::copy(u.begin(), u.begin() + u_size, load.begin());
}
void FiniteElement::solve_frictionless_displ(const Meshing* const mesh, std::vector<double>& load, std::vector<double>& lambda, const bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u0){
    const size_t vec_size = this->u_size + this->l_num;
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();

    double E = 0;
    size_t it = 0;

    std::vector<double> f(vec_size, 0);
    std::vector<double> u(vec_size, 0);
    std::vector<double> r(vec_size);
    std::vector<double> r2(vec_size);
    std::vector<double> dr(vec_size);
    std::vector<double> u_test(vec_size);

    mesh->de_extend_vector(0, u0, u);

    //std::fill(u.begin(), u.begin() + u_size, 0.0);
    //std::copy(lambda.begin(), lambda.begin() + 2*l_num, u.begin() + u_size);
    std::vector<double> Ku(vec_size, 0);
    std::vector<double> Kd(vec_size, 0);
    std::vector<double> u_ext(mesh->global_load_vector.size(), 0);
    std::vector<double> f_ext(mesh->global_load_vector.size(), 0);
    std::vector<double> Klag(l_num, 0);

    std::copy(load.begin(), load.end(), f.begin());

    mesh->extend_vector(0, load, f_ext);

    std::copy(lambda.begin(), lambda.end(), u.begin() + u_size);

    double rnorm = 0;
    //double old_rnorm = 0;

    //std::fill(Ku.begin() + u_size + 2*l_num, Ku.end(), 0);
    //double old_step = 1;
    double step = 1.0;
    //std::fill(Ku.begin(), Ku.end(), 0);
    //this->matrix->dot_vector(u, Ku);
    //std::fill(Ku.begin() + u_size, Ku.end(), 0);
    mesh->extend_vector(0, u, u_ext);

    this->apply_lambda(mesh, lambda, u_ext);
    this->matrix->append_Ku_frictionless(mesh, u, Ku, D_cache, topopt);
    for(size_t i = 0; i < vec_size; ++i){
        r[i] = -(Ku[i] - f[i]);
    }
    this->matrix->add_frictionless_part2(mesh, mesh->node_positions[0], u_ext, lambda, D_cache, topopt);
    rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));

    //++it;
    do{
        //this->reset_hessian();
        //std::copy(u.begin() + u_size, u.begin() + u_size + 2*l_num, lambda.begin());


        // Solve linear subproblem
        logger::quick_log("r min max", *std::min_element(r.begin(), r.end()), *std::max_element(r.begin(), r.end()));
        logger::quick_log("u min max", *std::min_element(u.begin(), u.end()), *std::max_element(u.begin(), u.end()));
        this->solve(r);

        step = this->max_step;
        
        this->reset_hessian();
        std::copy(r.begin(), r.end(), dr.begin());

        double drnorm = std::sqrt(cblas_ddot(dr.size(), dr.data(), 1, dr.data(), 1));
        logger::quick_log("drnorm", drnorm);
        logger::quick_log("dr min max", *std::min_element(dr.begin(), dr.end()), *std::max_element(dr.begin(), dr.end()));


        // Damped Newton-Rhapson
            const auto get_g = [&](double alpha)->double{
                for(size_t i = 0; i < u.size(); ++i){
                    u_test[i] = u[i] + alpha*dr[i];
                }
                std::fill(Ku.begin(), Ku.end(), 0);
                std::fill(Kd.begin(), Kd.end(), 0);
                this->matrix->dot_vector(u_test, Ku);
                this->matrix->dot_vector(dr, Kd);
                this->matrix->append_Ku_frictionless(mesh, u_test, Ku, D_cache, topopt);
                this->matrix->append_dKu_frictionless(mesh, u_test, dr, alpha, Kd, D_cache, topopt);
                for(size_t i = 0; i < vec_size; ++i){
                    r[i] = (Ku[i] - f[i]);
                }
                return cblas_ddot(vec_size, Kd.data(), 1, r.data(), 1);
            };
            double g1 = 0, g2 = 1;
            double a1 = 0, a2 = this->max_step;
            g1 = get_g(a1);
            double g0 = g1;
            g2 = get_g(a2);
            const double DIFF = 1e-2;
            logger::quick_log("g0", g0);
            logger::quick_log("g2", g2);
            double stop_criterium = 0;
            stop_criterium = DIFF*std::abs(g0);
            if(g2*g1 < 0 && std::abs(g2) > stop_criterium){
                while(std::abs(g2) > stop_criterium){
                    double a2_old = a2;
                    double diff = g2 - g1;
                    if(std::abs(diff) < 1e-13){
                        break;
                    }
                    a2 -= g2*(a2 - a1)/diff;;
                    a1 = a2_old;
                    g1 = g2;
                    g2 = get_g(a2);
                }
            } else if(g1 > 0 && g2 > 0){
                a2 = 0;
                g2 = get_g(a2);
            } else {
                a2 = 1;
                g2 = get_g(a2);
            }
            logger::quick_log("g2_final", g2);
            step = a2;
            if(step < 0){
                step = 1e-5;
                g2 = get_g(step);
            }
            for(size_t i = 0; i < vec_size; ++i){
                r[i] *= -1.0;
            }

        std::copy(u_test.begin(), u_test.end(), u.begin());
        std::copy(u.begin() + u_size, u.end(), lambda.begin());

        mesh->extend_vector(0, u, u_ext);
        this->apply_lambda(mesh, lambda, u_ext);

        this->matrix->add_frictionless_part2(mesh, mesh->node_positions[0], u_ext, lambda, D_cache, topopt);
        E = 0;
        for(size_t i = 0; i < vec_size; ++i){
            E += u[i]*(Ku[i]/2 - f[i]);
        }

        logger::quick_log("lambda n");
        const size_t l_num_orig = l_num/3;
        for(size_t i = 0; i < l_num_orig; ++i){
            std::cout << lambda[i] << " ";
        }
        std::cout << std::endl;

        logger::quick_log("Iteration:", it);
        logger::quick_log("Potential energy:", E);
        logger::quick_log("||r|| expected:", rnorm);
        rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
        logger::quick_log("||r|| actual:", rnorm);
        logger::quick_log("Step:", step);
        double max_gp = -1e100;
        double min_gp = 1e100;
        for(const auto& e:mesh->paired_boundary){
            for(size_t i = 0; i < bnum; ++i){
                const auto n1 = e.b1->nodes[i];
                const auto n2 = e.b2->nodes[i];
                const auto normal = e.b2->normal;
                double gp = 0;
                for(size_t j = 0; j < dof; ++j){
                    double u1 = 0;
                    double u2 = 0;
                    u1 = u_ext[n1->u_pos[j]];
                    u2 = u_ext[n2->u_pos[j]];
                    gp += (u2 - u1)*normal.Coord(1+j);
                }
                max_gp = std::max(gp, max_gp);
                min_gp = std::min(gp, min_gp);
            }
        }
        logger::quick_log("max_gp", max_gp);
        logger::quick_log("min_gp", min_gp);
        logger::quick_log("");
        //if(it == 10) exit(0);

        logger::log_assert(rnorm < 1e30, logger::ERROR, "Newton process diverged, increase rtol_abs and/or decrease mult");
        ++it;
    } while((rnorm > this->rtol_abs || it < 1) && step > 0);

    std::copy(u.begin(), u.begin() + u_size, load.begin());
}
void FiniteElement::solve_frictionless_displ_simple(const Meshing* const mesh, std::vector<double>& load, std::vector<double>& lambda, const std::vector<double>& u0){
    const size_t vec_size = this->u_size + this->l_num;
    const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();
    const size_t kw = mesh->elem_info->get_k_dimension();
    const auto& node_positions = mesh->node_positions[0];
    (void) u0;

    double E = 0;
    size_t it = 0;

    std::vector<double> f(vec_size, 0);
    std::vector<double> u(vec_size, 0);
    std::vector<double> r(vec_size);
    std::vector<double> G(vec_size);
    std::vector<double> HG(vec_size);
    std::vector<double> r2(vec_size);
    std::vector<double> dr(vec_size);
    std::vector<double> u_test(vec_size);
    std::vector<double> u1(vec_size);
    std::vector<double> u2(vec_size);
    std::vector<double> v0(vec_size);
    std::vector<double> v1(vec_size);
    std::vector<double> v2(vec_size);
    std::vector<double> UL(vec_size);
    std::vector<double> uold1(vec_size);
    std::vector<double> uold2(vec_size);

    //std::copy(u0.begin(), u0.end(), u.begin());

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

    // // Initiate interior point
    // for(const auto& e:mesh->paired_boundary){
    //     for(size_t i = 0; i < bnum; ++i){
    //         const auto n2 = e.b2->nodes[i];
    //         const auto normal = e.b1->normal;
    //         for(size_t j = 0; j < dof; ++j){
    //             auto ni2 = mesh->node_positions[0][n2->u_pos[j]];
    //             if(std::abs(normal.Coord(1+j)) > 1e-10){
    //                 const double sign = (normal.Coord(1+j) > 0) ? 1 : -1;
    //                 u[ni2] = -sign*1e-5;
    //             }
    //         }

    //         // If I ever start with non zero `u`, edit this
    //         /*
    //         for(size_t j = 0; j < dof; ++j){
    //             auto ni1 = mesh->node_positions[0][n1->u_pos[j]];
    //             auto ni2 = mesh->node_positions[0][n2->u_pos[j]];
    //             if(ni1 > -1){
    //                 du1 = dr[ni1]*normal.Coord(1+j);
    //             }
    //             if(ni2 > -1){
    //                 du2 = dr[ni2]*normal.Coord(1+j);
    //             }
    //         }
    //         gp = (du2 - du1);
    //         if(gp < -1e-10){
    //             if(std::abs(du2) < std::abs(du1)){
    //                 double a = du2/du1;
    //                 for(size_t j = 0; j < dof; ++j){
    //                     auto ni1 = mesh->node_positions[0][n1->u_pos[j]];
    //                     if(ni1 > -1 && std::abs(normal.Coord(1+j)) > 1e-10){
    //                         dr[ni1] *= a;
    //                     }
    //                 }
    //             } else {
    //                 double a = du1/du2;
    //                 for(size_t j = 0; j < dof; ++j){
    //                     auto ni2 = mesh->node_positions[0][n2->u_pos[j]];
    //                     if(ni2 > -1 && std::abs(normal.Coord(1+j)) > 1e-10){
    //                         dr[ni2] *= a;
    //                     }
    //                 }
    //             }
    //         }
    //         */
    //     }
    // }

    mesh->extend_vector(0, load, f_ext);
    mesh->extend_vector(0, u, u_ext);

    std::copy(lambda.begin(), lambda.end(), u.begin() + u_size);
    double rnorm = 0;
    //double old_rnorm = 0;

    double step = 1.0;
    double LAG = this->matrix->get_lag_displ_simple();
    //double gamma = 0;

    std::vector<gp_Pnt> points(bnum);
    std::vector<long> u_pos(2*kw);
    std::vector<long> l_pos(bnum);
    std::vector<double> u_e1(kw), u_e2(kw);
    std::vector<double> l_e(bnum);
    std::vector<double> du_e1(kw), du_e2(kw);
    std::vector<double> dl_e(bnum);

    typedef Eigen::SparseMatrix<double, Eigen::ColMajor, std::ptrdiff_t> Mat;
    typedef Eigen::VectorXd Vec;

    size_t rel_id = 0;
    std::unordered_map<size_t, size_t> global_to_local_id;
    {
        std::set<size_t> global_to_local_id_set;
        for(const auto& e:mesh->paired_boundary){
            for(size_t n = 0; n < node_num; ++n){
                global_to_local_id_set.insert(e.b1->parent->nodes[n]->id);
                global_to_local_id_set.insert(e.b2->parent->nodes[n]->id);
            }
        }
        for(const auto i:global_to_local_id_set){
            global_to_local_id[i] = rel_id;
            ++rel_id;
        }
    }

    const size_t ls_max = mesh->lag_node_map.size();
    const size_t uls_max = rel_id*dof*2;
    Mat Gl(uls_max, ls_max);
    Vec Ll(ls_max);
    Vec Rl(ls_max);
    Vec Rl0(uls_max);
    math::Matrix uL(2*kw, bnum, 0);
    for(const auto& e:mesh->paired_boundary){
        //for(size_t i = 0; i < bnum; ++i){
        //    points[i] = e.b1->nodes[i]->point;
        //}
        for(size_t i = 0; i < node_num; ++i){
            const auto n1 = e.b1->parent->nodes[i];
            const auto n2 = e.b2->parent->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                u_pos[dof*i + j] = global_to_local_id[n1->id]*dof + j;
                u_pos[dof*(node_num + i) + j] = global_to_local_id[n2->id]*dof + j;
            }
        }
        for(size_t i = 0; i < bnum; ++i){
            const auto b1_id = e.b1->nodes[i]->id;
            l_pos[i] = mesh->lag_node_map.at(b1_id);
            //l_e[i] = lambda[l_pos[i]];
            //l_pos[i] += u_size;
        }
        //for(size_t i = 0; i < node_num; ++i){
        //    const auto n1 = e.b1->parent->nodes[i];
        //    const auto n2 = e.b2->parent->nodes[i];
        //    for(size_t j = 0; j < dof; ++j){
        //        u1[dof*i + j] = u_ext[n1->u_pos[j]];
        //        u2[dof*i + j] = u_ext[n2->u_pos[j]];
        //    }
        //}
        uL = e.elem->fl2_uL(l_e, u_e1, u_e2);
        for(size_t i = 0; i < 2*kw; ++i){
            for(size_t j = 0; j < bnum; ++j){
                const double m = uL(i, j);
                if(std::abs(m) > 1e-15){
                    if(u_pos[i] >= 0){
                        Gl.coeffRef(u_pos[i], l_pos[j]) = m;
                    }
                }
            }
        }
    }

    const double L_MIN = -200;
    const double L_MAX = 0;
    const double L_STEP = 10;

    const auto apply_constraint = [&](){
        for(const auto& e:mesh->paired_boundary){
            for(size_t i = 0; i < node_num; ++i){
                const auto n1 = e.b1->parent->nodes[i];
                const auto n2 = e.b2->parent->nodes[i];
                for(size_t j = 0; j < dof; ++j){
                    const auto u_id1 = node_positions[dof*n1->id+i];
                    Rl0[global_to_local_id[n1->id]*dof + j] = -r[u_id1];
                    const auto u_id2 = node_positions[dof*n2->id+i];
                    Rl0[global_to_local_id[n2->id]*dof + j] = -r[u_id2];
                }
            }
        }
        std::vector<double> LS(ls_max, L_STEP);

        //Ll.fill(0);
        auto Ll_old = Ll;
        auto Ll_old2 = Ll;
        double dl = 1;
        Rl = (Rl0 + Gl*Ll).transpose()*Gl;
        logger::quick_log("Rl0", std::sqrt(Rl0.transpose()*Rl0));
        logger::quick_log("Rl1", std::sqrt(Rl.transpose()*Rl));
        while(dl > 1e-5){
            dl = 0;
            Rl = (Rl0 + Gl*Ll).transpose()*Gl;

            for(size_t i = 0; i < ls_max; ++i){
                const double s = LS[i];
                const double U = std::min(Ll[i] + s, L_MAX);
                const double UX = U - Ll[i];
                const double L = std::max(Ll[i] - s, L_MIN);
                const double LX = Ll[i] - L;
                const double dgdxp = std::max(0.0,  Rl[i]);
                const double dgdxm = std::max(0.0, -Rl[i]);
                const double P = (UX*UX)*(1.001*dgdxp + 0.001*dgdxm) + 1e-13;
                const double Q = (LX*LX)*(0.001*dgdxp + 1.001*dgdxm) + 1e-13;
                const bool P_nonzero = std::abs(P) > 1e-15;
                const bool Q_nonzero = std::abs(Q) > 1e-15;
                logger::log_assert(P_nonzero && Q_nonzero, logger::ERROR, "P or Q became zero somehow");
                //double R = 0;
                //if(UX > 0){
                //    R -= P/UX;
                //}
                //if(LX > 0){
                //    R -= Q/LX;  
                //}
                const double sqp = std::sqrt(P);
                const double sqq = std::sqrt(Q);
                Ll[i] = (U*sqq + L*sqp)/(sqp + sqq);

                if((Ll[i] - Ll_old[i])*(Ll_old[i] - Ll_old2[i]) < 0){
                    LS[i] *= 0.7;
                }
            }

            for(size_t i = 0; i < ls_max; ++i){
                dl = std::max(dl, std::abs(Ll[i] - Ll_old[i]));
                Ll_old2[i] = Ll_old[i];
                Ll_old[i] = Ll[i];
            }
            //for(auto& li:Ll){
            //    std::cout << li << " ";
            //    //if(li < -1e-10){
            //    //    li -= 10;
            //    //}
            //}
            //std::cout << std::endl;
            //logger::quick_log(dl);
        }
        logger::quick_log("Rl2", std::sqrt(Rl.transpose()*Rl));
        for(size_t i = 0; i < ls_max; ++i){
            std::cout << Ll[i] << " ";
            //Ll[i] *= 5;
        }
        std::cout << std::endl;
        const Eigen::VectorXd result = Gl*Ll;
        for(const auto& pair:global_to_local_id){
            const auto n_id = pair.first;
            //const auto l_id = pair.second;
            for(size_t i = 0; i < dof; ++i){
                const auto u_id = node_positions[n_id];
                if(u_id >= 0){
                    r[u_id] -= result[n_id];
                }
            }
        }
        // for(const auto& e:mesh->paired_boundary){
        //     for(size_t i = 0; i < node_num; ++i){
        //         const auto n1 = e.b1->parent->nodes[i];
        //         const auto n2 = e.b2->parent->nodes[i];
        //         for(size_t j = 0; j < dof; ++j){
        //             u_pos[dof*i + j] = node_positions[n1->u_pos[j]];
        //             u_pos[dof*(node_num + i) + j] = node_positions[n2->u_pos[j]];
        //         }
        //     }
        //     for(size_t i = 0; i < bnum; ++i){
        //         l_pos[i] = mesh->lag_node_map.at(e.b1->nodes[i]->id);
        //     }
        //     uL = e.elem->fl2_uL(l_e, u_e1, u_e2);
        //     for(size_t i = 0; i < 2*kw; ++i){
        //         if(u_pos[i] >= 0){
        //             for(size_t j = 0; j < bnum; ++j){
        //                 r[u_pos[i]] -= uL(i, j)*Ll[l_pos[j]];
        //             }
        //         }
        //     }
        // }
    };


    this->matrix->dot_vector(u, Ku);
    this->matrix->append_Ku_frictionless_simple(mesh, u, Ku);
    for(size_t i = 0; i < vec_size; ++i){
        r[i] = -(Ku[i] - f[i]);
    }
    this->matrix->add_frictionless_simple(mesh, mesh->node_positions[0], u_ext, lambda);

    const double u_max =  0.1;
    const double u_min = -0.1;
    const double l_max =  100;
    const double l_min =  0;
    //const double mma_step = 0.01;
    const double mma_step = 0.0001;
    const double mma_step_l = 10;
    //const double K_L = 10;
    const double K_L = 10;
    std::vector<double> S(vec_size, mma_step);
    std::fill(S.begin() + u_size, S.end(), mma_step_l);

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
        logger::quick_log("!!!!!!!!!!!!!!!!!!!!!!", u_size, vec_size);
        logger::quick_log("dr min max constr", *std::min_element(dr.begin()+u_size, dr.end()), *std::max_element(dr.begin()+u_size, dr.end()));


        const auto update_vars = [&](const size_t i, const double min, const double max){
            //const double s = std::min(std::abs(dr[i]), 1.0)*S[i];
            const double s = S[i];
            const double U = std::min(u[i] + s, max);
            const double UX = U - u[i];
            const double L = std::max(u[i] - s, min);
            const double LX = u[i] - L;
            const double P =  (UX*UX*UX)/(UX + LX) + 1e-13/(max - min);
            const double Q = -(LX*LX*LX)/(UX + LX) - 1e-13/(max - min);
            const double V = dr[i];
            const bool P_nonzero = std::abs(P) > 1e-15;
            const bool Q_nonzero = std::abs(Q) > 1e-15;
            logger::log_assert(P_nonzero && Q_nonzero, logger::ERROR, "P or Q became zero somehow");
            double R = 0;
            if(UX > 0){
                R -= P/UX;
            }
            if(LX > 0){
                R -= Q/LX;  
            }

            const double a = -(V - R);
            const double b = ((V - R)*(U + L) - (P - Q));
            const double c = -(V - R)*U*L - Q*U + P*L;

            double delta = b*b - 4*a*c;
            if(delta >= 0 && std::abs(a) > 1e-15){
                const double xt1 = (-b + std::sqrt(delta))/(2*a);
                const double xt2 = (-b - std::sqrt(delta))/(2*a);
                //const bool xt1_within = (L < xt1 || std::abs(xt1 - L) < 1e-5) && (xt1 < U || std::abs(U - xt1) < 1e-5);
                //const bool xt2_within = (L < xt2 || std::abs(xt2 - L) < 1e-5) && (xt2 < U || std::abs(U - xt2) < 1e-5);
                if((V-R) > 0){
                    u[i] = std::max(xt1, xt2);
                    u[i] = std::min(u[i], U);
                } else if(std::abs(V) > 1e-15){
                    u[i] = std::min(xt1, xt2);
                    u[i] = std::max(u[i], L);
                }
                //if(std::abs(V) > 0.001){
                //if(!xt1_within && !xt2_within){
                //    if(std::abs(xt1 - L) < 1e-10 || std::abs(xt2 - L) < 1e-10){
                //        u[i] = L;
                //    } else if(std::abs(xt1 - U) < 1e-10 || std::abs(xt2 - U) < 1e-10) {
                //        u[i] = U;
                //    }
                //} else {
                //    logger::log_assert(xt1_within != xt2_within || std::abs(xt1 - xt2) < 1e-5, logger::ERROR, "incorrect assumption about test points placement, L: {} xt1: {} xt2: {} U: {} xi: {} P: {} Q: {} R: {}, V: {}", L, xt1, xt2, U, u[i], P, Q, R, V);
                //    if(xt1_within){
                //        u[i] = xt1;
                //    } else {
                //        u[i] = xt2;
                //    }
                //}
            }
            // Otherwise, x == x0
            // logger::log_assert(L - 1e-15 < u[i] && u[i] < U + 1e-15, logger::ERROR, "Update violated asymptotes, L: {} U: {} xi: {} P: {} Q: {} R: {} V: {} P_nonzero: {} Q_nonzero: {}", L, U, u[i], P, Q, R, V, P_nonzero, Q_nonzero);
            // logger::log_assert(min - 1e-15 < u[i] && u[i] < max + 1e-15, logger::ERROR, "Update violated boundaries, L: {} U: {} xi: {} P: {} Q: {} R: {} V: {} P_nonzero: {} Q_nonzero: {}", L, U, u[i], P, Q, R, V, P_nonzero, Q_nonzero);
            
            if(it > 1){
                if((uold2[i] - uold1[i])*(uold1[i] - u[i]) < 0){
                    S[i] *= 0.7;
                    //S[i] = std::abs(uold1[i] - u[i]);
                //} else {
                //    S[i] = std::min(1.05*S[i], mma_step);
                }
            }
            uold2[i] = uold1[i];
            uold1[i] = u[i];
        };
        
        #pragma omp parallel
        {
            #pragma omp for
            for(size_t i = 0; i < u_size; ++i){
                update_vars(i, u_min, u_max);
            }
            //#pragma omp for
            //for(size_t i = u_size; i < vec_size; ++i){
            //    S[i] = mma_step_l;
            //    update_vars(i, l_min, l_max);
            //}
        }

        //for(size_t i = 0; i < vec_size; ++i){
        for(size_t i = 0; i < u_size; ++i){
            u_var = std::max(std::abs(u[i] - uold2[i]), u_var);
        }


        step = 1;
        
        this->reset_hessian();

        double drnorm = std::sqrt(cblas_ddot(dr.size(), dr.data(), 1, dr.data(), 1));
        logger::quick_log("drnorm", drnorm);
        logger::quick_log("dr min max", *std::min_element(dr.begin(), dr.end()), *std::max_element(dr.begin(), dr.end()));

        // 2-point MMA
        /*
            const double eta_0 = 0;
            const double eta_1 = this->max_step;
            const double L = -0.1*eta_1; 
            const double U = 1.1*eta_1;

            double g1 = 0, dg1 = 0;
            double g2 = 0, dg2 = 0;
            double LAG_1 = 0, LAG_2 = 0;

            for(size_t i = 0; i < u.size(); ++i){
                u1[i] = u[i] + eta_0*dr[i];
                u2[i] = u[i] + eta_1*dr[i];
            }
            //if(dLAG > 0){
                LAG_1 = LAG + eta_0*dLAG;
                LAG_2 = LAG + eta_1*dLAG;
            //} else {
            //    LAG_1 = LAG;
            //    LAG_2 = LAG;
            //    dLAG = 0;
            //}
            get_g_dg(g1, dg1, eta_0, u1, dr);
            get_g_dg(g2, dg2, eta_1, u2, dr);

            std::fill(Ku.begin(), Ku.end(), 0);
            std::fill(Ku2.begin(), Ku2.end(), 0);
            std::fill(Kd1.begin(), Kd1.end(), 0);
            std::fill(Kd2.begin(), Kd2.end(), 0);

            this->matrix->dot_vector(u1, Ku);
            this->matrix->dot_vector(u2, Ku2);
            this->matrix->dot_vector(dr, Kd1);
            std::copy(Kd1.begin(), Kd1.end(), Kd2.begin());

            this->matrix->set_lag_displ_simple(1);

            std::fill(Ku_tmp.begin(), Ku_tmp.end(), 0);
            this->matrix->append_Ku_frictionless_simple(mesh, u1, Ku_tmp);
            this->matrix->append_dKu_frictionless_simple(mesh, u1, dr, eta_0, Kd1);
            for(size_t i = 0; i < vec_size; ++i){
                Ku[i] += LAG_1*Ku_tmp[i];
                Kd1[i] += dLAG*Ku_tmp[i];
            }

            std::fill(Ku_tmp.begin(), Ku_tmp.end(), 0);
            this->matrix->append_Ku_frictionless_simple(mesh, u2, Ku_tmp);
            this->matrix->append_dKu_frictionless_simple(mesh, u2, dr, eta_1, Kd2);
            for(size_t i = 0; i < vec_size; ++i){
                Ku2[i] += LAG_2*Ku_tmp[i];
                Kd2[i] += dLAG*Ku_tmp[i];
            }

            for(size_t i = 0; i < vec_size; ++i){
                r[i] = (Ku[i] - f[i]);
                r2[i] = (Ku2[i] - f[i]);
            }
            rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1) + g1*g1);
            const double rnorm2 = std::sqrt(cblas_ddot(r2.size(), r2.data(), 1, r2.data(), 1) + g2*g2);

            const double b = 2*(cblas_ddot(r.size(), r.data(), 1, Kd1.data(), 1) + g1*dg1)/rnorm;
            const double a = 2*(cblas_ddot(r2.size(), r2.data(), 1, Kd2.data(), 1) + g2*dg2)/rnorm2;
            const double dU0 = U - eta_0;
            const double dL0 = eta_0 - L;
            const double dU1 = U - eta_1;
            const double dL1 = eta_1 - L;
            math::Matrix PQM({1.0/(dU0*dU0), -1.0/(dL0*dL0),
                              1.0/(dU1*dU1), -1.0/(dL1*dL1)}, 2, 2);
            math::Vector PQV({b, a});
            math::LU PQLU(PQM);
            PQLU.solve(PQV);
            const double P = PQV[0];
            const double Q = PQV[1];
            logger::quick_log("a b P Q", a, b, P, Q);
            logger::log_assert(P >= 0 || Q >= 0,
                    logger::WARNING,
                    "P and Q are less than zero, there might be a mistake in the problem's definition");
            if(P < 0){
                // Minimum is beyond the approximation's bounds
                step = eta_1;
            } else if(Q < 0){
                // Hessian is not SPD? Meaning this is a saddle point?
                step = 0;
            } else {
                // Minimum is within the approximation's bounds
                const double sqrtP = std::sqrt(P);
                const double sqrtQ = std::sqrt(Q);
                step = (U*sqrtQ + L*sqrtP)/(sqrtP + sqrtQ);
            }
            logger::quick_log("step", step);
            
            for(size_t i = 0; i < u.size(); ++i){
                u[i] += step*dr[i];
            }
            LAG += step*dLAG;
            mesh->extend_vector(0, u, u_ext);

            std::fill(Ku.begin(), Ku.end(), 0);
            this->matrix->dot_vector(u, Ku);
            std::fill(Ku.begin() + u_size, Ku.end(), 0);
            this->matrix->set_lag_displ_simple(LAG);
            this->matrix->append_Ku_frictionless_simple(mesh, u, Ku);
            for(size_t i = 0; i < vec_size; ++i){
                r[i] = -(Ku[i] - f[i]);
            }
            get_g_dg(g1, dg1, step, u, dr);
            rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1) + g1*g1);
            gamma = g1;
        */
        
        // Damped Newton-Rhapson
        /*
            const auto get_g = [&](double alpha)->double{
                for(size_t i = 0; i < u.size(); ++i){
                    u1[i] = u[i] + alpha*dr[i];
                }
                // double LAG = 0;
                // for(size_t i = u_size; i < vec_size; ++i){
                //     if(u1[i] < LAG){
                //         LAG = u1[i];
                //     }
                // }
                // LAG = 1e-2*std::exp(-LAG);
                // this->matrix->set_lag_displ_simple(LAG);

                std::fill(Ku.begin(), Ku.end(), 0);
                std::fill(Kd1.begin(), Kd1.end(), 0);
                this->matrix->dot_vector(u1, Ku);
                this->matrix->dot_vector(dr, Kd1);
                this->matrix->append_Ku_frictionless_simple(mesh, u1, Ku);
                this->matrix->append_dKu_frictionless_simple(mesh, u1, dr, alpha, Kd1);
                for(size_t i = 0; i < vec_size; ++i){
                    r[i] = (Ku[i] - f[i]);
                }

                rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));// + gamma*gamma);
                //double gamma = 0, dgamma = 0;
                //get_g_dg(gamma, dgamma, alpha, u1, dr);
                //if(dLAG != 0){
                //    return cblas_ddot(vec_size, Kd1.data(), 1, r.data(), 1) + dgamma*gamma;
                //} else {
                    return cblas_ddot(vec_size, Kd1.data(), 1, r.data(), 1)/rnorm;
                //}
            };
            const auto get_g_old = [&](double alpha)->double{
                for(size_t i = 0; i < u.size(); ++i){
                    u1[i] = u[i] + alpha*dr[i];
                }

                std::fill(Ku.begin(), Ku.end(), 0);
                this->matrix->dot_vector(u1, Ku);
                for(size_t i = 0; i < vec_size; ++i){
                    r[i] = (Ku[i] - f[i]);
                }
                const double dE = cblas_ddot(r.size(), r.data(), 1, dr.data(), 1);
                const double dg = get_dg_PI(alpha, u1, dr); 

                return dE + dg;
            };
            const auto init_step = [&]()->double{
                double alpha = this->max_step;
                for(const auto& e:mesh->paired_boundary){
                    for(size_t i = 0; i < bnum; ++i){
                        const auto n1 = e.b1->nodes[i];
                        const auto n2 = e.b2->nodes[i];
                        const auto normal = e.b1->normal;
                        double gp = 0;
                        double dgp = 0;
                        for(size_t j = 0; j < dof; ++j){
                            auto ni1 = mesh->node_positions[0][n1->u_pos[j]];
                            auto ni2 = mesh->node_positions[0][n2->u_pos[j]];
                            double u1 = 0;
                            double u2 = 0;
                            double du1 = 0;
                            double du2 = 0;
                            if(ni1 > -1){
                                u1 = u[ni1];
                                du1 = dr[ni1];
                            }
                            if(ni2 > -1){
                                u2 = u[ni2];
                                du2 = dr[ni2];
                            }
                            gp += (u2 - u1)*normal.Coord(1+j);
                            dgp += (du2 - du1)*normal.Coord(1+j);
                        }
                        double test_alpha = -gp/dgp;
                        if(test_alpha > 0 && test_alpha < alpha){
                            alpha = 0.5*test_alpha;
                        }
                    }
                }
                return alpha;
            };

            //get_LAG();
            double g1 = 0, g2 = 1;
            double a1 = 0, a2 = this->max_step;//init_step();//this->max_step;
            //if(drnorm > 1e3){
            //    a2 = this->max_step/std::pow(drnorm, 3.0/2.0);
            //}
            g1 = get_g(a1);
            double g0 = g1;
            g2 = get_g(a2);
            logger::quick_log("g0", g0);
            logger::quick_log("g2", g2);
            logger::quick_log("Improving starting points...");
            if(g0 > 0){
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
            if(a1 > -1){
                if(std::abs(g2/g1) > 1e4 && g2*g1 < 0){
                    a2 = std::abs(g1/g2);
                    g2 = get_g(a2);
                    while(g2*g1 > 0){
                        a2 *= 2;
                        g2 = get_g(a2);
                    }
                    logger::quick_log(a2, g2);
                }
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
            for(size_t i = 0; i < u.size(); ++i){
                u[i] += step*dr[i];
            }
            for(size_t i = u_size; i < vec_size; ++i){
                if(std::abs(u[i]) > 0.1){
                    if(u[i] > 0){
                        u[i] = 0.1;
                    } else {
                        u[i] = -0.1;
                    }
                }
            }
        */

            std::fill(Ku.begin(), Ku.end(), 0);
            this->matrix->dot_vector(u, Ku);
            E = 0;
            // This is definitely wrong, but works as an estimate
            for(size_t i = 0; i < vec_size; ++i){
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
                if(gp < 0){
                    u[u_size + mesh->lag_node_map.at(n1->id)] -= gp*K_L;
                }
                max_gp = std::max(gp, max_gp);
                min_gp = std::min(gp, min_gp);
            }
        }
        logger::quick_log("max_gp ", max_gp);
        logger::quick_log("min_gp", min_gp);

        logger::log_assert(rnorm < 1e30, logger::ERROR, "Newton process diverged, increase rtol_abs and/or decrease mult");
        ++it;
        //if(b > 1e-6 || step <= 0){
        //    break;
        //}
        if(step == 0){
            break;
        }
    } while((rnorm > this->rtol_abs && std::abs(step) > this->step_tol && u_var > 1e-6) || it < 1);

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
                    l[mesh->lag_node_map.at(n1->id)] += 1e-3*gp*this->start_lag_simple;
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

    const double u_max =  0.1;
    const double u_min = -0.1;
    const double mma_step = 0.01;
    const double mma_step_ctc = 0.0001;
    std::vector<double> S(vec_size, mma_step);

    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < bnum; ++i){
            const auto n1 = e.b1->nodes[i];
            const auto n2 = e.b2->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                const auto ni1 = mesh->node_positions[0][n1->u_pos[j]];
                const auto ni2 = mesh->node_positions[0][n2->u_pos[j]];
                is_contact_node[ni1] = true;
                is_contact_node[ni2] = true;
                S[ni1] = mma_step_ctc;
                S[ni2] = mma_step_ctc;
            }
        }
    }

    const auto update_vars = [&](const size_t i, const double min, const double max){
        const double s = S[i];
        const double U = std::min(u[i] + s, max);
        const double UX = U - u[i];
        const double L = std::max(u[i] - s, min);
        const double LX = u[i] - L;
        const double P =  (UX*UX*UX)/(UX + LX) + 1e-13;///(max - min);
        const double Q = -(LX*LX*LX)/(UX + LX) - 1e-13;///(max - min);
        const double V = dr[i];
        const bool P_nonzero = std::abs(P) > 1e-15;
        const bool Q_nonzero = std::abs(Q) > 1e-15;
        logger::log_assert(P_nonzero && Q_nonzero, logger::ERROR, "P or Q became zero somehow");
        double R = 0;
        if(UX > 0){
            R -= P/UX;
        }
        if(LX > 0){
            R -= Q/LX;  
        }

        const double a = -(V - R);
        const double b = ((V - R)*(U + L) - (P - Q));
        const double c = -(V - R)*U*L - Q*U + P*L;

        double delta = b*b - 4*a*c;
        if(delta >= 0 && std::abs(a) > 1e-15){
            const double xt1 = (-b + std::sqrt(delta))/(2*a);
            const double xt2 = (-b - std::sqrt(delta))/(2*a);
            const bool xt1_within = L < xt1 && xt1 < U;
            const bool xt2_within = L < xt2 && xt2 < U;
            //if(std::abs(V) > 0.001){
            if(!xt1_within && !xt2_within){
                if(std::abs(xt1 - L) < 1e-15 || std::abs(xt2 - L) < 1e-15){
                    u[i] = L;
                } else if(std::abs(xt1 - U) < 1e-15 || std::abs(xt2 - U) < 1e-15) {
                    u[i] = U;
                }
            } else {
                //logger::log_assert(xt1_within != xt2_within || std::abs(xt1 - xt2) < 1e-5, logger::ERROR, "incorrect assumption about test points placement, L: {} xt1: {} xt2: {} U: {} xi: {} P: {} Q: {} R: {}, V: {}", L, xt1, xt2, U, u[i], P, Q, R, V);
                if(xt1_within != xt2_within){
                    if(xt1_within){
                        u[i] = xt1;
                    } else {
                        u[i] = xt2;
                    }
                } else {
                    if(-r[i] > 0){
                        u[i] = std::min(xt1, xt2);
                    } else {
                        u[i] = std::max(xt1, xt2);
                    }
                }
            }
        }
        // Otherwise, x == x0
        logger::log_assert(L - 1e-15 < u[i] && u[i] < U + 1e-15, logger::ERROR, "Update violated asymptotes, L: {} U: {} xi: {} P: {} Q: {} R: {} V: {} P_nonzero: {} Q_nonzero: {}", L, U, u[i], P, Q, R, V, P_nonzero, Q_nonzero);
        logger::log_assert(min - 1e-15 < u[i] && u[i] < max + 1e-15, logger::ERROR, "Update violated boundaries, L: {} U: {} xi: {} P: {} Q: {} R: {} V: {} P_nonzero: {} Q_nonzero: {}", L, U, u[i], P, Q, R, V, P_nonzero, Q_nonzero);
        
        if(it > 1){
            if((uold2[i] - uold1[i])*(uold1[i] - u[i]) < 0){
                S[i] *= 0.7;
                //S[i] = std::abs(uold1[i] - u[i]);
            } else {
                //S[i] = std::min(1.05*S[i], mma_step);
            }
        }
        uold2[i] = uold1[i];
        uold1[i] = u[i];
    };

    double u_var = 0;
    step = this->max_step;
    
    do{
        u_var = 0;

        // Solve linear subproblem
        std::copy(r.begin(), r.end(), dr.begin());
        this->solve(dr);

        //#pragma omp parallel
        //{
        //    #pragma omp for
        //    for(size_t i = 0; i < u_size; ++i){
        //        update_vars(i, u_min, u_max);
        //    }
        //}

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
            //step = this->max_step;
            //for(size_t i = 0; i < u.size(); ++i){
            //    u[i] += this->max_step*dr[i];
            //}
            apply_multiplier(u, lambda);
            


        for(size_t i = 0; i < vec_size; ++i){
            u_var = std::max(std::abs(u[i] - uold2[i]), u_var);
            uold2[i] = u[i];
        }

        //for(size_t i = 0; i < u.size(); ++i){
        //    u[i] += step*dr[i];
        //}
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
        logger::quick_log("du", u_var);
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
    } while((rnorm > this->rtol_abs || it < 1) && step > 0 && u_var > 1e-8 && it < 50);

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

void FiniteElement::apply_constr_force(const Meshing* mesh, const std::vector<double>& u_ext, const std::vector<double>& lambda, std::vector<double>& b) const{
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    const size_t num_bound_nodes = mesh->elem_info->get_boundary_nodes_per_element();
    const size_t dof = mesh->elem_info->get_dof_per_node();
    const size_t kw = mesh->elem_info->get_k_dimension();
    const size_t l_num = mesh->lag_node_map.size();

    math::Vector ln_e(num_bound_nodes);
    math::Vector lp1_e(num_bound_nodes);
    math::Vector lp2_e(num_bound_nodes);
    math::Vector u_e(kw);
    std::vector<long> u_pos(kw);
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < num_bound_nodes; ++i){
            const long p1 = mesh->lag_node_map.at(e.b1->nodes[i]->id);
            const long p2 = p1 + l_num;
            const long p3 = p2 + l_num;
            ln_e[i] = lambda[p1];
            lp1_e[i] = lambda[p2];
            lp2_e[i] = lambda[p3];
        }
        for(size_t i = 0; i < num_nodes; ++i){
            const auto n1 = e.elem->e1->nodes[i];
            const auto n2 = e.elem->e2->nodes[i];
            const long id = this->constr_id_map.at(n2);
            if(id >= 0){
                for(size_t j = 0; j < dof; ++j){
                    u_pos[dof*i + j] = id*dof + j;
                    u_e[dof*i + j] = u_ext[n1->u_pos[j]];
                }
            } else {
                for(size_t j = 0; j < dof; ++j){
                    u_pos[dof*i + j] = -1;
                    u_e[dof*i + j] = 0;
                }
            }
        }
        const auto be(e.elem->fl3_eq(ln_e, lp1_e, lp2_e, u_e));
        for(size_t i = 0; i < kw; ++i){
            if(u_pos[i] >= 0){
                b[u_pos[i]] += be[i];
            }
        }
    }
}
void FiniteElement::apply_lambda(const Meshing* mesh, const std::vector<double>& lambda, std::vector<double>& u_ext){
    const size_t dof = mesh->elem_info->get_dof_per_node();

    std::fill(this->constr_force.begin(), this->constr_force.end(), 0);
    this->apply_constr_force(mesh, u_ext, lambda, this->constr_force);
    this->constr_solver.solve(this->constr_force);
    for(const auto n:this->boundary_reference_nodes){
        const long id = this->constr_id_map.at(n);
        for(size_t j = 0; j < dof; ++j){
            u_ext[n->u_pos[j]] += this->constr_force[id*dof + j];
        }
    }
}
