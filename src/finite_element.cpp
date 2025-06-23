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
#include "finite_element.hpp"
#include "logger.hpp"
#include "global_stiffness_matrix.hpp"
#include "math/matrix.hpp"

FiniteElement::FiniteElement(ContactType contact_type, double rtol_abs, double max_step)
    :contact_type(contact_type), rtol_abs(rtol_abs), max_step(max_step){

}

void FiniteElement::generate_matrix(const Meshing* const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda){
    this->u_size = u_size;
    this->l_num = l_num;
    this->generate_matrix_base(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, lambda, this->contact_type);

    if(this->contact_type == FRICTIONLESS_DISPL_CONSTR){
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
            this->solve_frictionless_penalty(mesh, load, topopt, D_cache, u0);
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
void FiniteElement::solve_frictionless_displ(const Meshing* const mesh, std::vector<double>& load, std::vector<double>& lambda, const bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u0){
    const size_t vec_size = this->u_size + this->l_num;

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
        const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
        const size_t dof = mesh->elem_info->get_dof_per_node();
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
    const size_t dof = mesh->elem_info->get_dof_per_node();
    const size_t kw = mesh->elem_info->get_k_dimension();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
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

    std::copy(lambda.begin(), lambda.end(), u.begin() + u_size);
    double rnorm = 0;
    //double old_rnorm = 0;

    double step = 1.0;
    double LAG = this->matrix->get_lag_displ_simple();
    double gamma = 0;

    std::vector<gp_Pnt> points(bnum);
    std::vector<long> u_pos(2*kw);
    std::vector<long> l_pos(bnum);
    std::vector<double> u_e1(kw), u_e2(kw);
    std::vector<double> l_e(bnum);
    std::vector<double> du_e1(kw), du_e2(kw);
    std::vector<double> dl_e(bnum);
    auto& node_positions = mesh->node_positions[0];

    const auto get_g_dg_ddg = [&](double& gamma, double& dgamma, double& ddgamma, double eta, const std::vector<double>& u, const std::vector<double>& du){
        gamma = 0;
        dgamma = 0;
        
        for(const auto& e:mesh->paired_boundary){
            for(size_t i = 0; i < bnum; ++i){
                const size_t pos = mesh->lag_node_map.at(e.b1->nodes[i]->id);
                l_e[i] = u[pos + u_size];
                dl_e[i] = du[pos + u_size];
            }
            for(size_t i = 0; i < node_num; ++i){
                const auto n1 = e.b1->parent->nodes[i];
                const auto n2 = e.b2->parent->nodes[i];
                for(size_t j = 0; j < dof; ++j){
                    u_e1[dof*i + j] = u[node_positions[n1->u_pos[j]]];
                    u_e2[dof*i + j] = u[node_positions[n2->u_pos[j]]];
                    du_e1[dof*i + j] = du[node_positions[n1->u_pos[j]]];
                    du_e2[dof*i + j] = du[node_positions[n2->u_pos[j]]];
                }
            }
            const double integ(e.elem->fl2_int(l_e, u_e1, u_e2));
            const double integ_deriv(e.elem->fl2_int_deriv(l_e, u_e1, u_e2, dl_e, du_e1, du_e2, eta));
            const double integ_deriv2(e.elem->fl2_int_deriv2(l_e, u_e1, u_e2, dl_e, du_e1, du_e2, eta));

            gamma += integ;
            dgamma += integ_deriv;
            ddgamma += integ_deriv2;
        }
    };


    //this->matrix->dot_vector(u, Ku);
    //std::fill(Ku.begin() + u_size, Ku.end(), 0);
    this->matrix->append_Ku_frictionless_simple(mesh, u, Ku);
    for(size_t i = 0; i < vec_size; ++i){
        r[i] = -(Ku[i] - f[i]);
    }
    this->matrix->add_frictionless_simple(mesh, mesh->node_positions[0], u_ext, lambda);
    double unused = 0;
    get_g_dg_ddg(gamma, unused, unused, step, u, dr);
    step = 1;
    rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1) + gamma*gamma);

    do{
        std::fill(G.begin(), G.end(), 0);
        for(const auto& e:mesh->paired_boundary){
            for(size_t i = 0; i < bnum; ++i){
                points[i] = e.b1->nodes[i]->point;
            }
            for(size_t i = 0; i < node_num; ++i){
                const auto n1 = e.b1->parent->nodes[i];
                const auto n2 = e.b2->parent->nodes[i];
                for(size_t j = 0; j < dof; ++j){
                    u_pos[dof*i + j] = node_positions[n1->u_pos[j]];
                    u_pos[kw + dof*i + j] = node_positions[n2->u_pos[j]];
                }
            }
            for(size_t i = 0; i < bnum; ++i){
                l_pos[i] = mesh->lag_node_map.at(e.b1->nodes[i]->id);
                l_e[i] = lambda[l_pos[i]];
                l_pos[i] += u_size;
            }
            for(size_t i = 0; i < node_num; ++i){
                const auto n1 = e.b1->parent->nodes[i];
                const auto n2 = e.b2->parent->nodes[i];
                for(size_t j = 0; j < dof; ++j){
                    u_e1[dof*i + j] = u_ext[n1->u_pos[j]];
                    u_e2[dof*i + j] = u_ext[n2->u_pos[j]];
                }
            }
            const auto LAGuL(e.elem->fl2_LAG(l_e, u_e1, u_e2));

            for(size_t i = 0; i < 2*kw; ++i){
                G[u_pos[i]] += LAGuL[i];
            }
            for(size_t i = 0; i < bnum; ++i){
                G[l_pos[i]] += LAGuL[i + 2*kw];
            }
        }
        std::copy(G.begin(), G.end(), HG.begin());
        // Solve linear subproblem
        logger::quick_log("r min max", *std::min_element(r.begin(), r.end()), *std::max_element(r.begin(), r.end()));
        logger::quick_log("u min max", *std::min_element(u.begin(), u.end()), *std::max_element(u.begin(), u.end()));
        this->solve(r);
        this->solve(HG);

        double GHr = cblas_ddot(vec_size, G.data(), 1, r.data(), 1);
        double GHG = cblas_ddot(vec_size, G.data(), 1, HG.data(), 1);
        // `r` is already negative
        //double dLAG = (gamma + GHr)/GHG;
        double dLAG = 0;
        //double dLAG = 0.1*LAG;
        if(dLAG < 0){
            dLAG = 0;
        }

        step = 1;
        
        this->reset_hessian();
        std::copy(r.begin(), r.end(), dr.begin());
        for(size_t i = 0; i < vec_size; ++i){
            dr[i] -= dLAG*HG[i];
        }
        double drnorm = std::sqrt(cblas_ddot(dr.size(), dr.data(), 1, dr.data(), 1) + dLAG*dLAG);
        logger::quick_log("gamma GHr GHG", gamma, GHr, GHG);
        logger::quick_log("dLag", dLAG);
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
            const auto get_g = [&](double alpha)->double{
                for(size_t i = 0; i < u.size(); ++i){
                    u1[i] = u[i] + alpha*dr[i];
                }
                double LAG_1 = LAG + alpha*dLAG;
                std::fill(Ku.begin(), Ku.end(), 0);
                std::fill(Kd1.begin(), Kd1.end(), 0);
                this->matrix->set_lag_displ_simple(LAG_1);
                this->matrix->dot_vector(u1, Ku);
                this->matrix->dot_vector(dr, Kd1);

                std::fill(Ku_tmp.begin(), Ku_tmp.end(), 0);
                this->matrix->append_Ku_frictionless_simple(mesh, u1, Ku_tmp);
                this->matrix->append_dKu_frictionless_simple(mesh, u1, dr, alpha, Kd1);
                for(size_t i = 0; i < vec_size; ++i){
                    Ku[i] += Ku_tmp[i];
                    Kd1[i] += dLAG*Ku_tmp[i];
                }

                for(size_t i = 0; i < vec_size; ++i){
                    r[i] = (Ku[i] - f[i]);
                }
                double gamma = 0, dgamma = 0, unused = 0;
                get_g_dg_ddg(gamma, dgamma, unused, alpha, u1, dr);
                if(dLAG != 0){
                    return cblas_ddot(vec_size, Kd1.data(), 1, r.data(), 1) + dgamma*gamma;
                } else {
                    return cblas_ddot(vec_size, Kd1.data(), 1, r.data(), 1);
                }
            };

            /*
            const auto get_a = [&](double alpha, double& gf)->double{
                for(size_t i = 0; i < u.size(); ++i){
                    u1[i] = u[i] + alpha*dr[i];
                }
                double LAG_1 = LAG + alpha*dLAG;
                std::fill(Ku.begin(), Ku.end(), 0);
                std::fill(Kd1.begin(), Kd1.end(), 0);
                std::fill(Kdd1.begin(), Kdd1.end(), 0);
                this->matrix->set_lag_displ_simple(LAG_1);

                std::fill(Ku_tmp.begin(), Ku_tmp.end(), 0);
                this->matrix->append_Ku_frictionless_simple(mesh, u1, Ku_tmp);
                this->matrix->append_dKu_frictionless_simple(mesh, u1, dr, alpha, Kd1);
                this->matrix->append_ddKu_frictionless_simple(mesh, u1, dr, alpha, Kdd1);
                for(size_t i = 0; i < vec_size; ++i){
                    Ku[i] += Ku_tmp[i];
                    Kdd1[i] += 2*dLAG*Kd1[i]/LAG_1;
                    Kd1[i] += dLAG*Ku_tmp[i];
                }
                this->matrix->dot_vector(u1, Ku);
                this->matrix->dot_vector(dr, Kd1);

                for(size_t i = 0; i < vec_size; ++i){
                    r[i] = (Ku[i] - f[i]);
                }
                double gamma = 0, dgamma = 0, ddgamma = 0;
                get_g_dg_ddg(gamma, dgamma, ddgamma, alpha, u1, dr);
                const double g1 = cblas_ddot(vec_size, Kd1.data(), 1, r.data(), 1) + dgamma*gamma;
                const double g2 = cblas_ddot(vec_size, Kd1.data(), 1, Kd1.data(), 1) + dgamma*dgamma +
                                  cblas_ddot(vec_size, Kdd1.data(), 1, r.data(), 1) + ddgamma*gamma;
                gf = g1;
                return alpha - g1/g2;
            };
            */
            double g1 = 0, g2 = 1;
            double a1 = 0, a2 = this->max_step;
            g1 = get_g(a1);
            double g0 = g1;
            g2 = get_g(a2);
            logger::quick_log("g0", g0);
            logger::quick_log("g2", g2);
            logger::quick_log("Improving starting points...");
            while(g0 > 0){
                a1 -= 0.01;
                g1 = get_g(a1);
                g0 = g1;
                if(g0 < 0) break;
                g1 = get_g(std::abs(a1));
                g0 = g1;
                if(g0 < 0){
                    a1 = -a1;
                }
            }
            if(std::abs(g2/g1) > 1e4 && g2*g1 < 0){
                //a2 = std::abs(g1/g2);
                //g2 = get_g(a2);
                a2 = std::abs(g1/g2);
                g2 = get_g(a2);
                while(g2*g1 > 0){
                    a2 *= 2;
                    g2 = get_g(a2);
                    //a2 = std::sqrt(a2);
                    //a2 = std::min(this->max_step, std::max(0.0, a2));
                    //g2 = get_g(a2);
                }
                logger::quick_log(a2, g2);
            }
            const double DIFF = 1e-2;
            const double stop = DIFF*std::min(std::abs(g0), std::abs(g2));
            //if(std::abs(g1) < std::abs(g2)){
            //    std::swap(a1, a2);
            //    std::swap(g1, g2);
            //}
            if(g2*g1 < 0){
                logger::quick_log("Regula falsi...");
                if(g1 > g2){
                    std::swap(g1, g2);
                    std::swap(a1, a2);
                }
                // Regula falsi
                //while(std::abs(g2) > stop){
                //    double diff = g2 - g1;
                //    if(std::abs(diff) < 1e-13){
                //        break;
                //    }
                //    double c = a2 - g2*(a2 - a1)/diff;
                //    if(c > 0){
                //        a2 = c;
                //        g2 = get_g(a2);
                //    } else {
                //        a1 = c;
                //        g1 = get_g(a1);
                //    }
                //}
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
                // Newton 
                // logger::quick_log(a2);
                // while(std::abs(a2 - a1) > 1e-2 || std::abs(g2) > stop){
                //     a1 = a2;
                //     a2 = get_a(a1, g2);
                // }
            } else if(g2 > 0 && g1 > 0){
                a2 = 0;
            }
            logger::quick_log("Applying step...");
            step = this->max_step;
            //logger::quick_log("g2_final", g2);
            if(a2 != this->max_step){
                step = 0.9*a2;
                //step = a2;
            }
            //if(step < 0){
            //    step = 1e-5;
            //    g2 = get_g(step);
            //}
            for(size_t i = 0; i < u.size(); ++i){
                u[i] += step*dr[i];
            }
            LAG += step*dLAG;
            this->matrix->set_lag_displ_simple(LAG);
            std::fill(Ku.begin(), Ku.end(), 0);
            this->matrix->dot_vector(u, Ku);
            this->matrix->append_Ku_frictionless_simple(mesh, u, Ku);
            for(size_t i = 0; i < vec_size; ++i){
                r[i] = -(Ku[i] - f[i]);
            }

            // for(size_t i = 0; i < vec_size; ++i){
            //     r[i] *= -1;
            // }
            // LAG += step*dLAG;
            // std::copy(u1.begin(), u1.end(), u.begin());
            mesh->extend_vector(0, u, u_ext);
            double dg1 = 0;
            get_g_dg_ddg(g1, dg1, dg1, step, u, dr);
            gamma = g1;

        rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1) + gamma*gamma);
       
        std::copy(u.begin() + u_size, u.end(), lambda.begin());

        this->matrix->set_lag_displ_simple(LAG);
        this->matrix->add_frictionless_simple(mesh, mesh->node_positions[0], u_ext, lambda);
        E = 0;
        // This is definitely wrong, but works as an estimate
        for(size_t i = 0; i < vec_size; ++i){
            E += u[i]*(Ku[i]/2 - f[i]);
        }

        logger::quick_log("Iteration:", it);
        logger::quick_log("Potential energy:", E);
        logger::quick_log("L", LAG);
        logger::quick_log("||r||", rnorm);
        const double Kd1norm = std::sqrt(cblas_ddot(Kd1.size(), Kd1.data(), 1, Kd1.data(), 1));
        logger::quick_log("||Kd1||:", Kd1norm);
        logger::quick_log("Step:", step);
        logger::quick_log("");
        const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
        const size_t dof = mesh->elem_info->get_dof_per_node();
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
        logger::quick_log("max_gp", max_gp);
        logger::quick_log("min_gp", min_gp);

        logger::log_assert(rnorm < 1e30, logger::ERROR, "Newton process diverged, increase rtol_abs and/or decrease mult");
        ++it;
        //if(b > 1e-6 || step <= 0){
        //    break;
        //}
        if(step == 0){
            break;
        }
    } while((rnorm > this->rtol_abs || it < 1));

    std::copy(u.begin(), u.begin() + u_size, load.begin());
}
void FiniteElement::solve_frictionless_penalty(const Meshing* const mesh, std::vector<double>& load, const bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u0){
    const size_t vec_size = this->u_size;

    double E = 0;
    size_t it = 0;

    std::vector<double> f(vec_size, 0);
    std::vector<double> u(vec_size, 0);
    std::vector<double> r(vec_size);

    mesh->de_extend_vector(0, u0, u);

    std::vector<double> Ku(vec_size, 0);
    std::vector<double> Kd(vec_size, 0);
    std::vector<double> u_ext(mesh->global_load_vector.size(), 0);
    std::vector<double> f_ext(mesh->global_load_vector.size(), 0);
    std::vector<double> Klag(l_num, 0);

    std::copy(load.begin(), load.end(), f.begin());

    mesh->extend_vector(0, load, f_ext);
    mesh->extend_vector(0, u, u_ext);
    double rnorm = 0;

    double step = 1.0;
    std::fill(Ku.begin(), Ku.end(), 0);
    this->matrix->dot_vector(u, Ku);
    for(size_t i = 0; i < u_size; ++i){
        r[i] = -(Ku[i] - f[i]);
    }
    rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));

    do{

        // Solve linear subproblem
        this->solve(r);

        step = 1;

        for(size_t i = 0; i < u.size(); ++i){
            u[i] += step*r[i];
        }
        mesh->extend_vector(0, u, u_ext);
        this->generate_matrix(mesh, u.size(), 0, mesh->node_positions[0], topopt, D_cache, u_ext, std::vector<double>());
        std::fill(Ku.begin(), Ku.end(), 0);
        this->matrix->dot_vector(u, Ku);

        for(size_t i = 0; i < u_size; ++i){
            r[i] = -(Ku[i] - f[i]);
        }
        E = 0;
        for(size_t i = 0; i < u_size; ++i){
            E += u[i]*(Ku[i]/2 - f[i]);
        }

        logger::quick_log("Iteration:", it);
        logger::quick_log("Potential energy:", E);
        rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
        logger::quick_log("||r||:", rnorm);
        logger::quick_log("Step:", step);
        logger::quick_log("");
        const size_t bnum = mesh->elem_info->get_boundary_nodes_per_element();
        const size_t dof = mesh->elem_info->get_dof_per_node();
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
        logger::quick_log("max_gp", max_gp);
        logger::quick_log("min_gp", min_gp);

        logger::log_assert(rnorm < 1e30, logger::ERROR, "Newton process diverged, increase rtol_abs and/or decrease mult");
        ++it;
    } while((rnorm > this->rtol_abs || it < 1) && step > 0);

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
void FiniteElement::apply_constr_force(const Meshing* mesh, const std::vector<double>& u_ext, const std::vector<double>& lambda, std::vector<double>& b, const size_t dof) const{
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    const size_t num_bound_nodes = mesh->elem_info->get_boundary_nodes_per_element();
    //const size_t max_dof = mesh->elem_info->get_dof_per_node();
    //const size_t kw = mesh->elem_info->get_k_dimension();
    const size_t l_num = mesh->lag_node_map.size();

    math::Vector ln_e(num_bound_nodes);
    math::Vector lp1_e(num_bound_nodes);
    math::Vector lp2_e(num_bound_nodes);
    math::Vector u_e(num_nodes);
    std::vector<long> u_pos(num_nodes);
    for(const auto& e:mesh->paired_boundary){
        for(size_t i = 0; i < num_bound_nodes; ++i){
            const long p1 = mesh->lag_node_map.at(e.b1->nodes[i]->id);
            const long p2 = p1 + l_num;
            const long p3 = p2 + l_num;
            ln_e[i]  = lambda[p1];
            lp1_e[i] = lambda[p2];
            lp2_e[i] = lambda[p3];
        }
        for(size_t i = 0; i < num_nodes; ++i){
            const auto n1 = e.elem->e1->nodes[i];
            const auto n2 = e.elem->e2->nodes[i];
            const long id = this->constr_id_map.at(n2);
            u_pos[i] = id;
            if(id >= 0){
                u_e[i] = u_ext[n1->u_pos[dof]];
            } else {
                u_e[i] = 0;
            }
        }
        const auto be(e.elem->fl3_eq(ln_e, lp1_e, lp2_e, u_e, dof));
        for(size_t i = 0; i < num_nodes; ++i){
            if(u_pos[i] >= 0){
                b[u_pos[i]] += be[i];
            }
        }
    }
}
void FiniteElement::apply_lambda_old(const Meshing* mesh, const std::vector<double>& lambda, std::vector<double>& u_ext){
    const size_t dof = mesh->elem_info->get_dof_per_node();

    /*
    const size_t num_bound_nodes = mesh->elem_info->get_boundary_nodes_per_element();
    std::vector<double> D(this->constr_force.size(), 0);
    for(const auto& e:mesh->paired_boundary){
        const double A = e.elem->get_area()/num_bound_nodes;
        for(size_t i = 0; i < num_bound_nodes; ++i){
            const auto id = this->constr_id_map.at(e.elem->nodes[i]);
            D[id] += A;
        }
    }
    math::Vector l_e(num_bound_nodes);
    const size_t l_num_orig = this->l_num/3;
    for(size_t k = 0; k < dof; ++k){
        std::fill(this->constr_force.begin(), this->constr_force.end(), 0);
        for(const auto& e:mesh->paired_boundary){
            const double A = e.elem->get_area()/num_bound_nodes;
            for(size_t i = 0; i < num_bound_nodes; ++i){
                const long p1 = mesh->lag_node_map.at(e.b1->nodes[i]->id);
                const long p2 = p1 + l_num_orig;
                const long p3 = p2 + l_num_orig;
                l_e[0] = lambda[p3];
                l_e[1] = lambda[p2];
                l_e[2] = lambda[p1]*lambda[p1];
                const auto id = this->constr_id_map.at(e.elem->nodes[i]);
                const auto& R = e.elem->get_R();
                const double l = l_e[0]*R(k, 0) + l_e[1]*R(k, 1) + l_e[2]*R(k, 2);
                this->constr_force[id] += A*l;
                if(k == 2 && id == 118){
                    logger::quick_log(118, A, l);
                    logger::quick_log(R);
                }
            }
        }
        for(const auto n:this->boundary_reference_nodes){
            const long id = this->constr_id_map.at(n);
            if(id == 118){
                logger::quick_log(this->constr_force[id], u_ext[n->u_pos[k]]);
            }
            u_ext[n->u_pos[k]] += this->constr_force[id]/D[id];
            if(id == 118){
                logger::quick_log(this->constr_force[id], u_ext[n->u_pos[k]]);
            }
        }
    }
    */

    //for(size_t i = 0; i < dof; ++i){
    //    std::fill(this->constr_force.begin(), this->constr_force.end(), 0);
    //    this->apply_constr_force(mesh, u_ext, lambda, this->constr_force, i);
    //    this->constr_solver.solve(this->constr_force);
    //    for(const auto n:this->boundary_reference_nodes){
    //        const long id = this->constr_id_map.at(n);
    //        u_ext[n->u_pos[i]] += this->constr_force[id];
    //    }
    //}

    std::fill(this->constr_force.begin(), this->constr_force.end(), 0);
    this->apply_constr_force(mesh, u_ext, lambda, this->constr_force);
    this->constr_solver.solve(this->constr_force);
    for(const auto n:this->boundary_reference_nodes){
        const long id = this->constr_id_map.at(n);
        for(size_t j = 0; j < dof; ++j){
            u_ext[n->u_pos[j]] += this->constr_force[id*dof + j];
        }
    }
    
    //const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    //for(const auto& e:mesh->paired_boundary){
    //    for(size_t i = 0; i < num_nodes; ++i){
    //        //const auto n1 = e.elem->e1->nodes[i];
    //        const auto n2 = e.elem->e2->nodes[i];
    //        const long id = this->constr_id_map.at(n2);
    //        if(id >= 0){
    //            for(size_t j = 0; j < dof; ++j){
    //                u_ext[n2->u_pos[j]] = this->constr_force[id*dof + j];
    //            }
    //        }
    //    }
    //}
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
