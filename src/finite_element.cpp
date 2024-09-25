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
#include <numeric>
#include <set>
#include <cblas.h>
#include "finite_element.hpp"
#include "logger.hpp"
#include "project_data.hpp"
#include "global_stiffness_matrix.hpp"

FiniteElement::FiniteElement(NonlinearSolver* nl, GlobalStiffnessMatrix* m)
    :prob_type(nl->get_class()), nl_type(nl->get_type()), matrix(m), nl_solver(nl){

}

void FiniteElement::generate_matrix(const Meshing* const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache, const std::vector<double>& u_ext){
    MatrixType mtype = MatrixType::RIGID;
    this->u_size = u_size;
    this->l_num = l_num;
    if(l_num != 0){
        mtype = MatrixType::RIGID;
    } else {
        mtype = MatrixType::FRICTIONLESS;
    }
    this->generate_matrix_base(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, mtype);
}

void FiniteElement::calculate_displacements(const Meshing* const mesh, std::vector<double>& load, const std::vector<double>& u0, std::vector<double>& lambda, const bool topopt, const std::vector<std::vector<double>>& D_cache){
    if(this->l_num != 0){
        this->solve_rigid(load);
    } else {
        this->solve_newton(mesh, load, topopt, D_cache, u0);
        //this->solve_frictionless(mesh, load, lambda);
    }
    /*
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
                    u1 = load[ni1];
                }
                if(ni2 > -1){
                    u2 = load[ni2];
                }
                gp += (u1 - u2)*normal.Coord(1+j);
            }
            max_gp = std::max(gp, max_gp);
            min_gp = std::min(gp, min_gp);
        }
    }
    logger::quick_log("max_gp", max_gp);
    logger::quick_log("min_gp", min_gp);
    */
    //logger::quick_log(lambda);
}

void FiniteElement::solve_rigid(std::vector<double>& load){
    this->solve(load);
}
void FiniteElement::solve_frictionless(const Meshing* const mesh, std::vector<double>& load, std::vector<double>& lambda){
    const size_t vec_size = this->u_size + this->l_num;

    this->nl_solver->setup(vec_size);

    double E = 0;
    size_t it = 0;

    std::vector<double> f(vec_size, 0);
    std::vector<double> u(vec_size, 0);
    std::vector<double> r(vec_size);
    std::vector<double> dr(vec_size);
    std::vector<double> u_test(vec_size);

    //std::fill(u.begin(), u.begin() + u_size, 0.0);
    //std::copy(lambda.begin(), lambda.begin() + 2*l_num, u.begin() + u_size);
    std::vector<double> Ku(vec_size, 0);
    std::vector<double> Kd(vec_size, 0);
    std::vector<double> u_ext(mesh->global_load_vector.size(), 0);
    std::vector<double> f_ext(mesh->global_load_vector.size(), 0);
    std::vector<double> Klag(l_num, 0);

    std::copy(load.begin(), load.end(), f.begin());

    mesh->extend_vector(0, load, f_ext);
    mesh->extend_vector(0, u, u_ext);
    //std::copy(load.begin(), load.end(), u.begin());
    //this->solve(u);
    //this->reset_hessian();
    std::copy(lambda.begin(), lambda.end(), u.begin() + u_size);
    //E = -0.5*cblas_ddot(u_ext.size(), f_ext.data(), 1, u_ext.data(), 1); 
    double rnorm = 0;
    double old_rnorm = 0;

    //logger::quick_log("Iteration:", it);
    //logger::quick_log("Potential energy:", E);
    //logger::quick_log("");

    //std::fill(Ku.begin() + u_size + 2*l_num, Ku.end(), 0);
    //double old_step = 1;
    double step = 1.0;
    //std::fill(Ku.begin(), Ku.end(), 0);
    this->dot_vector(u, Ku);
    std::fill(Ku.begin() + u_size, Ku.end(), 0);
    this->matrix->append_Ku_frictionless(mesh, u, Ku);
    for(size_t i = 0; i < vec_size; ++i){
        r[i] = -(Ku[i] - f[i]);
    }
    this->matrix->add_frictionless_part2(mesh, mesh->node_positions[0], u_ext, lambda);
    rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));


    //++it;
    do{

        //this->reset_hessian();
        //std::copy(u.begin() + u_size, u.begin() + u_size + 2*l_num, lambda.begin());


        // Solve linear subproblem
        logger::quick_log("r min max", *std::min_element(r.begin(), r.end()), *std::max_element(r.begin(), r.end()));
        logger::quick_log("u min max", *std::min_element(u.begin(), u.end()), *std::max_element(u.begin(), u.end()));
        this->solve(r);

        step = 1;//this->get_newton_step(r, lambda, Ku);
        //this->nl_solver->update(u.data(), step, r.data());
        
        this->reset_hessian();
        std::copy(r.begin(), r.end(), dr.begin());
        double drnorm = std::sqrt(cblas_ddot(dr.size(), dr.data(), 1, dr.data(), 1));
        logger::quick_log("drnorm", drnorm);
        logger::quick_log("dr min max", *std::min_element(dr.begin(), dr.end()), *std::max_element(dr.begin(), dr.end()));
        old_rnorm = rnorm;
        bool should_break = false;
        do {
            std::copy(u.begin(), u.end(), u_test.begin());
            //if(it != 1){
            //    this->nl_solver->update(u_test.data(), step, dr.data());
            //}
            this->nl_solver->update(u_test.data(), step, dr.data());
            mesh->extend_vector(0, u_test, u_ext);

            std::fill(Ku.begin(), Ku.end(), 0);
            this->dot_vector(u_test, Ku);
            std::fill(Ku.begin() + u_size, Ku.end(), 0);
            this->matrix->append_Ku_frictionless(mesh, u_test, Ku);

            for(size_t i = 0; i < vec_size; ++i){
                r[i] = -(Ku[i] - f[i]);
            }
            rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
            //if(it == 1){
            //    std::fill(Kd.begin(), Kd.end(), 0);
            //    this->dot_vector(dr, Kd);
            //    double dmr = cblas_ddot(r.size(), r.data(), 1, Kd.data(), 1);
            //    double dmrnorm = -dmr/rnorm;
            //    logger::quick_log("dmr", dmr);
            //    logger::quick_log("rnorm", rnorm);
            //    logger::quick_log("dmrnorm", dmrnorm);
            //    exit(0);
            //}
            if(rnorm < old_rnorm){
                should_break = true;
            } else {
                step *= 0.9;
                //old_rnorm = rnorm;
            }
            if(step < 1e-10){
                break;
            }
        } while(!should_break);
        std::copy(u_test.begin(), u_test.end(), u.begin());
        
        //for(size_t i = u_size; i < r.size(); ++i){
        //    std::cout << r[i] << " ";
        //}
        //std::cout << std::endl;
        //std::cout << std::endl;

        //this->nl_solver->update(u.data(), step, r.data());
        std::copy(u.begin() + u_size, u.end(), lambda.begin());

        logger::quick_log(lambda);
        mesh->extend_vector(0, u, u_ext);
        this->reset_hessian();
        //this->generate_matrix(mesh, u.size(), 0, mesh->node_positions[0], topopt, D_cache, u_ext);
        std::fill(Ku.begin(), Ku.end(), 0);
        this->dot_vector(u, Ku);
        std::fill(Ku.begin() + u_size, Ku.end(), 0);
        this->matrix->append_Ku_frictionless(mesh, u, Ku);
        for(size_t i = 0; i < vec_size; ++i){
            r[i] = -(Ku[i] - f[i]);
        }
        this->matrix->add_frictionless_part2(mesh, mesh->node_positions[0], u_ext, lambda);
        //rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
        E = 0;
        for(size_t i = 0; i < vec_size; ++i){
            E += u[i]*(Ku[i]/2 - f[i]);
        }
        //for(size_t i = 0; i < u_size; ++i){
        //    r[i] = -(Ku[i] - f[i]);
        //}
        //rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));

        //this->reset_hessian();

        logger::quick_log("Iteration:", it);
        logger::quick_log("Potential energy:", E);
        logger::quick_log("||r|| expected:", rnorm);
        rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
        logger::quick_log("||r|| actual:", rnorm);
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
                    gp += (u1 - u2)*normal.Coord(1+j);
                }
                max_gp = std::max(gp, max_gp);
                min_gp = std::min(gp, min_gp);
            }
        }
        logger::quick_log("max_gp", max_gp);
        logger::quick_log("min_gp", min_gp);
        //if(it == 10) exit(0);

        logger::log_assert(rnorm < 1e30, logger::ERROR, "Newton process diverged, increase rtol_abs and/or decrease mult");
        ++it;
    } while(rnorm > this->nl_solver->rtol_abs || it < 1);// && step > 1e-6);

    std::copy(u.begin(), u.begin() + u_size, load.begin());
}
void FiniteElement::solve_opt(const Meshing* const mesh, std::vector<double>& load, std::vector<double>& lambda, const bool topopt, const std::vector<std::vector<double>>& D_cache){
    const size_t l_num = mesh->lambda_elements.size();
    const size_t u_size = load.size();

    this->nl_solver->setup(l_num);
    std::vector<double> grad(l_num, 0);
    double E = 0;
    size_t it = 0;

    std::vector<double> f(load.size() + 2*l_num);
    std::vector<double> u_ext(mesh->global_load_vector.size(), 0);
    std::vector<double> f_ext(mesh->global_load_vector.size(), 0);
    mesh->extend_vector(0, load, f_ext);

    do{
        // Solve linear subproblem
        std::copy(load.begin(), load.end(), f.begin());
        std::fill(f.begin() + load.size(), f.end(), 0);
        this->apply_lambda_force(mesh, f, lambda, topopt, D_cache);

        this->solve(f);
        std::copy(f.begin() + u_size, f.end(), lambda.begin());
        mesh->extend_vector(0, f, u_ext);
        //this->calculate_gradient2(mesh, grad, u_ext, f_ext, lambda, topopt, D_cache);
        mesh->apply_lambda(lambda, u_ext);

        E = -0.5*cblas_ddot(u_ext.size(), f_ext.data(), 1, u_ext.data(), 1); 

        this->calculate_gradient(mesh, grad, u_ext, f_ext, lambda, topopt, D_cache);
        double gnorm = cblas_ddot(grad.size(), grad.data(), 1, grad.data(), 1);

        logger::quick_log("Iteration:", it);
        logger::quick_log("Potential energy:", E);
        logger::quick_log("grad norm:", gnorm);
        logger::quick_log("");

        ++it;
        //if(it == 50){
        //    exit(0);
        //}

        // Update
    } while(!this->nl_solver->update(lambda.data() + 2*l_num, E, grad.data()));

    std::copy(load.begin(), load.end(), f.begin());
    std::fill(f.begin() + load.size(), f.end(), 0);
    this->apply_lambda_force(mesh, f, lambda, topopt, D_cache);
    this->solve(f);
    std::copy(f.begin(), f.begin() + u_size, load.begin());
    std::copy(f.begin() + u_size, f.end(), lambda.begin());
        mesh->extend_vector(0, f, u_ext);
        //this->calculate_gradient2(mesh, grad, u_ext, f_ext, lambda, topopt, D_cache);
        mesh->apply_lambda(lambda, u_ext);

        E = -0.5*cblas_ddot(u_ext.size(), f_ext.data(), 1, u_ext.data(), 1); 

        this->calculate_gradient(mesh, grad, u_ext, f_ext, lambda, topopt, D_cache);
        double gnorm = cblas_ddot(grad.size(), grad.data(), 1, grad.data(), 1);
        logger::quick_log("Iteration:", it);
        logger::quick_log("Potential energy:", E);
        logger::quick_log("grad norm:", gnorm);
        logger::quick_log("");
}
void FiniteElement::solve_newton(const Meshing* const mesh, std::vector<double>& load, std::vector<double>& lambda, const bool topopt, const std::vector<std::vector<double>>& D_cache){
    const size_t l_num = mesh->lambda_elements.size();
    const size_t u_size = load.size();
    const size_t vec_size = u_size + 3*l_num;

    this->nl_solver->setup(vec_size);

    double E = 0;
    size_t it = 0;

    std::vector<double> f(vec_size, 0);
    std::vector<double> u(vec_size, 0);
    std::vector<double> r(vec_size);
    std::vector<double> dr(vec_size);
    std::vector<double> u_test(vec_size);
    std::vector<double> l0(lambda.size(), 0);

    std::copy(load.begin(), load.end(), f.begin());
    std::copy(load.begin(), load.end(), u.begin());

    //std::fill(u.begin(), u.begin() + u_size, 0.0);
    //std::copy(lambda.begin(), lambda.begin() + 2*l_num, u.begin() + u_size);
    std::vector<double> Ku(vec_size, 0);
    std::vector<double> Ku0(vec_size, 0);
    std::vector<double> u_ext(mesh->global_load_vector.size(), 0);
    std::vector<double> f_ext(mesh->global_load_vector.size(), 0);
    std::vector<double> Klag(l_num, 0);

    std::copy(load.begin(), load.end(), u.begin());

    this->apply_lambda_force(mesh, u, lambda, topopt, D_cache);

    //std::fill(Ku.begin() + u_size + 2*l_num, Ku.end(), 0.5*1e-9);
    this->generate_hessian(l0, Ku);
    this->solve(u);
    bool Ku_neg = false;
    do{
        Ku_neg = false;
        std::copy(u.begin() + u_size, u.begin() + u_size + 2*l_num, lambda.begin());
        for(size_t i = 0; i < l_num; ++i){
            const double l = lambda[2*l_num + i];
            u[u_size + 2*l_num + i] = l*l;
        }
        this->reset_hessian();
        this->dot_vector(u, Ku);
        //for(size_t i = 0; i < l_num; ++i){
        //    const size_t j = u_size + 2*l_num + i;
        //    double& Kuj = Ku[j];
        //    if(Kuj < 0){
        //        Ku_neg = true;
        //        lambda[2*l_num + i] *= 1.2;
        //    }
        //}
        if(Ku_neg){
            std::fill(Ku.begin(), Ku.end(), 0);
            std::fill(u.begin(), u.end(), 0);
            std::copy(load.begin(), load.end(), u.begin());
            this->apply_lambda_force(mesh, u, lambda, topopt, D_cache);
            this->generate_hessian(l0, Ku);
            this->solve(u);
        }
    } while(Ku_neg);

    mesh->extend_vector(0, load, f_ext);
    mesh->extend_vector(0, u, u_ext);
    mesh->apply_lambda(lambda, u_ext);
    E = -0.5*cblas_ddot(u_ext.size(), f_ext.data(), 1, u_ext.data(), 1); 
    double rnorm = 0;
    double old_rnorm = 0;

    logger::quick_log("Iteration:", it);
    logger::quick_log("Potential energy:", E);
    logger::quick_log("");

    //std::fill(Ku.begin() + u_size + 2*l_num, Ku.end(), 0);
    //double old_step = 1;
    double step = 1.0;

    ++it;
    do{
        //bool Ku_neg = false;
        //do{
        //    Ku_neg = false;
        //    for(size_t i = 0; i < l_num; ++i){
        //        const double l = lambda[2*l_num + i];
        //        u[u_size + 2*l_num + i] = l*l;
        //    }
        //    std::fill(Ku.begin(), Ku.end(), 0);
        //    this->dot_vector(u, Ku);
        //    for(size_t i = 0; i < l_num; ++i){
        //        const size_t j = u_size + 2*l_num + i;
        //        double& Kuj = Ku[j];
        //        if(Kuj < 0){
        //            Ku_neg = true;
        //            lambda[2*l_num + i] *= 1.2;
        //        }
        //    }
        //} while(Ku_neg);

        //std::fill(Ku.begin(), Ku.end(), 0);
        //this->generate_hessian(l0, Ku);
        //std::fill(u.begin(), u.end(), 0);
        //std::copy(load.begin(), load.end(), u.begin());
        //this->apply_lambda_force(mesh, u, lambda, topopt, D_cache);
        //this->solve(u);
        //bool Ku_neg = false;
        //do{
        //    Ku_neg = false;
        //    std::copy(u.begin() + u_size, u.begin() + u_size + 2*l_num, lambda.begin());
        //    for(size_t i = 0; i < l_num; ++i){
        //        const double l = lambda[2*l_num + i];
        //        u[u_size + 2*l_num + i] = l*l;
        //    }
        //    this->reset_hessian();
        //    this->dot_vector(u, Ku);
        //    //for(size_t i = 0; i < l_num; ++i){
        //    //    const size_t j = u_size + 2*l_num + i;
        //    //    double& Kuj = Ku[j];
        //    //    if(Kuj < 0){
        //    //        Ku_neg = true;
        //    //        lambda[2*l_num + i] *= 1.2;
        //    //    }
        //    //}
        //    if(Ku_neg){
        //        std::fill(Ku.begin(), Ku.end(), 0);
        //        std::fill(u.begin(), u.end(), 0);
        //        std::copy(load.begin(), load.end(), u.begin());
        //        this->apply_lambda_force(mesh, u, lambda, topopt, D_cache);
        //        this->generate_hessian(l0, Ku);
        //        this->solve(u);
        //    }
        //} while(Ku_neg);

        for(size_t i = 0; i < l_num; ++i){
            const double l = lambda[2*l_num + i];
            u[u_size + 2*l_num + i] = l*l;
        }
        std::fill(Ku.begin(), Ku.end(), 0);
        this->dot_vector(u, Ku);

        //this->reset_hessian();
        //std::copy(u.begin() + u_size, u.begin() + u_size + 2*l_num, lambda.begin());

        this->generate_hessian(lambda, Ku0);
        for(size_t i = 0; i < u_size + 2*l_num; ++i){
            r[i] = -(Ku[i] - f[i]);
        }
        for(size_t i = 0; i < l_num; ++i){
            const size_t j = u_size + 2*l_num + i;
            const size_t li = 2*l_num + i;
            const double l = 0.5;//lambda[li];
            r[j] = -(2*l*(Ku[j] - f[j]));///std::sqrt(l*l+1e-14);
        }
        rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
        // Solve linear subproblem
        this->solve(r);
        std::copy(r.begin(), r.end(), dr.begin());
        logger::quick_log("-r.dr", -cblas_ddot(r.size(), r.data(), 1, dr.data(), 1));
        double max_dri = 0;
        for(auto& dri:dr){
            //dri *= -1;
            const double abs_dri = std::abs(dri);
            if(abs_dri > max_dri){
                max_dri = abs_dri;
            }
        }
        for(auto& dri:dr){
            dri /= max_dri;
        }
        //exit(0);

        std::copy(lambda.begin() + 2*l_num, lambda.end(), u.begin() + u_size + 2*l_num);
        step = 1;//this->get_newton_step(r, lambda, Ku);
        //if(step > old_step){
        //    step = 0.1*old_step;
        //}
        //old_step = step;
        old_rnorm = rnorm;
        this->reset_hessian();
        bool should_break = false;
        do {
            std::copy(u.begin(), u.end(), u_test.begin());
            this->nl_solver->update(u_test.data(), step, dr.data());
            std::copy(u_test.begin() + u_size, u_test.end(), lambda.begin());
            for(size_t i = 0; i < l_num; ++i){
                const double l = lambda[2*l_num + i];
                u_test[u_size + 2*l_num + i] = l*l;//std::sqrt(l*l+1e-14) - 1e-7;
            }
            std::fill(Ku.begin(), Ku.end(), 0);
            this->dot_vector(u_test, Ku);

            for(size_t i = 0; i < u_size + 2*l_num; ++i){
                r[i] = -(Ku[i] - f[i]);
            }
            for(size_t i = 0; i < l_num; ++i){
                const size_t j = u_size + 2*l_num + i;
                const size_t li = 2*l_num + i;
                const double l = 0.5;//lambda[li];
                r[j] = -(2*l*(Ku[j] - f[j]));///std::sqrt(l*l+1e-14);
                //if(Ku[j] < 0){
                //    Ku_pos = false;
                //}
            }
            rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
            if(should_break){
                break;
            }
            if(rnorm < old_rnorm){
                should_break = true;
            } else {
                step *= 0.9;
                //old_rnorm = rnorm;
            }
            if(step < 1e-10){
                break;
            }
        } while(!should_break);
        //std::copy(u_test.begin(), u_test.end(), u.begin());
        //std::copy(u.begin() + u_size, u.end(), lambda.begin());
        //this->nl_solver->update(u.data(), step, dr.data());
        mesh->extend_vector(0, u, u_ext);
        mesh->apply_lambda(lambda, u_ext);
        E = -0.5*cblas_ddot(u_ext.size(), f_ext.data(), 1, u_ext.data(), 1); 

        //this->reset_hessian();
        std::fill(Ku.begin(), Ku.end(), 0);

        logger::quick_log("Iteration:", it);
        logger::quick_log("Potential energy:", E);
        logger::quick_log("||r||:", rnorm);
        logger::quick_log("Step:", step);
        logger::quick_log("");
        //if(it == 10) exit(0);

        logger::log_assert(rnorm < 1e20, logger::ERROR, "Newton process diverged, increase rtol_abs and/or decrease mult");
        ++it;
    } while(rnorm > this->nl_solver->rtol_abs);// && step > 1e-6);

    std::copy(u.begin(), u.begin() + u_size, load.begin());
}
void FiniteElement::solve_newton(const Meshing* const mesh, std::vector<double>& load, const bool topopt, const std::vector<std::vector<double>>& D_cache, const std::vector<double>& u0){
    const size_t l_num = mesh->lambda_elements.size();
    const size_t u_size = load.size();
    const size_t vec_size = u_size;

    this->nl_solver->setup(vec_size);

    double E = 0;
    size_t it = 0;

    std::vector<double> f(vec_size, 0);
    std::vector<double> u(vec_size, 0);
    std::vector<double> r(vec_size);
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
    mesh->extend_vector(0, u, u_ext);
    //E = -0.5*cblas_ddot(u_ext.size(), f_ext.data(), 1, u_ext.data(), 1); 
    double rnorm = 0;
    double old_rnorm = 0;

    //logger::quick_log("Iteration:", it);
    //logger::quick_log("Potential energy:", E);
    //logger::quick_log("");

    //std::fill(Ku.begin() + u_size + 2*l_num, Ku.end(), 0);
    //double old_step = 1;
    double step = 1.0;
    std::fill(Ku.begin(), Ku.end(), 0);
    this->dot_vector(u, Ku);
    for(size_t i = 0; i < u_size; ++i){
        r[i] = -(Ku[i] - f[i]);
    }
    rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));

    //++it;
    do{

        //this->reset_hessian();
        //std::copy(u.begin() + u_size, u.begin() + u_size + 2*l_num, lambda.begin());

        // Solve linear subproblem
        logger::quick_log("r min max", *std::min_element(r.begin(), r.end()), *std::max_element(r.begin(), r.end()));
        logger::quick_log("u min max", *std::min_element(u.begin(), u.end()), *std::max_element(u.begin(), u.end()));
        this->solve(r);
        //double max_dri = 0;
        //for(auto& dri:dr){
        //    //dri *= -1;
        //    const double abs_dri = std::abs(dri);
        //    if(abs_dri > max_dri){
        //        max_dri = abs_dri;
        //    }
        //}
        //for(auto& dri:dr){
        //    dri /= max_dri;
        //}
        //exit(0);

        step = 1;//this->get_newton_step(r, lambda, Ku);
        //this->nl_solver->update(u.data(), step, r.data());
        std::copy(r.begin(), r.end(), dr.begin());
        double drnorm = std::sqrt(cblas_ddot(dr.size(), dr.data(), 1, dr.data(), 1));
        //for(auto& dri:dr){
        //    if(std::abs(dri) > drnorm*1e-5){
        //        dri /= drnorm*1e2;
        //    }
        //}
        logger::quick_log("drnorm", drnorm);
        logger::quick_log("dr min max", *std::min_element(dr.begin(), dr.end()), *std::max_element(dr.begin(), dr.end()));
        //if(step > old_step){
        //    step = 0.1*old_step;
        //}
        //old_step = step;
        old_rnorm = rnorm;
        bool should_break = false;
        this->generate_matrix(mesh, u.size(), 0, mesh->node_positions[0], topopt, D_cache, u_ext);
        do {
            std::copy(u.begin(), u.end(), u_test.begin());
            //if(it != 1){
            //    this->nl_solver->update(u_test.data(), step, dr.data());
            //}
            this->nl_solver->update(u_test.data(), step, dr.data());
            mesh->extend_vector(0, u_test, u_ext);
            std::fill(Ku.begin(), Ku.end(), 0);
            this->dot_vector(u_test, Ku);

            for(size_t i = 0; i < u_size; ++i){
                r[i] = -(Ku[i] - f[i]);
            }
            rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
            break;
            //if(it == 1){
            //    std::fill(Kd.begin(), Kd.end(), 0);
            //    this->dot_vector(dr, Kd);
            //    double dmr = cblas_ddot(r.size(), r.data(), 1, Kd.data(), 1);
            //    double dmrnorm = -dmr/rnorm;
            //    logger::quick_log("dmr", dmr);
            //    logger::quick_log("rnorm", rnorm);
            //    logger::quick_log("dmrnorm", dmrnorm);
            //    exit(0);
            //}
            if(rnorm < old_rnorm){
                should_break = true;
            } else {
                step *= 0.9;
                //old_rnorm = rnorm;
            }
            if(step < 1e-10){
                break;
            }
        } while(!should_break);
        std::copy(u_test.begin(), u_test.end(), u.begin());
        //std::copy(u.begin() + u_size, u.end(), lambda.begin());
        //this->nl_solver->update(u.data(), step, dr.data());
        mesh->extend_vector(0, u, u_ext);
        this->generate_matrix(mesh, u.size(), 0, mesh->node_positions[0], topopt, D_cache, u_ext);
        std::fill(Ku.begin(), Ku.end(), 0);
        this->dot_vector(u, Ku);
        for(size_t i = 0; i < u_size; ++i){
            r[i] = -(Ku[i] - f[i]);
        }
        //rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
        E = 0;
        for(size_t i = 0; i < u_size; ++i){
            E += u[i]*(Ku[i]/2 - f[i]);
        }
        //for(size_t i = 0; i < u_size; ++i){
        //    r[i] = -(Ku[i] - f[i]);
        //}
        //rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));

        //this->reset_hessian();

        logger::quick_log("Iteration:", it);
        logger::quick_log("Potential energy:", E);
        logger::quick_log("||r|| expected:", rnorm);
        rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
        logger::quick_log("||r|| actual:", rnorm);
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
        //if(it == 10) exit(0);

        logger::log_assert(rnorm < 1e30, logger::ERROR, "Newton process diverged, increase rtol_abs and/or decrease mult");
        ++it;
    } while(rnorm > this->nl_solver->rtol_abs || it < 1);// && step > 1e-6);

    std::copy(u.begin(), u.begin() + u_size, load.begin());
}

void FiniteElement::apply_lambda_force(const Meshing* const mesh, std::vector<double>& f, const std::vector<double>& lambda, bool topopt, const std::vector<std::vector<double>>& D_cache) const{
    const size_t l_num = mesh->lambda_elements.size();
    const size_t kw = mesh->elem_info->get_k_dimension();
    const size_t lw = mesh->elem_info->get_dof_per_node();
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t t = mesh->thickness;
    const size_t u_size = mesh->load_vector[0].size();

    std::vector<size_t> lv_it(node_num);
    std::vector<const std::vector<size_t>*> lv_vec(node_num);

    const std::vector<long>& npos = mesh->node_positions[0];
    const auto& lambda_list = mesh->lambda_elements;
    std::vector<double> k;
    std::vector<double> R(kw);
    std::vector<double> u_e(kw);
    size_t geom = 0;
    const Geometry* g = mesh->geometries[geom];
    bool use_D_cache = (topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous();
    size_t D_offset = 0;
    std::vector<double> D;
    if(!use_D_cache){
        D = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
    }
    for(const auto& l_e:mesh->lambda_affected_elements){
        if(l_e.geom_id != geom){
            if(!use_D_cache){
                D_offset += g->mesh.size();
            }
            geom = l_e.geom_id;
            const Geometry* g = mesh->geometries[geom];
            use_D_cache = (topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous();
            if(!use_D_cache){
                D = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
            }
        }
        if(use_D_cache){
            k = l_e.e->get_k(D_cache[l_e.e->id - D_offset], t);
        } else {
            k = l_e.e->get_k(D, t);
        }
        const size_t ln = l_e.lambdas.size();
        std::vector<double> R(kw, 0);
        std::vector<double> Rp(kw*ln*2, 0);
        std::vector<double> KR(kw, 0);
        std::vector<double> RpKR(ln*2, 0);

        std::fill(lv_it.begin(), lv_it.end(), 0);
        for(size_t i = 0; i < node_num; ++i){
            const auto& vi = mesh->node_lambda_map.find(l_e.e->nodes[i]->id);
            if(vi != mesh->node_lambda_map.end()){
                lv_vec[i] = &vi->second;
            } else {
                lv_vec[i] = nullptr;
            }
        }
        for(size_t l_it = 0; l_it < ln; ++l_it){
            const size_t l_i = l_e.lambdas[l_it];
            const auto& ll_e = lambda_list[l_i];
            const double ld = lambda[2*l_num + l_i];
            const std::vector<double> R0 =
            {
                ll_e.n.X()*(ld*ld),
                ll_e.n.Y()*(ld*ld),
                ll_e.n.Z()*(ld*ld)
            };
            const std::vector<double> Rp0 =
            {
                ll_e.p1.X(), ll_e.p2.X(),
                ll_e.p1.Y(), ll_e.p2.Y(),
                ll_e.p1.Z(), ll_e.p2.Z()
            };

            for(size_t i = 0; i < node_num; ++i){
                if(lv_vec[i] != nullptr && lv_it[i] < lv_vec[i]->size() && lv_vec[i]->at(lv_it[i]) == l_i){
                    for(size_t ii = 0; ii < lw; ++ii){
                        R[i*lw + ii] += R0[ii];
                        for(size_t jj = 0; jj < 2; ++jj){
                            Rp[(i*lw + ii)*2*ln + (l_it*2 + jj)] = -Rp0[ii*2 + jj];
                        }
                    }
                    ++lv_it[i];
                }
            }
        }

        for(size_t i = 0; i < kw; ++i){
            for(size_t j = 0; j < kw; ++j){
                KR[i] += k[i*kw + j]*R[j];
            }
        }
        for(size_t i = 0; i < 2*ln; ++i){
            for(size_t j = 0; j < kw; ++j){
                RpKR[i] += Rp[j*2*ln + i]*KR[j];
            }
        }

        for(size_t i = 0; i < node_num; ++i){
            const auto& n = l_e.e->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                const size_t p = n->id*dof + j;
                const long fpos = npos[p];
                if(fpos > -1){
                    f[fpos] += KR[(i*dof + j)];
                }
            }
        }
        for(size_t l_it2 = 0; l_it2 < ln; ++l_it2){
            const size_t l_i = l_e.lambdas[l_it2];
            f[u_size + l_i] += RpKR[2*l_it2];
            f[u_size + l_i + l_num] += RpKR[(2*l_it2 + 1)];
        }
    }
}

void FiniteElement::calculate_gradient(const Meshing* const mesh, std::vector<double>& grad, const std::vector<double>& u_ext, const std::vector<double>& f_ext, const std::vector<double>& lambda, bool topopt, const std::vector<std::vector<double>>& D_cache) const{
    const size_t l_num = mesh->lambda_elements.size();
    const size_t kw = mesh->elem_info->get_k_dimension();
    const size_t lw = mesh->elem_info->get_dof_per_node();
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t t = mesh->thickness;

    std::vector<size_t> lv_it(node_num);
    std::vector<const std::vector<size_t>*> lv_vec(node_num);

    const auto& lambda_list = mesh->lambda_elements;
    std::vector<double> k;
    std::vector<double> R(kw);
    std::vector<double> u_e(kw);
    std::vector<double> f_e(kw);
    size_t geom = 0;
    const Geometry* g = mesh->geometries[geom];
    bool use_D_cache = (topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous();
    size_t D_offset = 0;
    std::vector<double> D;
    if(!use_D_cache){
        D = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
    }
    std::fill(grad.begin(), grad.end(), 0);
    for(const auto& l_e:mesh->lambda_affected_elements){
        if(l_e.geom_id != geom){
            if(!use_D_cache){
                D_offset += g->mesh.size();
            }
            geom = l_e.geom_id;
            const Geometry* g = mesh->geometries[geom];
            use_D_cache = (topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous();
            if(!use_D_cache){
                D = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
            }
        }
        if(use_D_cache){
            k = l_e.e->get_k(D_cache[l_e.e->id - D_offset], t);
        } else {
            k = l_e.e->get_k(D, t);
        }
        const size_t ln = l_e.lambdas.size();
        std::vector<double> R(kw*ln, 0);
        std::vector<double> KR(kw*ln, 0);

        for(size_t i = 0; i < node_num; ++i){
            const auto& n = l_e.e->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                const size_t p = n->id*dof + j;
                u_e[i*dof + j] = u_ext[p];
                f_e[i*dof + j] = f_ext[p];
            }
        }

        std::fill(lv_it.begin(), lv_it.end(), 0);
        for(size_t i = 0; i < node_num; ++i){
            const auto& vi = mesh->node_lambda_map.find(l_e.e->nodes[i]->id);
            if(vi != mesh->node_lambda_map.end()){
                lv_vec[i] = &vi->second;
            } else {
                lv_vec[i] = nullptr;
            }
        }
        for(size_t l_it = 0; l_it < ln; ++l_it){
            const size_t l_i = l_e.lambdas[l_it];
            const auto& ll_e = lambda_list[l_i];
            const double ld = lambda[2*l_num + l_i];
            const std::vector<double> R0 =
            {
                ll_e.n.X()*(2.0*ld),
                ll_e.n.Y()*(2.0*ld),
                ll_e.n.Z()*(2.0*ld)
            };

            for(size_t i = 0; i < node_num; ++i){
                if(lv_vec[i] != nullptr && lv_it[i] < lv_vec[i]->size() && lv_vec[i]->at(lv_it[i]) == l_i){
                    for(size_t ii = 0; ii < lw; ++ii){
                        // Transpose R0
                        R[(i*lw + ii)*ln + l_it] = R0[ii];
                    }
                    ++lv_it[i];
                }
            }
        }

        for(size_t i = 0; i < kw; ++i){
            for(size_t j = 0; j < ln; ++j){
                for(size_t l = 0; l < kw; ++l){
                    KR[i*ln + j] += k[i*kw + l]*R[l*ln + j];
                }
            }
        }

        for(size_t l_it = 0; l_it < ln; ++l_it){
            const size_t l_i = l_e.lambdas[l_it];
            double sum = 0;
            for(size_t i = 0; i < kw; ++i){
                sum += f_e[i]*R[i*ln + l_it] - u_e[i]*KR[i*ln + l_it];
            }
            grad[l_i] += sum;
        }
    }

    /* Parallelizable version (incomplete)
       (hopefully this is unnecessary)

    for(size_t i = 0; i < lambda_list.size(); ++i){
        const auto& l_e = lambda_list[i];
        const double li = lambda[2*l_num + i];
        if(l_e.parent->geom_id != geom){
            if(!use_D_cache){
                D_offset += g->mesh.size();
            }
            geom = l_e.parent->geom_id;
            const Geometry* g = mesh->geometries[geom];
            use_D_cache = (topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous();
            if(!use_D_cache){
                D = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
            }
        }
        if(use_D_cache){
            k = l_e.parent->parent->get_k(D_cache[l_e.parent->parent->id - D_offset], t);
        } else {
            k = l_e.parent->parent->get_k(D, t);
        }

        std::vector<double> R0 =
        {
            l_e.n.X()*(-2.0*li),
            l_e.p1.X()*(-2.0*li),
            l_e.p2.X()*(-2.0*li)
        };
        std::set<MeshElement*> affected_elements;


    }
    */
}

void FiniteElement::calculate_gradient2(const Meshing* const mesh, std::vector<double>& grad, const std::vector<double>& u_ext_pure, const std::vector<double>& f_ext, const std::vector<double>& lambda, bool topopt, const std::vector<std::vector<double>>& D_cache) const{
    const size_t l_num = mesh->lambda_elements.size();
    const size_t kw = mesh->elem_info->get_k_dimension();
    const size_t lw = mesh->elem_info->get_dof_per_node();
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t t = mesh->thickness;

    std::vector<size_t> lv_it(node_num);
    std::vector<const std::vector<size_t>*> lv_vec(node_num);

    const auto& lambda_list = mesh->lambda_elements;
    std::vector<double> k;
    std::vector<double> R(kw);
    std::vector<double> u_e(kw);
    std::vector<double> f_e(kw);
    size_t geom = 0;
    const Geometry* g = mesh->geometries[geom];
    bool use_D_cache = (topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous();
    size_t D_offset = 0;
    std::vector<double> D;
    if(!use_D_cache){
        D = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
    }
    std::fill(grad.begin(), grad.end(), 0);
    for(const auto& l_e:mesh->lambda_affected_elements){
        if(l_e.geom_id != geom){
            if(!use_D_cache){
                D_offset += g->mesh.size();
            }
            geom = l_e.geom_id;
            const Geometry* g = mesh->geometries[geom];
            use_D_cache = (topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous();
            if(!use_D_cache){
                D = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
            }
        }
        if(use_D_cache){
            k = l_e.e->get_k(D_cache[l_e.e->id - D_offset], t);
        } else {
            k = l_e.e->get_k(D, t);
        }
        const size_t ln = l_e.lambdas.size();
        std::vector<double> R(kw*ln, 0);
        std::vector<double> KR(kw*ln, 0);
        std::vector<double> RKR(ln*ln, 0);
        std::vector<double> ll(ln, 0);

        for(size_t i = 0; i < node_num; ++i){
            const auto& n = l_e.e->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                const size_t p = n->id*dof + j;
                u_e[i*dof + j] = u_ext_pure[p];
                f_e[i*dof + j] = f_ext[p];
            }
        }

        std::fill(lv_it.begin(), lv_it.end(), 0);
        for(size_t i = 0; i < node_num; ++i){
            const auto& vi = mesh->node_lambda_map.find(l_e.e->nodes[i]->id);
            if(vi != mesh->node_lambda_map.end()){
                lv_vec[i] = &vi->second;
            } else {
                lv_vec[i] = nullptr;
            }
        }
        for(size_t l_it = 0; l_it < ln; ++l_it){
            const size_t l_i = l_e.lambdas[l_it];
            const auto& ll_e = lambda_list[l_i];
            const double ld = lambda[2*l_num + l_i];
            const std::vector<double> R0 =
            {
                ll_e.n.X(),
                ll_e.n.Y(),
                ll_e.n.Z(),
            };

            ll[l_it] = ld*ld;

            for(size_t i = 0; i < node_num; ++i){
                if(lv_vec[i] != nullptr && lv_it[i] < lv_vec[i]->size() && lv_vec[i]->at(lv_it[i]) == l_i){
                    for(size_t ii = 0; ii < lw; ++ii){
                        // Transpose R0
                        R[(i*lw + ii)*ln + l_it] = R0[ii];
                    }
                    ++lv_it[i];
                }
            }
        }

        for(size_t i = 0; i < kw; ++i){
            for(size_t j = 0; j < ln; ++j){
                for(size_t l = 0; l < kw; ++l){
                    KR[i*ln + j] += k[i*kw + l]*R[l*ln + j];
                }
            }
        }
        for(size_t i = 0; i < ln; ++i){
            for(size_t j = 0; j < ln; ++j){
                for(size_t l = 0; l < kw; ++l){
                    RKR[i*ln + j] += R[l*ln + i]*KR[l*ln + j];
                }
            }
        }

        for(size_t l_it = 0; l_it < ln; ++l_it){
            const size_t l_i = l_e.lambdas[l_it];
            double sum = 0;
            for(size_t i = 0; i < kw; ++i){
                sum += f_e[i]*R[i*ln + l_it] - u_e[i]*KR[i*ln + l_it];
            }
            for(size_t i = 0; i < ln; ++i){
                sum -= RKR[l_it*ln + i]*ll[i];
            }
            grad[l_i] -= sum;
        }
    }

    /* Parallelizable version (incomplete)
       (hopefully this is unnecessary)

    for(size_t i = 0; i < lambda_list.size(); ++i){
        const auto& l_e = lambda_list[i];
        const double li = lambda[2*l_num + i];
        if(l_e.parent->geom_id != geom){
            if(!use_D_cache){
                D_offset += g->mesh.size();
            }
            geom = l_e.parent->geom_id;
            const Geometry* g = mesh->geometries[geom];
            use_D_cache = (topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous();
            if(!use_D_cache){
                D = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
            }
        }
        if(use_D_cache){
            k = l_e.parent->parent->get_k(D_cache[l_e.parent->parent->id - D_offset], t);
        } else {
            k = l_e.parent->parent->get_k(D, t);
        }

        std::vector<double> R0 =
        {
            l_e.n.X()*(-2.0*li),
            l_e.p1.X()*(-2.0*li),
            l_e.p2.X()*(-2.0*li)
        };
        std::set<MeshElement*> affected_elements;


    }
    */
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

