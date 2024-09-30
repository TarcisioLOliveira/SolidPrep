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
#include "finite_element.hpp"
#include "logger.hpp"
#include "global_stiffness_matrix.hpp"

FiniteElement::FiniteElement(ContactType contact_type, double rtol_abs, GlobalStiffnessMatrix* m)
    :contact_type(contact_type), rtol_abs(rtol_abs), matrix(m){

}

void FiniteElement::generate_matrix(const Meshing* const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache, const std::vector<double>& u_ext){
    this->u_size = u_size;
    this->l_num = l_num;
    this->generate_matrix_base(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, this->contact_type);
}

void FiniteElement::calculate_displacements(const Meshing* const mesh, std::vector<double>& load, const std::vector<double>& u0, std::vector<double>& lambda, const bool topopt, const std::vector<std::vector<double>>& D_cache){
    switch(this->contact_type){
        case RIGID:
            this->solve_rigid(load);
            break;
        case FRICTIONLESS_PENALTY:
            this->solve_frictionless_penalty(mesh, load, topopt, D_cache, u0);
            break;
        case FRICTIONLESS_DISPL:
            this->solve_frictionless_displ(mesh, load, lambda);
            break;
    }
}

void FiniteElement::solve_rigid(std::vector<double>& load){
    this->solve(load);
}
void FiniteElement::solve_frictionless_displ(const Meshing* const mesh, std::vector<double>& load, std::vector<double>& lambda){
    const size_t vec_size = this->u_size + this->l_num;

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
    this->matrix->dot_vector(u, Ku);
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
            for(size_t i = 0; i < u.size(); ++i){
                u_test[i] += step*dr[i];
            }
            mesh->extend_vector(0, u_test, u_ext);

            std::fill(Ku.begin(), Ku.end(), 0);
            this->matrix->dot_vector(u_test, Ku);
            std::fill(Ku.begin() + u_size, Ku.end(), 0);
            this->matrix->append_Ku_frictionless(mesh, u_test, Ku);

            for(size_t i = 0; i < vec_size; ++i){
                r[i] = -(Ku[i] - f[i]);
            }
            rnorm = std::sqrt(cblas_ddot(r.size(), r.data(), 1, r.data(), 1));
            //if(it == 1){
            //    std::fill(Kd.begin(), Kd.end(), 0);
            //    this->matrix->dot_vector(dr, Kd);
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
        //this->generate_matrix(mesh, u.size(), 0, mesh->node_positions[0], topopt, D_cache, u_ext);
        std::fill(Ku.begin(), Ku.end(), 0);
        this->matrix->dot_vector(u, Ku);
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
    } while(rnorm > this->rtol_abs || it < 1);

    std::copy(u.begin(), u.begin() + u_size, load.begin());
}
void FiniteElement::solve_frictionless_penalty(const Meshing* const mesh, std::vector<double>& load, const bool topopt, const std::vector<std::vector<double>>& D_cache, const std::vector<double>& u0){
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
        this->generate_matrix(mesh, u.size(), 0, mesh->node_positions[0], topopt, D_cache, u_ext);
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
    } while(rnorm > this->rtol_abs || it < 1);

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

