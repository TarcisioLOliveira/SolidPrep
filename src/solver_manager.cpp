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

#include "solver_manager.hpp"
#include "math/matrix.hpp"
#include "project_data.hpp"
#include "logger.hpp"

void SolverManager::generate_matrix(const Meshing* const mesh, const std::vector<double>& density, double pc, double psi){
    this->topopt = density.size() > 0;
    this->update_D_matrices(mesh, density, pc, psi);
    size_t l_num = 0;
    bool needs_lambda = mesh->proj_data->contact_data.contact_type > FiniteElement::ContactType::FRICTIONLESS_PENALTY;
    if(mesh->proj_data->contact_data.contact_type == FiniteElement::ContactType::FRICTIONLESS_DISPL_SIMPLE){
        l_num = mesh->lag_node_map.size();
    } else if(mesh->proj_data->contact_data.contact_type == FiniteElement::ContactType::FRICTIONLESS_DISPL_LOG){
        l_num = mesh->lag_node_map.size();
    }
    if(this->iteration == 0){
        this->lambdas.resize(mesh->sub_problems->size());
    }
    if(this->iteration == 0){
        this->split_u.resize(mesh->sub_problems->size());
        this->split_u[0].resize(mesh->max_dofs, 0);
    }
    for(size_t i = 0; i < mesh->sub_problems->size(); ++i){
        if(this->iteration == 0){
            if(needs_lambda){
                this->lambdas[i].emplace_back(l_num, 0.00);
                std::fill(this->lambdas[i][0].begin(), this->lambdas[i][0].end(), 0);
            } else {
                this->lambdas[i].emplace_back();
            }
        } else {
            if(needs_lambda){
                std::fill(this->split_u[i].begin(), this->split_u[i].end(), 0.0);
                std::fill(this->lambdas[i][0].begin(), this->lambdas[i][0].end(), 0.0);
            }
        }
        auto& n = mesh->node_positions[i];
        auto& l = mesh->load_vector[i];
        this->solvers[i]->generate_matrix(mesh, l.size(), l_num, n, this->topopt, this->D_matrices, this->split_u[i], this->lambdas[i][0]);
    }
}

void SolverManager::calculate_displacements_global(const Meshing* const mesh, const std::vector<std::vector<double>>& load, std::vector<double>& u){
    ++this->iteration;
    u.resize(mesh->max_dofs, 0);
    std::fill(u.begin(), u.end(), 0);
    std::vector<double> l(mesh->load_vector[0].size());

    this->split_u.resize(mesh->sub_problems->size());
    for(size_t i = 0; i < mesh->sub_problems->size(); ++i){
        std::copy(load[i].begin(), load[i].end(), l.begin());
        this->split_u[i].resize(mesh->max_dofs, 0);


        auto& ui = this->split_u[i];
        this->solvers[i]->calculate_displacements(mesh, l, ui, this->lambdas[i][0]);
        mesh->extend_vector(i, l, ui);
        #pragma omp parallel for
        for(size_t j = 0; j < u.size(); ++j){
            u[j] += ui[j];
        }
    }
    this->solve_step = 1;
}

void SolverManager::calculate_displacements_adjoint(const Meshing* const mesh, const std::vector<std::vector<double>>& load, std::vector<double>& u){
    logger::log_assert(this->iteration > 0, logger::ERROR, "calculate_displacements_global() must be called once before calculate_displacements_adjoint() each iteration");
    size_t l_num = 0;
    bool needs_lambda = mesh->proj_data->contact_data.contact_type > FiniteElement::ContactType::FRICTIONLESS_PENALTY;
    if(mesh->proj_data->contact_data.contact_type == FiniteElement::ContactType::FRICTIONLESS_DISPL_SIMPLE){
        l_num = mesh->lag_node_map.size();
    } else if(mesh->proj_data->contact_data.contact_type == FiniteElement::ContactType::FRICTIONLESS_DISPL_LOG){
        l_num = mesh->lag_node_map.size();
    }
    std::vector<double> u_tmp(mesh->max_dofs, 0);
    for(size_t i = 0; i < mesh->sub_problems->size(); ++i){
        if(this->lambdas[i].size() <= this->solve_step){
            if(needs_lambda){
                this->lambdas[i].emplace_back(l_num, 0);
            } else {
                this->lambdas[i].emplace_back();
            }
        }
        std::vector<double> l(mesh->load_vector[i].size());
        auto& ui = this->split_u[i];
        if(load[i].size() == mesh->max_dofs){
            mesh->de_extend_vector(i, load[i], l);
            this->solvers[i]->calculate_adjoint(mesh, l, ui, this->lambdas[i][this->solve_step]);
        } else if(load[i].size() == mesh->dofs_per_subproblem[i]){
            std::copy(load[i].begin(), load[i].end(), l.begin());
            this->solvers[i]->calculate_adjoint(mesh, l, ui, this->lambdas[i][this->solve_step]);
        } else {
            logger::log_assert(false, logger::ERROR, "incorrect load vector size, must be equal to mesh->max_dofs OR mesh->dofs_per_subproblem[i]");
        }
        mesh->extend_vector(i, l, u_tmp);
        #pragma omp parallel for
        for(size_t j = 0; j < u.size(); ++j){
            u[j] += u_tmp[j];
        }
    }
    ++this->solve_step;
}

void SolverManager::update_D_matrices(const Meshing* const mesh, const std::vector<double>& density, double pc, double psi){
    logger::quick_log("Generating constitutive matrices...");
    if(density.size() == 0 && this->D_matrices.size() > 0 && !this->force_update_materials){
        return;
    } else if(density.size() == 0){
        // Cache non-homogeneous matrices only
        // Makes it faster to generate K when using non-homogeneous matrices,
        // especially when using multiple subproblems
        size_t cache_size = 0;
        for(const auto& g:mesh->geometries){
            const auto& mat = g->materials.get_materials()[0];
            if(!mat->is_homogeneous()){
                cache_size += g->mesh.size();
            }
        }
        this->D_matrices.resize(cache_size);
        size_t offset = 0;
        for(const auto& g:mesh->geometries){
            const auto& mat = g->materials.get_materials()[0];
            if(!mat->is_homogeneous()){
                #pragma omp parallel for
                for(size_t i = 0; i < g->mesh.size(); ++i){
                    const auto& e = g->mesh[i];
                    const auto c = e->get_centroid();
                    this->D_matrices[offset+i] = g->materials.get_D(e.get(), c);
                }
                offset += g->mesh.size();
            }
        }
    } else if(this->old_densities.size() == 0){
        // Cache non-homogeneous materials and multimaterial elements 
        const size_t s_size = mesh->elem_info->get_D_dimension();
        size_t cache_size = 0;
        for(const auto& g:mesh->geometries){
            const auto& mat = g->materials.get_materials()[0];
            if(g->do_topopt || !mat->is_homogeneous()){
                cache_size += g->mesh.size();
            }
        }
        this->D_matrices.resize(cache_size, math::Matrix(s_size, s_size));
        size_t D_offset = 0;
        size_t x_offset = 0;
        for(const auto& g:mesh->geometries){
            const size_t num_den = g->number_of_densities_needed();
            const auto& mat = g->materials.get_materials()[0];
            if(g->do_topopt || !mat->is_homogeneous()){
                if(g->do_topopt){
                    if(g->with_void){
                        #pragma omp parallel for
                        for(size_t i = 0; i < g->mesh.size(); ++i){
                            const auto& e = g->mesh[i];
                            const auto c = e->get_centroid();
                            auto xv = density[x_offset+i];
                            auto rho = density.cbegin() + x_offset + i;
                            g->materials.get_D(rho, psi, e.get(), c, this->D_matrices[D_offset+i]);
                            this->D_matrices[D_offset+i] *= std::pow(xv, pc);
                        }
                    } else {
                        #pragma omp parallel for
                        for(size_t i = 0; i < g->mesh.size(); ++i){
                            const auto& e = g->mesh[i];
                            const auto c = e->get_centroid();
                            auto rho = density.cbegin() + x_offset + i;
                            g->materials.get_D(rho, psi, e.get(), c, this->D_matrices[D_offset+i]);
                        }
                    }
                    x_offset += g->mesh.size()*num_den;
                } else {
                    #pragma omp parallel for
                    for(size_t i = 0; i < g->mesh.size(); ++i){
                        const auto& e = g->mesh[i];
                        const auto c = e->get_centroid();
                        this->D_matrices[D_offset+i] = g->materials.get_D(e.get(), c);
                    }
                }
                D_offset += g->mesh.size();
            }
        }
        this->old_densities = density;
    } else {
        // Only update elements whose densities were modified.
        size_t D_offset = 0;
        size_t x_offset = 0;
        for(const auto& g:mesh->geometries){
            const size_t num_den = g->number_of_densities_needed();
            const auto& mat = g->materials.get_materials()[0];
            if(g->do_topopt || !mat->is_homogeneous()){
                if(g->do_topopt){
                    if(g->with_void){
                        #pragma omp parallel for
                        for(size_t i = 0; i < g->mesh.size(); ++i){
                            if(modified_densities(density, x_offset+i, num_den)){
                                const auto& e = g->mesh[i];
                                const auto c = e->get_centroid();
                                auto xv = density[x_offset+i];
                                auto rho = density.cbegin() + x_offset + i;
                                g->materials.get_D(rho, psi, e.get(), c, this->D_matrices[D_offset+i]);
                                this->D_matrices[D_offset+i] *= std::pow(xv, pc);
                                this->update_densities(density, x_offset+i, num_den);
                            }
                        }
                    } else {
                        #pragma omp parallel for
                        for(size_t i = 0; i < g->mesh.size(); ++i){
                            if(modified_densities(density, x_offset+i, num_den)){
                                const auto& e = g->mesh[i];
                                const auto c = e->get_centroid();
                                auto rho = density.cbegin() + x_offset + i;
                                g->materials.get_D(rho, psi, e.get(), c, this->D_matrices[D_offset+i]);
                                this->update_densities(density, x_offset+i, num_den);
                            }
                        }
                    }
                    x_offset += g->mesh.size()*num_den;
                }
                D_offset += g->mesh.size();
            }
        }
    }
    this->force_update_materials = false;
    logger::quick_log("Done.");
}

std::vector<double> SolverManager::calculate_reactions(const Meshing* const mesh, const std::vector<double>& u) const{
    const double t = mesh->thickness;
    const size_t dof       = mesh->elem_info->get_dof_per_node();
    const size_t node_num  = mesh->elem_info->get_nodes_per_element();
    const size_t k_size    = mesh->elem_info->get_k_dimension();
    std::vector<double> f(u.size(),0);

    size_t D_offset = 0;
    std::vector<long> u_pos(k_size);
    std::vector<double> F(dof);
    for(auto& g : mesh->geometries){
        if(!g->materials.get_materials()[0]->is_homogeneous()){
            size_t it = 0;
            for(const auto& e:g->mesh){
                for(size_t i = 0; i < node_num; ++i){
                    const auto& n = e->nodes[i];
                    for(size_t j = 0; j < dof; ++j){
                        u_pos[i*dof + j] = n->u_pos[j];
                    }
                }
                const auto& D = this->D_matrices[D_offset + it];
                const math::Matrix k = e->get_k(D, t);
                ++it;
                for(size_t i = 0; i < k_size; ++i){
                    for(size_t j = 0; j < k_size; ++j){
                        if(u_pos[i] > -1 && u_pos[j] > -1){
                            const double val = k(i, j)*u[u_pos[j]];
                            f[u_pos[i]] += val;
                            F[i % dof] += val;
                        }
                    }                    
                }
            }
            D_offset += g->mesh.size();
        } else {
            const auto D = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
            for(auto& e : g->mesh){
                for(size_t i = 0; i < node_num; ++i){
                    const auto& n = e->nodes[i];
                    for(size_t j = 0; j < dof; ++j){
                        u_pos[i*dof + j] = n->u_pos[j];
                    }
                }
                const math::Matrix k = e->get_k(D, t);
                for(size_t i = 0; i < k_size; ++i){
                    for(size_t j = 0; j < k_size; ++j){
                        if(u_pos[i] > -1 && u_pos[j] > -1){
                            const double val = k(i, j)*u[u_pos[j]];
                            f[u_pos[i]] += val;
                            F[i % dof] += val;
                        }
                    }                    
                }
            }
        }
    }
    //if(mesh->springs->size() > 0){
    //    std::vector<gp_Pnt> points(bnode_num);
    //    for(size_t it = 0; it < mesh->springs->size(); ++it){
    //        for(const auto& b : mesh->springs->at(it).submesh){
    //            const auto c = b->get_centroid(bnode_num);
    //            const auto& K = mesh->springs->at(it).get_K(b->parent, c);
    //            const auto& e = b->parent;
    //            for(size_t i = 0; i < bnode_num; ++i){
    //                points[i] = b->nodes[i]->point;
    //            }
    //            for(size_t i = 0; i < node_num; ++i){
    //                const auto& n = e->nodes[i];
    //                for(size_t j = 0; j < dof; ++j){
    //                    u_pos[i*dof + j] = n->u_pos[j];
    //                }
    //            }
    //            const std::vector<double> R = e->get_R(K, t, points);
    //            for(size_t i = 0; i < k_size; ++i){
    //                for(size_t j = 0; j < k_size; ++j){
    //                    if(u_pos[i] > -1 && u_pos[j] > -1){
    //                        const double val = R[i*k_size + j]*u[u_pos[j]];
    //                        f[u_pos[i]] += val;
    //                        F[i % dof] += val;
    //                    }
    //                }                    
    //            }
    //        }
    //    }
    //}
    logger::quick_log("Force equilibrium:");
    logger::quick_log(F);
    logger::quick_log("");

    return f;
}
