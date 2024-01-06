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
#include "logger.hpp"

void SolverManager::generate_matrix(const Meshing* const mesh, const std::vector<double>& density, double pc, double psi){
    this->update_D_matrices(mesh, density, pc, psi);
    for(size_t i = 0; i < mesh->sub_problems->size(); ++i){
        auto& n = mesh->node_positions[i];
        auto& l = mesh->load_vector[i];
        this->solvers[i]->generate_matrix(mesh, l.size(), n, density.size() > 0, this->D_matrices);
    }
}

void SolverManager::calculate_displacements_global(const Meshing* const mesh, std::vector<std::vector<double>>& load, std::vector<double>& u){
    u.resize(mesh->max_dofs, 0);
    std::fill(u.begin(), u.end(), 0);
    this->split_u.resize(mesh->sub_problems->size());
    for(size_t i = 0; i < mesh->sub_problems->size(); ++i){
        auto& l = load[i];
        auto& n = mesh->node_positions[i];
        this->split_u[i].resize(mesh->max_dofs, 0);
        auto& ui = this->split_u[i];
        this->solvers[i]->calculate_displacements(l);
        #pragma omp parallel for
        for(size_t j = 0; j < n.size(); ++j){
            const long p_new = n[j];
            if(p_new >= 0){
                u[j] += l[p_new];
                ui[j] = l[p_new];
            }
        }
    }
}

void SolverManager::calculate_displacements_adjoint(const Meshing* const mesh, std::vector<std::vector<double>>& load, std::vector<double>& u){
    for(size_t i = 0; i < mesh->sub_problems->size(); ++i){
        auto& l = load[i];
        auto& n = mesh->node_positions[i];
        this->split_u[i].resize(mesh->max_dofs, 0);
        if(l.size() == mesh->max_dofs){
            std::vector<double> load_pos(mesh->dofs_per_subproblem[i]);
            #pragma omp parallel for
            for(size_t j = 0; j < n.size(); ++j){
                const long p_new = n[j];
                if(p_new >= 0){
                    load_pos[p_new] = l[j];
                }
            }
            this->solvers[i]->calculate_displacements(load_pos);
            #pragma omp parallel for
            for(size_t j = 0; j < n.size(); ++j){
                const long p_new = n[j];
                if(p_new >= 0){
                    u[j] += load_pos[p_new];
                    l[p_new] = load_pos[p_new];
                }
            }
        } else if(l.size() == mesh->dofs_per_subproblem[i]){
            this->solvers[i]->calculate_displacements(l);
            #pragma omp parallel for
            for(size_t j = 0; j < n.size(); ++j){
                const long p_new = n[j];
                if(p_new >= 0){
                    u[j] += l[p_new];
                }
            }
        } else {
            logger::log_assert(false, logger::ERROR, "incorrect load vector size, must be equal to mesh->max_dofs OR mesh->dofs_per_subproblem[i]");
        }
    }
}

void SolverManager::update_D_matrices(const Meshing* const mesh, const std::vector<double>& density, double pc, double psi){
    logger::quick_log("Generating constitutive matrices...");
    if(density.size() == 0 && this->D_matrices.size() > 0){
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
        const size_t DN = s_size*s_size;
        size_t cache_size = 0;
        for(const auto& g:mesh->geometries){
            const auto& mat = g->materials.get_materials()[0];
            if(g->do_topopt || !mat->is_homogeneous()){
                cache_size += g->mesh.size();
            }
        }
        this->D_matrices.resize(cache_size);
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
                            this->D_matrices[D_offset+i].resize(DN);
                            auto xv = density[x_offset+i];
                            auto rho = density.cbegin() + x_offset + i;
                            g->materials.get_D(rho, psi, e.get(), c, this->D_matrices[D_offset+i]);
                            cblas_dscal(DN, std::pow(xv, pc), this->D_matrices[D_offset+i].data(), 1);
                        }
                    } else {
                        #pragma omp parallel for
                        for(size_t i = 0; i < g->mesh.size(); ++i){
                            const auto& e = g->mesh[i];
                            const auto c = e->get_centroid();
                            this->D_matrices[D_offset+i].resize(DN);
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
        const size_t s_size = mesh->elem_info->get_D_dimension();
        const size_t DN = s_size*s_size;
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
                                cblas_dscal(DN, std::pow(xv, pc), this->D_matrices[D_offset+i].data(), 1);
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
    logger::quick_log("Done.");
}
