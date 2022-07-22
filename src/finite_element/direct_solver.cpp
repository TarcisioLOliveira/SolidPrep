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


#include "finite_element/direct_solver.hpp"
#include "lapacke.h"
#include "logger.hpp"
#include "utils.hpp"
#include <limits>
#include <cblas.h>
#include "project_data.hpp"

namespace finite_element{

std::vector<double> DirectSolver::calculate_displacements(ProjectData* data, Meshing* mesh, const std::vector<double>& density, double pc, bool use_stored_matrix, const std::vector<double>& virtual_load){

    if(!use_stored_matrix || this->K.size() == 0){
        this->K.clear();

        const size_t k_dim = data->topopt_element->get_k_dimension();
        const size_t dof = data->topopt_element->get_dof_per_node();

        W = mesh->load_vector.size();
        N = k_dim;
        double rho_min = 1e-6;

        size_t ei = 0;
        for(auto& e:mesh->element_list){
            size_t min_i = 0;
            size_t max_i = 0;
            long min = std::numeric_limits<long>::max();
            long max = -1;
            std::vector<long> pos;
            for(auto& n : e->nodes){
                for(size_t i = 0; i < dof; ++i){
                    pos.push_back(n->u_pos[i]);
                }
            }
            for(size_t i = 0; i < pos.size(); ++i){
                if(pos[i] > -1){
                    if(pos[i] < min){
                        min = pos[i];
                        min_i = i;
                    }
                }
                if(pos[i] > max){
                    max = pos[i];
                    max_i = i;
                }
            }
            size_t N_candidate = pos[max_i] - pos[min_i] + 1;
            if(N_candidate > N){
                N = N_candidate;
            }
            ++ei;
        }


        this->K = std::vector<double>(W*N, 0);


        const auto D = data->geometries[0]->material->stiffness_2D();
        const double t = data->thickness;
        if(density.size() == 0){
            logger::quick_log("Generating stiffness matrix...");
            for(auto& e : mesh->element_list){
                std::vector<long> u_pos;
                for(auto& n : e->nodes){
                    for(size_t i = 0; i < dof; ++i){
                        u_pos.push_back(n->u_pos[i]);
                    }
                }
                std::vector<double> k = e->get_k(D, t);
                this->insert_element_matrix(K, k, u_pos, N);
            }
        } else {
            logger::quick_log("Generating stiffness matrix...");
            auto rho = density.begin();
            for(auto& e : mesh->element_list){
                std::vector<long> u_pos;
                for(auto& n : e->nodes){
                    for(size_t i = 0; i < dof; ++i){
                        u_pos.push_back(n->u_pos[i]);
                    }
                }
                if(rho < density.end() && density.size() >= mesh->element_list.size()){
                    std::vector<double> k = e->get_k(D, t);
                    double rho_scal = rho_min + (1-rho_min)*std::pow(*rho, pc);
                    cblas_dscal(k.size(), rho_scal, k.data(), 1);
                    this->insert_element_matrix(K, k, u_pos, N);
                    ++rho;
                } else {
                    std::vector<double> k = e->get_k(D, t);
                    this->insert_element_matrix(K, k, u_pos, N);
                }
            }
        }
        int info = LAPACKE_dpbtrf_work(LAPACK_COL_MAJOR, 'L', W, N-1, K.data(), N);
        logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while factoring stiffness matrix.", info);
    }
    std::vector<double> U;
    if(virtual_load.size() > 0){
        U = virtual_load;
    } else {
        U = mesh->load_vector;
    }
    logger::quick_log("Done.");
    logger::quick_log("Calculating displacements...");
    logger::quick_log("W: ",W," N: ", N);
    logger::quick_log(U.size());

    int info = LAPACKE_dpbtrs_work(LAPACK_COL_MAJOR, 'L', W, N-1, 1, K.data(), N, U.data(), W);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating displacements.", info);
    // int info = LAPACKE_dpbsv_work(LAPACK_COL_MAJOR, 'L', W, N-1, 1, K.data(), N, U.data(), W);
    // logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating displacements.", info);

    logger::quick_log("Done.");
   
    return U; 

}


void DirectSolver::insert_element_matrix(std::vector<double>& K, const std::vector<double>& k, const std::vector<long>& pos, size_t n) const{
    size_t W = pos.size();
    for(size_t i = 0; i < W; ++i){
        for(size_t j = i; j < W; ++j){
            if(pos[i] > -1 && pos[j] > -1){
                K[utils::to_band(pos[j], pos[i], n)] += k[W*i + j];
            }
        }
    }
}

std::vector<double> DirectSolver::calculate_displacements_simple(ProjectData* data, Meshing* mesh, bool use_stored_matrix){
    if(!use_stored_matrix || this->K.size() == 0){
        this->K.clear();
        const size_t k_dim = data->topopt_element->get_k_dimension();
        const size_t dof = data->topopt_element->get_dof_per_node();

        W = mesh->load_vector.size();
        N = k_dim;

        size_t ei = 0;
        for(auto& e:mesh->element_list){
            size_t min_i = 0;
            size_t max_i = 0;
            long min = std::numeric_limits<long>::max();
            long max = -1;
            std::vector<long> pos;
            for(auto& n : e->nodes){
                for(size_t i = 0; i < dof; ++i){
                    pos.push_back(n->u_pos[i]);
                }
            }
            for(size_t i = 0; i < pos.size(); ++i){
                if(pos[i] > -1){
                    if(pos[i] < min){
                        min = pos[i];
                        min_i = i;
                    }
                }
                if(pos[i] > max){
                    max = pos[i];
                    max_i = i;
                }
            }
            size_t N_candidate = pos[max_i] - pos[min_i] + 1;
            if(N_candidate > N){
                N = N_candidate;
            }
            ++ei;
        }

        K = std::vector<double>(W*N, 0);

        logger::quick_log("Generating stiffness matrix...");
        const auto D = data->geometries[0]->material->stiffness_2D();
        const double t = data->thickness;
        for(auto& e : mesh->element_list){
            std::vector<long> u_pos;
            for(auto& n : e->nodes){
                for(size_t i = 0; i < dof; ++i){
                    u_pos.push_back(n->u_pos[i]);
                }
            }
            std::vector<double> k = e->get_k(D, t);
            this->insert_element_matrix(K, k, u_pos, N);
        }

        int info = LAPACKE_dpbtrf_work(LAPACK_COL_MAJOR, 'L', W, N-1, K.data(), N);
        logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while factoring matrix.", info);
    }

    std::vector<double> U;
    U = mesh->load_vector;
    logger::quick_log("Done.");
    logger::quick_log("Calculating displacements...");
    logger::quick_log("W: ",W," N: ", N);

    int info = LAPACKE_dpbtrs_work(LAPACK_COL_MAJOR, 'L', W, N-1, 1, K.data(), N, U.data(), W);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating displacements.", info);
    // int info = LAPACKE_dpbsv_work(LAPACK_COL_MAJOR, 'L', W, N-1, 1, K.data(), N, U.data(), W);
    // logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating displacements.", info);

    logger::quick_log("Done.");
   
    return U; 

}

}
