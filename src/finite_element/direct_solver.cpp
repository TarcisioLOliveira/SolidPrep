/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   Foobar is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Foobar is distributed in the hope that it will be useful,
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

namespace finite_element{

std::vector<double> DirectSolver::calculate_displacements(ProjectData* data, Meshing* mesh, const std::vector<double>& density, double pc, const std::vector<double>& virtual_load) const{

    if(density.size() == 0){
        return this->calculate_displacements_simple(data, mesh);
    }

    const size_t k_dim = std::sqrt(mesh->element_list[0]->get_k().size());
    const size_t dof = k_dim / mesh->element_list[0]->nodes.size();

    size_t W = mesh->load_vector.size();
    size_t N = k_dim;

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
        if(pos[max_i] - pos[min_i] + 1 > N){
            N = pos[max_i] - pos[min_i] + 1;
        }
        ++ei;
    }


    std::vector<double> K(W*N, 0);
    std::vector<double> U;
    if(virtual_load.size() > 0){
        U = virtual_load;
    } else {
        U = mesh->load_vector;
    }

    logger::quick_log("Generating stiffness matrix...");
    auto rho = density.begin();
    for(auto& e : mesh->element_list){
        std::vector<long> u_pos;
        for(auto& n : e->nodes){
            for(size_t i = 0; i < dof; ++i){
                u_pos.push_back(n->u_pos[i]);
            }
        }
        if(rho < density.end() && density.size() == mesh->element_list.size()){
            std::vector<double> k = e->get_k();
            cblas_dscal(k.size(), std::pow(*rho, pc), k.data(), 1);
            this->insert_element_matrix(K, k, u_pos, W, N);
            ++rho;
        } else {
            std::vector<double> k = e->get_k();
            this->insert_element_matrix(K, k, u_pos, W, N);
        }
    }
    logger::quick_log("Done.");
    logger::quick_log("Calculating displacements...");
    logger::quick_log("W: ",W," N: ", N);

    int info = LAPACKE_dpbsv_work(LAPACK_COL_MAJOR, 'L', W, N-1, 1, K.data(), N, U.data(), W);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating displacements.", info);

    logger::quick_log("Done.");
   
    return U; 

}

// Element removal method, currently doesn't work well
// std::vector<double> DirectSolver::calculate_displacements(ProjectData* data, Meshing* mesh, const std::vector<double>& density, double pc, const std::vector<double>& virtual_load) const{
// 
//     if(density.size() == 0){
//         return this->calculate_displacements_simple(data, mesh);
//     }
// 
//     const size_t k_dim = std::sqrt(mesh->element_list[0]->get_k().size());
//     const size_t dof = k_dim / mesh->element_list[0]->nodes.size();
// 
//     size_t W = mesh->load_vector.size();
//     size_t N = k_dim;
// 
//     // std::vector<bool> used(mesh->node_list.size()*dof, false);
//     std::vector<long> new_pos(mesh->node_list.size()*dof, -1);
//     if(density.size() == mesh->element_list.size()){
//         for(size_t i = 0; i < mesh->element_list.size(); ++i){
//             if(density[i] >= 1e-3){
//                 for(auto& n:mesh->element_list[i]->nodes){
//                     for(size_t j = 0; j < dof; ++j){
//                         if(n->u_pos[j] > -1){
//                             new_pos[n->u_pos[j]] = 0;
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     size_t new_W = 0;
//     for(auto& i:new_pos){
//         if(i == 0){
//             i = new_W;
//             ++new_W;
//         }
//     }
//     W = new_W;
// 
//     size_t ei = 0;
//     for(auto& e:mesh->element_list){
//         size_t min_i = 0;
//         size_t max_i = 0;
//         long min = std::numeric_limits<long>::max();
//         long max = -1;
//         std::vector<long> pos;
//         if(density.size() == mesh->element_list.size() && density[ei] < 1e-3){
//             ++ei;
//             continue;
//         }
//         for(auto& n : e->nodes){
//             for(size_t i = 0; i < dof; ++i){
//                 if(n->u_pos[i] > -1){
//                     pos.push_back(new_pos[n->u_pos[i]]);//n->u_pos[i]);
//                 } else {
//                     pos.push_back(-1);
//                 }
//             }
//         }
//         for(size_t i = 0; i < pos.size(); ++i){
//             if(pos[i] > -1){
//                 if(pos[i] < min){
//                     min = pos[i];
//                     min_i = i;
//                 }
//             }
//             if(pos[i] > max){
//                 max = pos[i];
//                 max_i = i;
//             }
//         }
//         if(pos[max_i] - pos[min_i] + 1 > N){
//             N = pos[max_i] - pos[min_i] + 1;
//         }
//         ++ei;
//     }
// 
// 
//     std::vector<double> K(W*N, 0);
//     std::vector<double> U;
//     if(virtual_load.size() > 0){
//         U = virtual_load;
//     } else {
//         U = mesh->load_vector;
//     }
// 
//     logger::quick_log("Generating stiffness matrix...");
//     auto rho = density.begin();
//     for(auto& e : mesh->element_list){
//         std::vector<long> u_pos;
//         for(auto& n : e->nodes){
//             for(size_t i = 0; i < dof; ++i){
//                 if(n->u_pos[i] > -1){
//                     u_pos.push_back(new_pos[n->u_pos[i]]);//n->u_pos[i]);
//                 } else {
//                     u_pos.push_back(-1);
//                 }
//             }
//         }
//         if(rho < density.end() && density.size() == mesh->element_list.size()){
//             if(*rho >= 1e-3){
//                 std::vector<double> k = e->get_k();
//                 cblas_dscal(k.size(), std::pow(*rho, pc), k.data(), 1);
//                 this->insert_element_matrix(K, k, u_pos, W, N);
//             }
//             ++rho;
//         } else {
//             std::vector<double> k = e->get_k();
//             this->insert_element_matrix(K, k, u_pos, W, N);
//         }
//     }
//     logger::quick_log("Done.");
//     logger::quick_log("Calculating displacements...");
//     logger::quick_log("W: ",W," N: ", N);
// 
//     int info = LAPACKE_dpbsv_work(LAPACK_COL_MAJOR, 'L', W, N-1, 1, K.data(), N, U.data(), W);
//     logger::log_assert(info >= 0, logger::ERROR, "LAPACKE returned {} while calculating displacements.", info);
//     if(info > 0){
//         logger::quick_log("Matrix is singular, penalizing...");
//         std::fill(U.begin(), U.end(), 1e6);
//     }
//     bool zerod = true;
//     for(auto& u:U){
//         if(std::abs(u) > 1e-7){
//             zerod = false;
//             break;
//         }
//     }
//     if(zerod){
//         logger::quick_log("Displacements are zero, penalizing...");
//         std::fill(U.begin(), U.end(), 1e6);
//     }
// 
//     logger::quick_log("Done.");
//    
//     return expand_U(new_pos, U, mesh, dof, mesh->load_vector.size());//U; 
// }

void DirectSolver::insert_element_matrix(std::vector<double>& K, const std::vector<double>& k, const std::vector<long>& pos, size_t w, size_t n) const{
    size_t W = pos.size();
    size_t N = 0;
    for(size_t i = 0; i < W; ++i){
        for(size_t j = i; j < W; ++j){
            if(pos[i] > -1 && pos[j] > -1){
                K[utils::to_band(pos[j], pos[i], n)] += k[W*i + j];
            }
        }
    }
}
std::vector<double> DirectSolver::calculate_displacements_simple(ProjectData* data, Meshing* mesh) const{
    const size_t k_dim = std::sqrt(mesh->element_list[0]->get_k().size());
    const size_t dof = k_dim / mesh->element_list[0]->nodes.size();

    size_t W = mesh->load_vector.size();
    size_t N = k_dim;

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
        if(pos[max_i] - pos[min_i] + 1 > N){
            N = pos[max_i] - pos[min_i] + 1;
        }
        ++ei;
    }


    std::vector<double> K(W*N, 0);
    std::vector<double> U;
    U = mesh->load_vector;

    logger::quick_log("Generating stiffness matrix...");
    for(auto& e : mesh->element_list){
        std::vector<long> u_pos;
        for(auto& n : e->nodes){
            for(size_t i = 0; i < dof; ++i){
                u_pos.push_back(n->u_pos[i]);
            }
        }
        std::vector<double> k = e->get_k();
        this->insert_element_matrix(K, k, u_pos, W, N);
    }
    logger::quick_log("Done.");
    logger::quick_log("Calculating displacements...");
    logger::quick_log("W: ",W," N: ", N);

    int info = LAPACKE_dpbsv_work(LAPACK_COL_MAJOR, 'L', W, N-1, 1, K.data(), N, U.data(), W);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating displacements.", info);

    logger::quick_log("Done.");
   
    return U; 

}

std::vector<double> DirectSolver::expand_U(const std::vector<long>& new_pos, std::vector<double>& U, Meshing* mesh, size_t dof, size_t size) const{
    std::vector<double> u(size, 0);
    // for(size_t i = 0; i < new_pos.size(); ++i){
    //     if(new_pos[i] > -1){
    //         u[new_pos[i]] = U[i];
    //     }
    // }
    for(size_t i = 0; i < mesh->node_list.size(); ++i){
        auto& n = mesh->node_list[i];
        for(size_t j = 0; j < dof; ++j){
            if(n->u_pos[j] > -1){
                if(new_pos[n->u_pos[j]] > -1){
                    u[n->u_pos[j]] = U[new_pos[n->u_pos[j]]];
                }
            }
        }
    }

    return u;
}

}
