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
    const size_t k_dim = std::sqrt(mesh->element_list[0]->get_k().size());
    const size_t dof = k_dim / mesh->element_list[0]->nodes.size();

    int W = mesh->load_vector.size();
    int N = k_dim;
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
        std::vector<double> k = e->get_k();
        if(rho < density.end()){
            cblas_dscal(k.size(), std::pow(*rho, pc), k.data(), 1);
            ++rho;
        }
        this->insert_element_matrix(K, k, u_pos, W, N);
    }
    logger::quick_log("Done.");
    logger::quick_log("Calculating displacements...");
    logger::quick_log("W: ",W," N: ", N);

    int info = LAPACKE_dpbsv(LAPACK_ROW_MAJOR, 'L', W, N-1, 1, K.data(), W, U.data(), 1);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating displacements.", info);

    logger::quick_log("Done.");
   
    return U; 
}

void DirectSolver::insert_element_matrix(std::vector<double>& K, const std::vector<double>& k, const std::vector<long>& pos, int w, int& n) const{
    size_t min_i = 0;
    size_t max_i = 0;
    long min = std::numeric_limits<long>::max();
    long max = -1;
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
    if(pos[max_i] - pos[min_i] + 1 > n){
        n = pos[max_i] - pos[min_i] + 1;
        K.resize(w*n);
    }

    size_t W = pos.size();
    
    for(size_t i = 0; i < W; ++i){
        for(size_t j = i; j < W; ++j){
            if(pos[i] > -1 && pos[j] > -1){
                K[utils::to_band(pos[j], pos[i], w)] += k[W*i + j];
            }
        }
    }
}

}