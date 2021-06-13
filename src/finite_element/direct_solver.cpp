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

namespace finite_element{

std::vector<double> DirectSolver::calculate_displacements(const std::vector<MeshElement*>& mesh, const std::vector<double>& loads) const{
    size_t k_dim = std::sqrt(mesh[0]->get_k().size());

    int W = loads.size();
    int N = k_dim;
    std::vector<double> K(W*N);
    std::vector<double> U(loads);

    for(auto& e : mesh){
        this->insert_element_matrix(K, e->get_k(), e->u_pos, W, N);
    }

    int info = LAPACKE_dpbsv(LAPACK_ROW_MAJOR, 'L', W, N-1, 1, K.data(), W, U.data(), 1);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating displacements.", info);
   
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
