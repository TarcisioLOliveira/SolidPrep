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
    long first = -1;
    long last = 0;
    for(size_t i = 0; i < pos.size(); ++i){
        if(first == -1 && pos[i] > -1){
            first = i;
        }
        if(pos[i] > -1){
            last = i;
        }
    }
    if(pos[last] - pos[first] > n){
        n = pos[last] - pos[first];
        K.resize(w*n);
    }

    int W = pos.size();
    
    for(long i = first; i < last+1; ++i){
        for(long j = i; j < last+1; ++j){
            if(pos[i] > -1 && pos[j] > -1){
                K[utils::to_band(pos[j], pos[i], w)] += k[W*i + j];
            }
        }
    }
}

}
