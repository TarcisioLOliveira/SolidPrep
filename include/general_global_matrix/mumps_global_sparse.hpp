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

#ifndef MUMPS_GLOBAL_SPARSE_HPP
#define MUMPS_GLOBAL_SPARSE_HPP

#include <vector>
#include "logger.hpp"
#include "utils/coo.hpp"

namespace general_global_matrix{

class MUMPSGlobalSparse{
    public:
    MUMPSGlobalSparse() = default;
    MUMPSGlobalSparse(bool spd, size_t L):
        L(L), spd(spd) {}
    ~MUMPSGlobalSparse() = default;

    inline void reinitialize(bool spd, size_t L){
        this->spd = spd;
        this->L = L;
        this->M.clear();
    }

    inline std::vector<int>& get_rows(){
        return this->M.rows;
    }
    inline std::vector<int>& get_cols(){
        return this->M.cols;
    }
    inline std::vector<double>& get_vals(){
        return this->M.vals;
    }

    inline void print_matrix() const{
        for(size_t i = 0; i < L; ++i){
            for(size_t j = 0; j < L; ++j){
                std::cout << this->M.get(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }

    inline void add_element(const std::vector<double>& matrix, const std::vector<long>& pos){
        if(spd){
            this->M.insert_matrix_symmetric(matrix, pos);
        } else {
            this->M.insert_matrix_general(matrix, pos);
        }
    }
    inline void add_element(const std::vector<double>& matrix, const std::vector<long>& pos_i, const std::vector<long>& pos_j){
        this->M.insert_block(matrix, pos_i, pos_j);
    }
    inline void add_value(size_t i, size_t j, double val){
        if(spd){
            if(j > i){
                std::swap(i, j);
            }
            this->M.add(i, j, val);
        } else {
            this->M.add(i, j, val);
        }
    }

    inline void make_zero(){
        this->M.zero();
    }

    inline void clear_matrix(){
        this->M.clear();
        setted = false;
    }

    inline void set_up(){
        if(!setted){
            this->M.generate_coo(L);
            setted = true;
        }
    }

    private:
    size_t L;
    bool spd;
    utils::COO<int> M = utils::COO<int>(1);
    bool setted = false;
};

}

#endif
