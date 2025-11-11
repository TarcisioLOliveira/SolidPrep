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

#ifndef PETSC_GLOBAL_SPARSE_HPP
#define PETSC_GLOBAL_SPARSE_HPP

#include "math/matrix.hpp"
#include "utils/coo.hpp"
#include <vector>
#include <petsc.h>

namespace general_global_matrix{

class PETScGlobalSparse{
    public:
    ~PETScGlobalSparse();

    void initialize(size_t L);

    inline void add_element(const math::Matrix& matrix, const std::vector<long>& pos){
        this->M.insert_matrix_general(matrix, pos);
    }
    inline void add_element(const math::Matrix& matrix, const std::vector<long>& pos_i, const std::vector<long>& pos_j){
        this->M.insert_block(matrix, pos_i, pos_j, false);
    }

    inline void print_matrix() const{
        for(size_t i = 0; i < L; ++i){
            for(size_t j = 0; j < L; ++j){
                std::cout << this->M.get(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }
    inline void add_value(size_t i, size_t j, double val){
        this->M.add(i, j, val);
    }

    inline void set_up(){
        if(!this->setted){
            if(M.vals.size() == 0){
                this->M.generate_coo(this->L);
                MatSetPreallocationCOO(this->PM, this->M.nvals, this->M.rows.data(), this->M.cols.data());
            }
            MatSetValuesCOO(this->PM, this->M.vals.data(), INSERT_VALUES);
            this->setted = true;
        }
    }

    inline void make_zero(){
        this->M.zero();
    }

    inline void clear_matrix(){
        this->M.clear();
        setted = false;
    }

    Mat get_matrix() const{
        return PM;
    }

    private:
    Mat PM;
    Mat tmp;
    utils::COO<PetscInt> M = utils::COO<PetscInt>(0);
    size_t L;
    bool setted = false;
};

}

#endif
