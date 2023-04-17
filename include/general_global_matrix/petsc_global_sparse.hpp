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

#include <vector>
#include <petsc.h>

namespace general_global_matrix{

class PETScGlobalSparse{
    public:
    ~PETScGlobalSparse();

    void initialize(size_t L);

    void begin_preallocation();

    void end_preallocation();

    inline void add_element(const std::vector<double> matrix, const std::vector<long> pos){
        MatSetValues(this->M, pos.size(), pos.data(), pos.size(), pos.data(), matrix.data(), ADD_VALUES);
    }

    inline void clear_matrix(){
        MatZeroEntries(this->M);
    }

    inline void set_up(){
        MatAssemblyBegin(this->M, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(this->M, MAT_FINAL_ASSEMBLY);
    }

    Mat get_matrix() const{
        return M;
    }

    private:
    Mat M;
    Mat tmp;
    size_t L;

};

}

#endif
