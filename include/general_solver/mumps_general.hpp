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

#ifndef MUMPS_GENERAL_HPP
#define MUMPS_GENERAL_HPP

#include <vector>
#include <dmumps_c.h>
#include "general_global_matrix/mumps_global_sparse.hpp"

// Recommended by the documentation
#define ICNTL( i ) icntl[ (i) - 1 ]
#define INFO( i ) info[ (i) - 1 ]

namespace general_solver{

class MUMPSGeneral{
    public:
    MUMPSGeneral() = default;
    ~MUMPSGeneral();

    void initialize_matrix(bool spd, size_t L);

    void compute();

    void solve(std::vector<double>& x);

    inline void add_element(const std::vector<double>& matrix, const std::vector<long>& pos){
        this->M.add_element(matrix, pos);
    }
    inline void add_element(const std::vector<double>& matrix, const std::vector<long>& pos_i, const std::vector<long>& pos_j){
        this->M.add_element(matrix, pos_i, pos_j);
    }
    inline void add_value(size_t i, size_t j, double val){
        this->M.add_value(i, j, val);
    }

    inline void make_zero(){
        this->M.make_zero();
    }

    inline void clear_matrix(){
        this->M.clear_matrix();
        this->setted = false;
    }

    private:
    general_global_matrix::MUMPSGlobalSparse M;
    DMUMPS_STRUC_C config;
    std::vector<double> buffer;
    size_t L = 0;
    bool factorized = false;
    bool setted = false;
};

}

#endif
