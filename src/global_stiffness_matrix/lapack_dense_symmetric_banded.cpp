/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
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

#include "global_stiffness_matrix/lapack_dense_symmetric_banded.hpp"
#include "logger.hpp"

namespace global_stiffness_matrix{

void LAPACKDenseSymmetricBanded::generate(const Meshing* const mesh, const std::vector<long>& node_positions, const size_t matrix_width, const std::vector<double>& density, const double pc, const double psi){
    if(this->first_time){
        this->K.clear();
        this->calculate_dimensions(mesh, node_positions, matrix_width);
        this->K.resize(W*N,0);
        this->first_time = false;
    } else {
        std::fill(this->K.begin(), this->K.end(), 0);
    }
    logger::quick_log("Generating stiffness matrix...");
    this->generate_base(mesh, node_positions, density, pc, psi);
    logger::quick_log("Done.");
}


}
