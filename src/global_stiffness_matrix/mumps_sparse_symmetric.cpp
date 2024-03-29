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

#include "logger.hpp"
#include "global_stiffness_matrix/mumps_sparse_symmetric.hpp"

namespace global_stiffness_matrix{

void MUMPSSparseSymmetric::generate(const Meshing* const mesh, const std::vector<long>& node_positions, const size_t matrix_width, bool topopt, const std::vector<std::vector<double>>& D_cache) {
    logger::quick_log("Generating stiffness matrix...");
    this->sK.zero();
    this->generate_base(mesh, node_positions, topopt, D_cache);
    if(this->first_time){
        this->sK.generate_coo(matrix_width);
        this->first_time = false;
    }
    logger::quick_log("Done.");
}

}
