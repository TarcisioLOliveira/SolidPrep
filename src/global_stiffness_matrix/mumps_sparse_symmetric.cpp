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

void MUMPSSparseSymmetric::generate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache, const std::vector<double>& u_ext, const FiniteElement::ContactType type){
    logger::quick_log("Generating stiffness matrix...");
    this->u_size = u_size;
    this->l_num = l_num;
    this->sK.zero();
    this->generate_base(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, type);
    size_t M = u_size;
    if(type != FiniteElement::ContactType::RIGID){
        M += l_num;
    }
    if(this->first_time){
        this->sK.generate_coo(M);
        this->first_time = false;
    }
    if(type == FiniteElement::ContactType::FRICTIONLESS_DISPL){
        this->sK.backup_matrix();
    }
    logger::quick_log("Done.");
}

void MUMPSSparseSymmetric::reset_hessian(){
    this->sK.restore_matrix();
}

}
