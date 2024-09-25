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

void MUMPSSparseSymmetric::generate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache, const std::vector<double>& u_ext, const FiniteElement::MatrixType type){
    logger::quick_log("Generating stiffness matrix...");
    this->u_size = u_size;
    this->l_num = l_num;
    this->sK.zero();
    this->generate_base(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, type);
    size_t M = u_size;
    if(type == FiniteElement::MatrixType::FRICTIONLESS){
        M += l_num;
    }
    if(this->first_time){
        this->sK.generate_coo(M);
        this->first_time = false;
    }
    if(type != FiniteElement::MatrixType::RIGID){
        this->sK.backup_matrix();
    }
    logger::quick_log("Done.");
}

bool MUMPSSparseSymmetric::generate_hessian(std::vector<double>& lambda, const std::vector<double>& Ku){
    return this->sK.generate_hessian(this->u_size + 2*this->l_num, this->l_num, lambda, Ku, false);
}

void MUMPSSparseSymmetric::reset_hessian(){
    this->sK.restore_matrix();
}

double MUMPSSparseSymmetric::get_newton_step(const std::vector<double>& delta, const std::vector<double>& lambda, const std::vector<double>& Ku){
    const size_t hoffset = this->u_size + 2*this->l_num;
    double M = 1.0;
    for(size_t i = 0; i < l_num; ++i){
        const size_t ui = hoffset + i;
        const size_t li = 2*l_num + i;
        if(Ku[ui] > 0 || std::abs(delta[ui]) < 1e-14){
            continue;
        }
        double M_test = (std::sqrt((this->K_MIN - 2*Ku[ui])/this->sK.get(ui,ui)) - lambda[li])/delta[ui];
        if(M_test > 0 && M_test < M){
            M = M_test;
        }
    }
    return M;
}

}
