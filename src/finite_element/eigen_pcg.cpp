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

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/src/Core/util/Constants.h>
#include <mpich-x86_64/mpi.h>
#include "finite_element/eigen_pcg.hpp"
#include "logger.hpp"

namespace finite_element{

EigenPCG::EigenPCG():gsm(){}

void EigenPCG::generate_matrix_base(const Meshing* const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache, const MatrixType type){

    this->l_num = l_num;
    this->gsm.generate(mesh, u_size, l_num, node_positions, topopt, D_cache, type);
    auto& K = this->gsm.get_K();
    this->cg.compute(K);
}

void EigenPCG::solve(std::vector<double>& load, std::vector<double>& lambda){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
    // TODO: lambda

    if(mpi_id != 0){
        return;
    }

    Eigen::VectorXd f(load.size() + 2*l_num);
    std::copy(load.begin(), load.end(), f.begin());
    Eigen::VectorXd u = cg.solve(f);

    std::copy(u.cbegin(), u.cbegin() + load.size(), load.begin());
    if(l_num > 0){
        std::copy(u.cbegin() + load.size(), u.cend(), lambda.begin());
    }
}

}
