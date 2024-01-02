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

void EigenPCG::generate_matrix(const Meshing* const mesh, const size_t L, const std::vector<long>& node_positions, const std::vector<double>& density, double pc, double psi){
    this->gsm.generate(mesh, node_positions, L, density, pc, psi);
    auto& K = this->gsm.get_K();
    this->cg.compute(K);
}

void EigenPCG::calculate_displacements(std::vector<double>& load){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return;
    }

    Eigen::VectorXd f = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(load.data(), load.size());
    Eigen::VectorXd u = cg.solve(f);

    std::copy(u.cbegin(), u.cend(), load.begin());
}

}
