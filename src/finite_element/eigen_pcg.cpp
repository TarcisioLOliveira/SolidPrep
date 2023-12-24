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

EigenPCG::EigenPCG():gsm(), u(1){}

std::vector<double> EigenPCG::calculate_displacements(const Meshing* const mesh, const std::vector<long>& node_positions, std::vector<double> load, const std::vector<double>& density, double pc, double psi){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return std::vector<double>();
    }

    if(this->current_step == 0){
        this->gsm.generate(mesh, node_positions, load.size(), density, pc, psi);
        auto& K = this->gsm.get_K();
        this->cg.compute(K);
    }

    Eigen::VectorXd f = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(load.data(), load.size());
    if(u[current_step].size() == 0){
        u[current_step] = f;
    }
    u[current_step] = cg.solveWithGuess(f, u[current_step]);

    std::copy(u[current_step].cbegin(), u[current_step].cend(), load.begin());

    this->current_step = (this->current_step + 1) % this->steps;
   
    return load; 
}

}
