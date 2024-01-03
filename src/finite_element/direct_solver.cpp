/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
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


#include "finite_element/direct_solver.hpp"
#include "lapacke.h"
#include "logger.hpp"
#include "utils.hpp"
#include <limits>
#include <cblas.h>
#include <mpich-x86_64/mpi.h>
#include "project_data.hpp"

namespace finite_element{

void DirectSolver::generate_matrix(const Meshing* const mesh, const size_t L, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache){
    this->gsm.generate(mesh, node_positions, L, topopt, D_cache);
    this->factorized = false;
}

void DirectSolver::calculate_displacements(std::vector<double>& load){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return;
    }

    const size_t& W = this->gsm.get_W();
    const size_t& N = this->gsm.get_N();
    std::vector<double>& K = this->gsm.get_K();

    if(!this->factorized){

        logger::quick_log("Decomposing...");
        int info = LAPACKE_dpbtrf_work(LAPACK_COL_MAJOR, 'L', W, N-1, K.data(), N);
        logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while factoring stiffness matrix.", info);
        logger::quick_log("Done.");
        factorized = true;
    }
    logger::quick_log("Calculating displacements...");
    logger::quick_log("W: ",W," N: ", N);

    int info = LAPACKE_dpbtrs_work(LAPACK_COL_MAJOR, 'L', W, N-1, 1, K.data(), N, load.data(), W);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating displacements.", info);
    // int info = LAPACKE_dpbsv_work(LAPACK_COL_MAJOR, 'L', W, N-1, 1, K.data(), N, U.data(), W);
    // logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating displacements.", info);

    logger::quick_log("Done.");
}

}
