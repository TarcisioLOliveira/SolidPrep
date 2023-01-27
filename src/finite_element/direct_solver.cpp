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

std::vector<double> DirectSolver::calculate_displacements(const Meshing* const mesh, std::vector<double> load, const std::vector<double>& density, double pc, double psi){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return std::vector<double>();
    }

    const size_t& W = this->gsm.get_W();
    const size_t& N = this->gsm.get_N();
    std::vector<double>& K = this->gsm.get_K();

    if(W == 0 || N == 0){
        this->gsm.calculate_dimensions(mesh, load);
    }

    if(this->current_step == 0){

        this->gsm.generate(mesh, density, pc, psi);

        logger::quick_log("Decomposing...");
        int info = LAPACKE_dpbtrf_work(LAPACK_COL_MAJOR, 'L', W, N-1, K.data(), N);
        logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while factoring stiffness matrix.", info);
        logger::quick_log("Done.");
    }
    logger::quick_log("Calculating displacements...");
    logger::quick_log("W: ",W," N: ", N);

    int info = LAPACKE_dpbtrs_work(LAPACK_COL_MAJOR, 'L', W, N-1, 1, K.data(), N, load.data(), W);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating displacements.", info);
    // int info = LAPACKE_dpbsv_work(LAPACK_COL_MAJOR, 'L', W, N-1, 1, K.data(), N, U.data(), W);
    // logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating displacements.", info);

    logger::quick_log("Done.");

    this->current_step = (this->current_step + 1) % this->steps;
   
    return load; 
}

}
