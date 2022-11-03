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

// TODO: add support for OpenMPI
#include <mpich-x86_64/mpi.h>
#include "finite_element/mumps_solver.hpp"
#include "logger.hpp"

namespace finite_element{

MUMPSSolver::MUMPSSolver(){
    this->config.sym = 1; // Hermitian matrix
    this->config.job = -1; // Configuration initialization
    this->config.par = 1; // Host process also does computations
    this->config.comm_fortran = -987654; // Default communicator
    // Matrix assembled in host
    this->config.ICNTL(5) = 0;
    this->config.ICNTL(18) = 0;
    // Diagonal only has non-zero elements
    this->config.ICNTL(6) = 0;
    this->config.ICNTL(8) = 0;
    // Parallel preprocessing
    this->config.ICNTL(28) = 2;
    this->config.ICNTL(7) = 7; // (if sequential)
    this->config.ICNTL(29) = 0;
    // No block format
    this->config.ICNTL(15) = 0;
    // No iterative refinement
    this->config.ICNTL(10) = 0;
    // No error analysis
    this->config.ICNTL(11) = 0;
    // In-core factorization and solution
    this->config.ICNTL(22) = 0;
    // No forward elimination
    this->config.ICNTL(32) = 0;
    // Right-hand side
    this->config.ICNTL(20) = 0; // dense
    this->config.ICNTL(21) = 0; // assembled

    dmumps_c(&this->config);
}

std::vector<double> MUMPSSolver::calculate_displacements(const Meshing* const mesh, std::vector<double> load, const std::vector<double>& density, double pc){
    if(this->current_step == 0){

        this->generate_K(mesh, density, pc);
        this->sK.to_mumps_format(this->rows, this->cols, this->vals);
        this->sK.clear(); // Spare some RAM
        this->config.job = 4; // Prepare matrix and decompose
        // Insert matrix data
        // Do this after every regeneration as the vectors may expand, which
        // will change their address.
        this->config.n = load.size();
        this->config.nz = this->vals.size();
        this->config.a = this->vals.data();
        this->config.irn = this->rows.data();
        this->config.jcn = this->cols.data();

        logger::quick_log("Decomposing...");
        dmumps_c(&this->config);

        logger::quick_log("Done.");
    }
    this->config.job = 3; // Solve using decomposed matrix
    logger::quick_log("Calculating displacements...");

    dmumps_c(&this->config);

    logger::quick_log("Done.");

    this->current_step = (this->current_step + 1) % this->steps;
   
    return load; 
}

void MUMPSSolver::insert_element_matrix(const std::vector<double>& k, const std::vector<long>& pos, const size_t n){
    (void) n;
    this->sK.insert_matrix_symmetric_mumps(k, pos);
}

}
