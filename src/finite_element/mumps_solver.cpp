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
#include "math/matrix.hpp"
#include "project_data.hpp"

namespace finite_element{

MUMPSSolver::MUMPSSolver(const projspec::DataMap& data):
    FiniteElement(data.proj->contact_data),
    gsm(data.proj->contact_data.EPS_DISPL_SIMPLE)
{
    this->set_global_matrix(&this->gsm);

    this->config.job = -1; // Configuration initialization
    this->config.sym = 1; // Hermitian matrix
    this->config.par = 1; // Host process also does computations
    this->config.comm_fortran = -987654; // Default communicator
    // No text output
    this->config.ICNTL(1) = 0;
    this->config.ICNTL(2) = 0;
    this->config.ICNTL(3) = 0;
    this->config.ICNTL(4) = 0;
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
    // Complete factorization
    this->config.ICNTL(19) = 0;
    dmumps_c(&this->config);
    // No text output, but for real now
    this->config.ICNTL(1) = 0;
    this->config.ICNTL(2) = 0;
    this->config.ICNTL(3) = 0;
    this->config.ICNTL(4) = 0;
}

MUMPSSolver::~MUMPSSolver(){
    // Clean up. Especially necessary if using out-of-core memory
    this->config.job = -2;

    dmumps_c(&this->config);
}

void MUMPSSolver::generate_matrix_base(const Meshing* const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const ContactType type){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    long M = u_size + l_num;
    if(mpi_id == 0){
        this->gsm.generate(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, lambda, type);

        std::vector<int>& rows = this->gsm.get_rows();
        std::vector<int>& cols = this->gsm.get_cols();
        std::vector<double>& vals = this->gsm.get_vals();
        // Insert matrix data
        // Do this after every regeneration as the vectors may expand, which
        // will change their address.
        this->config.n   = M;
        this->config.nnz = vals.size();
        this->config.a   = vals.data();
        this->config.irn = rows.data();
        this->config.jcn = cols.data();
        this->factorized = false;
    }
    if(this->first_time){
        const double buffer_mult = 1.35;
        this->config.job = 1; // Perform analysis
        dmumps_c(&this->config);
        int size = this->config.INFO(8);
        if(size > 0){
            this->buffer.resize(size*buffer_mult);
            if(this->buffer.size() < 1e6){
                this->config.lwk_user = buffer.size();
            } else {
                this->config.lwk_user = -double(buffer.size())/1e6;
            }
        } else {
            this->buffer.resize(-double(size)*1e6*buffer_mult);
            this->config.lwk_user = -double(this->buffer.size())/1e6;
        }
        this->config.ICNTL(14) = 35;
        this->config.ICNTL(23) = this->config.INFO(15)*1.2 - double(this->buffer.size())/1e6;;
        this->config.wk_user = this->buffer.data();
        this->first_time = false;
    }
}

void MUMPSSolver::solve(std::vector<double>& load){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(!this->factorized){
        if(mpi_id == 0){
            logger::quick_log("Factorizing...");
        }

        this->config.job = 2; // Factorize
        this->config.ICNTL(33) = 1;
        dmumps_c(&this->config);
        logger::quick_log("DET", this->config.rinfog[12-1], this->config.rinfog[34-1]);
        this->factorized = true;

        if(mpi_id == 0){
            logger::quick_log("Done.");
        }
    }

    this->config.job = 3; // Solve using decomposed matrix
    if(mpi_id == 0){
        logger::log_assert(load.size() == (size_t)this->config.n, logger::ERROR, "invalid vector input vector dimension, is {}, should be {}", load.size(), this->config.n);
        this->config.rhs = load.data(); // Set right-hand side
        logger::quick_log("Calculating displacements...");
    }

    dmumps_c(&this->config);

    if(mpi_id == 0){
        logger::quick_log("Done.");
    }
}

void MUMPSSolver::reset_hessian(){
    this->gsm.reset_hessian();
    this->factorized = false;
}

using namespace projspec;
const bool MUMPSSolver::reg = Factory<FiniteElement>::add(
    [](const DataMap& data){
        return std::make_unique<MUMPSSolver>(data);
    },
    ObjectRequirements{
        "mumps",
        {}
    }
);

}
