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
#include "utils/sparse_matrix.hpp"

namespace finite_element{

MUMPSSolver::MUMPSSolver(){
    this->config.sym = 1; // Hermitian matrix
    this->config.job = -1; // Configuration initialization
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
        this->config.nnz = this->vals.size();
        this->config.a = this->vals.data();
        this->config.irn = this->rows.data();
        this->config.jcn = this->cols.data();

        logger::quick_log("Decomposing...");
        dmumps_c(&this->config);

        logger::quick_log("Done.");
    }

    this->config.job = 3; // Solve using decomposed matrix
    this->config.rhs = load.data(); // Set right-hand side
    logger::quick_log("Calculating displacements...");

    dmumps_c(&this->config);

    logger::quick_log("Done.");

    this->current_step = (this->current_step + 1) % this->steps;
   
    return load; 
}

void MUMPSSolver::_add_geometry_to_K(const Meshing* const mesh, const Geometry* const g){
    const auto D = g->get_D(0);
    const double t = mesh->thickness;
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();

    utils::SparseMatrix tmp;

    #pragma omp declare reduction(merge : utils::SparseMatrix :\
       omp_out.merge(omp_in))\
       initializer (omp_priv=(omp_orig))

    #pragma omp parallel for reduction(merge:tmp)
    for(size_t id = 0; id < g->mesh.size(); ++id){
        const auto& e = g->mesh[id];
        std::vector<long> u_pos;
        u_pos.reserve(dof*node_num);
        for(size_t i = 0; i < node_num; ++i){
            const auto& n = e->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                u_pos.push_back(n->u_pos[j]);
            }
        }
        std::vector<double> k = e->get_k(D, t);
        tmp.insert_matrix_symmetric_mumps(k, u_pos);
    }
    this->sK = std::move(tmp);
}

void MUMPSSolver::_add_geometry_to_K(const Meshing* const mesh, const Geometry* const g, std::vector<double>::const_iterator& rho, const double pc){
    const double t = mesh->thickness;
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();

    const size_t num_den = g->number_of_densities_needed();
    // const size_t num_mat = g->number_of_materials();
    if(num_den == 1){
        for(const auto& e:g->mesh){
            std::vector<long> u_pos;
            u_pos.reserve(dof*node_num);
            for(size_t i = 0; i < node_num; ++i){
                const auto& n = e->nodes[i];
                for(size_t j = 0; j < dof; ++j){
                    u_pos.push_back(n->u_pos[j]);
                }
            }
            const auto D = g->get_D_topopt(*rho, pc, this->K_MIN);
            const std::vector<double> k = e->get_k(D, t);
            this->sK.insert_matrix_symmetric_mumps(k, u_pos);
            ++rho;
        }
    } else {
        logger::log_assert(num_den == 1, logger::ERROR, "FEA problems that require more than 1 design variables are currently not supported.");
    }
}

void MUMPSSolver::insert_element_matrix(const std::vector<double>& k, const std::vector<long>& pos, const size_t n){
    (void) n;
    this->sK.insert_matrix_symmetric_mumps(k, pos);
}

}
