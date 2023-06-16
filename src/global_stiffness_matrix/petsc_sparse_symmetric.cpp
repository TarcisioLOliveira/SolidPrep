/*
 *   Copyright (C) 2023 Tarcísio Ladeia de Oliveira.
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

#include <mpich-x86_64/mpi.h>
#include "global_stiffness_matrix/petsc_sparse_symmetric.hpp"
#include "logger.hpp"


namespace global_stiffness_matrix{

PETScSparseSymmetric::PETScSparseSymmetric(Backend backend){
    switch(backend){
        case Backend::CPU:
            this->mat_type = MATMPIAIJ;
            break;
        case Backend::CUDA:
            this->mat_type = MATAIJCUSPARSE;
            break;
    }
}


PETScSparseSymmetric::~PETScSparseSymmetric(){
    MatDestroy(&this->K);
}

void PETScSparseSymmetric::generate(const Meshing* const mesh, const std::vector<double>& density, const double pc, const double psi){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id == 0){
        logger::quick_log("Generating stiffness matrix...");
    }

    if(first_time){
        MatCreate(PETSC_COMM_WORLD, &this->K);
        MatSetType(this->K, this->mat_type.c_str());

        long M = mesh->load_vector.size();
        MPI_Bcast(&M, 1, MPI_LONG, 0, MPI_COMM_WORLD);

        MatSetSizes(this->K, PETSC_DECIDE, PETSC_DECIDE, M, M);
        MatSetOption(this->K, MAT_SPD, PETSC_TRUE);

        Mat tmp;
        MatCreate(PETSC_COMM_WORLD, &tmp);
        MatSetType(tmp, MATPREALLOCATOR);
        MatSetSizes(tmp, PETSC_DECIDE, PETSC_DECIDE, M, M);
        MatSetOption(tmp, MAT_SPD, PETSC_TRUE);
        MatSetUp(tmp);

        if(mpi_id == 0){
            std::swap(tmp, K);
            this->fill_matrix(mesh, density, pc, psi);
            std::swap(tmp, K);
        }

        MatAssemblyBegin(tmp, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(tmp, MAT_FINAL_ASSEMBLY);

        MatPreallocatorPreallocate(tmp, PETSC_TRUE, this->K);

        MatSetUp(this->K);

        MatDestroy(&tmp);
        //MatSetSizes(this->K, mesh->load_vector.size(), mesh->load_vector.size(), mesh->load_vector.size(), mesh->load_vector.size());
        //MatSetUp(this->K);
        this->first_time = false;
    } else {
        MatZeroEntries(this->K);
    }

    if(mpi_id == 0){
        this->fill_matrix(mesh, density, pc, psi);
    }

    MatAssemblyBegin(this->K, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(this->K, MAT_FINAL_ASSEMBLY);
    if(mpi_id == 0){
        logger::quick_log("Done.");
    }
}

void PETScSparseSymmetric::add_geometry(const Meshing* const mesh, const Geometry* const g){
    const auto D = g->materials.get_D();
    const double t = mesh->thickness;
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();

    for(auto& e : g->mesh){
        std::vector<long> u_pos;
        u_pos.reserve(dof*node_num);
        for(size_t i = 0; i < node_num; ++i){
            const auto& n = e->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                u_pos.push_back(n->u_pos[j]);
            }
        }
        std::vector<double> k = e->get_k(D, t);
        this->insert_element_matrix(k, u_pos);
    }
}

void PETScSparseSymmetric::add_geometry(const Meshing* const mesh, const Geometry* const g, std::vector<double>::const_iterator& rho, const double pc, const double psi){
    const double t = mesh->thickness;
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();

    //const size_t num_den = g->number_of_densities_needed();
    const size_t num_mat = g->number_of_materials();
    if(num_mat == 1){
        const auto D = g->materials.get_D();
        auto rhoD = D;
        for(const auto& e:g->mesh){
            std::vector<long> u_pos;
            u_pos.reserve(dof*node_num);
            for(size_t i = 0; i < node_num; ++i){
                const auto& n = e->nodes[i];
                for(size_t j = 0; j < dof; ++j){
                    u_pos.push_back(n->u_pos[j]);
                }
            }
            const double rhop = std::pow(*rho, pc);
            for(size_t i = 0; i < D.size(); ++i){
                rhoD[i] = rhop*D[i];
            }
            const std::vector<double> k = e->get_k(rhoD, t);
            this->insert_element_matrix(k, u_pos);
            ++rho;
        }
    } else {
        auto D = g->materials.get_D();
        for(const auto& e:g->mesh){
            std::vector<long> u_pos;
            u_pos.reserve(dof*node_num);
            for(size_t i = 0; i < node_num; ++i){
                const auto& n = e->nodes[i];
                for(size_t j = 0; j < dof; ++j){
                    u_pos.push_back(n->u_pos[j]);
                }
            }
            g->materials.get_D(rho, g->with_void, pc, this->K_MIN, psi, D);
            const std::vector<double> k = e->get_k(D, t);
            this->insert_element_matrix(k, u_pos);
        }
    }
}

}
