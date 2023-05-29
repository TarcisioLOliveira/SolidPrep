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

void MUMPSSparseSymmetric::generate(const Meshing* const mesh, const std::vector<double>& density, const double pc, const double psi){
    logger::quick_log("Generating stiffness matrix...");
    if(density.size() == 0){
        for(auto& g : mesh->geometries){
            this->add_geometry(mesh, g);
        }
    } else {
        auto rho = density.begin();
        for(auto& g : mesh->geometries){
            if(g->do_topopt){
                this->add_geometry(mesh, g, rho, pc, psi);
            } else {
                this->add_geometry(mesh, g);
            }
        }
    }
    this->sK.to_mumps_format(this->rows, this->cols, this->vals);
    // this->sK.clear(); // Spare some RAM
    this->sK.zero();
    logger::quick_log("Done.");
}

void MUMPSSparseSymmetric::add_geometry(const Meshing* const mesh, const Geometry* const g){
    const auto D = g->materials.get_D();
    const double t = mesh->thickness;
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();

    /*
     * utils::SparseMatrix tmp;

     * #pragma omp declare reduction(merge : utils::SparseMatrix :\
     *    omp_out.merge(omp_in))\
     *    initializer (omp_priv=(omp_orig))

     * #pragma omp parallel for reduction(merge:tmp)
     */
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
        this->sK.insert_matrix_symmetric_mumps(k, u_pos);
    }
}

void MUMPSSparseSymmetric::add_geometry(const Meshing* const mesh, const Geometry* const g, std::vector<double>::const_iterator& rho, const double pc, const double psi){
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
            const double rhop = K_MIN + (1-K_MIN)*std::pow(*rho, pc);
            for(size_t i = 0; i < D.size(); ++i){
                rhoD[i] = rhop*D[i];
            }
            const std::vector<double> k = e->get_k(rhoD, t);
            this->sK.insert_matrix_symmetric_mumps(k, u_pos);
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
            this->sK.insert_matrix_symmetric_mumps(k, u_pos);
        }
    }
}

}
