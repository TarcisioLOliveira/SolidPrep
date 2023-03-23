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

#ifndef PETSC_SPARSE_SYMMETRIC_HPP
#define PETSC_SPARSE_SYMMETRIC_HPP

#include <vector>
#include <petsc.h>
#include "meshing.hpp"

namespace global_stiffness_matrix{

class PETScSparseSymmetric{
    public:
    const double K_MIN = 1e-6;

    PETScSparseSymmetric();

    virtual ~PETScSparseSymmetric();

    virtual void generate(const Meshing * const mesh, const std::vector<double>& density, const double pc, const double psi);

    Mat get_K() const{
        return this->K;
    }

    protected:
    Mat K = 0;
    bool first_time = true;

    virtual void add_geometry(const Meshing * const mesh, const Geometry * const g);

    virtual void add_geometry(const Meshing * const mesh, const Geometry * const g, std::vector<double>::const_iterator& rho, const double pc, const double psi);

    inline void fill_matrix(const Meshing * const mesh, const std::vector<double>& density, const double pc, const double psi){
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
    }

    inline void insert_element_matrix(const std::vector<double>& k, const std::vector<long>& pos){
        // Requires 64-bit indices
        MatSetValues(this->K, pos.size(), pos.data(), pos.size(), pos.data(), k.data(), ADD_VALUES);
        //const size_t w = pos.size();
        //for(size_t i = 0; i < w; ++i){
        //    for(size_t j = i; j < w; ++j){
        //        if(pos[i] > -1 && pos[j] > -1){
        //           // if(pos[j] > pos[i]){
        //           //     K.coeffRef(pos[i], pos[j]) += k[w*i + j];
        //           // } else {
        //           //     K.coeffRef(pos[j], pos[i]) += k[w*i + j];
        //           // }
        //        }
        //    }
        //}
    }
};

}

#endif
