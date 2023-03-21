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

#ifndef EIGEN_SPARSE_SYMMETRIC_HPP
#define EIGEN_SPARSE_SYMMETRIC_HPP

#include "meshing.hpp"
#include <Eigen/SparseCore>
#include <Eigen/src/SparseCore/SparseMatrix.h>

namespace global_stiffness_matrix{

class EigenSparseSymmetric{
    public:
    const double K_MIN = 1e-6;
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor, std::ptrdiff_t> Mat;

    virtual ~EigenSparseSymmetric() = default;

    virtual void generate(const Meshing * const mesh, const std::vector<double>& density, const double pc, const double psi);

    Mat& get_K() {
        return K;
    }

    protected:
    bool first_time = true;
    Mat K;
    size_t W, N;


    virtual void calculate_dimensions(const Meshing * const mesh, const std::vector<double>& load);

    virtual void add_geometry(const Meshing * const mesh, const Geometry * const g);

    virtual void add_geometry(const Meshing * const mesh, const Geometry * const g, std::vector<double>::const_iterator& rho, const double pc, const double psi);

    inline void insert_element_matrix(const std::vector<double>& k, const std::vector<long>& pos){
        const size_t w = pos.size();
        for(size_t i = 0; i < w; ++i){
            for(size_t j = i; j < w; ++j){
                if(pos[i] > -1 && pos[j] > -1){
                    if(pos[j] > pos[i]){
                        K.coeffRef(pos[i], pos[j]) += k[w*i + j];
                    } else {
                        K.coeffRef(pos[j], pos[i]) += k[w*i + j];
                    }
                }
            }
        }
    }
};

}

#endif
