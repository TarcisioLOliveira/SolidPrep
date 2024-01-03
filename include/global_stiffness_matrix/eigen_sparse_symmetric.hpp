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

#include <Eigen/SparseCore>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include "meshing.hpp"
#include "global_stiffness_matrix.hpp"

namespace global_stiffness_matrix{

class EigenSparseSymmetric : public GlobalStiffnessMatrix{
    public:
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor, std::ptrdiff_t> Mat;

    virtual ~EigenSparseSymmetric() = default;

    virtual void generate(const Meshing * const mesh, const std::vector<long>& node_positions, const size_t matrix_width, bool topopt, const std::vector<std::vector<double>>& D_cache) override;

    Mat& get_K() {
        return K;
    }

    protected:
    bool first_time = true;
    Mat K;

    inline void insert_element_matrix(const std::vector<double>& k, const std::vector<long>& pos) override{
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
    inline virtual void add_to_matrix(size_t i, size_t j, double val) override{
        K.coeffRef(i, j) += val;
    }
};

}

#endif
