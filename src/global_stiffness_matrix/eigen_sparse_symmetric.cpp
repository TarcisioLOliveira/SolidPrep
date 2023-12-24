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
#include "global_stiffness_matrix/eigen_sparse_symmetric.hpp"
#include "utils/sparse_matrix.hpp"
#include <Eigen/src/SparseCore/SparseMatrix.h>

namespace global_stiffness_matrix{

namespace internal{

class EigenSparseSymmetricTriplets : public GlobalStiffnessMatrix{
    public:
    typedef Eigen::Triplet<double, std::ptrdiff_t> T;

    virtual ~EigenSparseSymmetricTriplets() = default;

    virtual void generate(const Meshing * const mesh, const std::vector<long>& node_positions, const size_t matrix_width, const std::vector<double>& density, const double pc, const double psi) override;

    std::vector<T> triplets;

    protected:
    bool first_time = true;
    utils::SparseMatrix K;

    virtual inline void insert_element_matrix(const std::vector<double>& k, const std::vector<long>& pos) override{
        this->K.insert_matrix_symmetric_mumps(k, pos);
    }
    inline virtual void add_to_matrix(size_t i, size_t j, double val) override{
        this->K.add(i, j, val);
    }
};

void EigenSparseSymmetricTriplets::generate(const Meshing * const mesh, const std::vector<long>& node_positions, const size_t matrix_width, const std::vector<double>& density, const double pc, const double psi){
    (void) matrix_width;
    this->generate_base(mesh, node_positions, density, pc, psi);
    this->triplets = K.get_eigen_triplets();
    this->K.clear();
}

}

void EigenSparseSymmetric::generate(const Meshing* const mesh, const std::vector<long>& node_positions, const size_t matrix_width, const std::vector<double>& density, const double pc, const double psi){
    logger::quick_log("Generating stiffness matrix...");
    if(this->first_time){
        this->K = Mat(matrix_width, matrix_width);
        internal::EigenSparseSymmetricTriplets trigen;
        trigen.generate(mesh, node_positions, matrix_width, density, pc, psi);
        this->K.setFromTriplets(trigen.triplets.begin(), trigen.triplets.end());
        this->first_time = false;
    } else {
        std::fill(this->K.valuePtr(), this->K.valuePtr() + this->K.nonZeros(), 0);
        this->generate_base(mesh, node_positions, density, pc, psi);
    }
    logger::quick_log("Done.");
}

}
