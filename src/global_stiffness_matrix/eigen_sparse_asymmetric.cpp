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

#include "global_stiffness_matrix/eigen_sparse_asymmetric.hpp"
#include "logger.hpp"
#include "utils/sparse_matrix.hpp"

namespace global_stiffness_matrix{

namespace internal{

class EigenSparseAsymmetricTriplets : public GlobalStiffnessMatrix{
    public:
    typedef Eigen::Triplet<double, std::ptrdiff_t> T;

    virtual ~EigenSparseAsymmetricTriplets() = default;

    virtual void generate(const Meshing * const mesh, const std::vector<long>& node_positions, const size_t matrix_width, bool topopt, const std::vector<std::vector<double>>& D_cache) override;

    std::vector<T> triplets;

    protected:
    bool first_time = true;
    utils::SparseMatrix K;

    virtual inline void insert_element_matrix(const std::vector<double>& k, const std::vector<long>& pos) override{
        this->K.insert_matrix_general_mumps(k, pos);
    }
    inline virtual void add_to_matrix(size_t i, size_t j, double val) override{
        this->K.add(i, j, val);
    }
};

void EigenSparseAsymmetricTriplets::generate(const Meshing * const mesh, const std::vector<long>& node_positions, const size_t matrix_width, bool topopt, const std::vector<std::vector<double>>& D_cache){
    (void)matrix_width;
    this->generate_base(mesh, node_positions, topopt, D_cache);
    this->triplets = K.get_eigen_triplets();
    this->K.clear();
}

}

void EigenSparseAsymmetric::generate(const Meshing* const mesh, const std::vector<long>& node_positions, const size_t matrix_width, bool topopt, const std::vector<std::vector<double>>& D_cache){
    logger::quick_log("Generating stiffness matrix...");
    if(this->first_time){
        this->K = Mat(matrix_width, matrix_width);
        internal::EigenSparseAsymmetricTriplets trigen;
        trigen.generate(mesh, node_positions, matrix_width, topopt, D_cache);
        this->K.setFromTriplets(trigen.triplets.begin(), trigen.triplets.end());
        this->first_time = false;
    } else {
        std::fill(this->K.valuePtr(), this->K.valuePtr() + this->K.nonZeros(), 0);
        this->generate_base(mesh, node_positions, topopt, D_cache);
    }
    logger::quick_log("Done.");
}

}
