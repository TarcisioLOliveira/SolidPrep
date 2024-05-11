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

    virtual void generate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache, const FiniteElement::MatrixType type) override;

    std::vector<T> triplets;

    protected:
    bool first_time = true;
    utils::SparseMatrix K;

    inline virtual void insert_block_symmetric(const std::vector<double>& k, const std::vector<long>& posi, const std::vector<long>& posj) override{
        this->K.insert_block(k, posi, posj);
        this->K.insert_block(k, posj, posi);
    }
    virtual inline void insert_element_matrix(const std::vector<double>& k, const std::vector<long>& pos) override{
        this->K.insert_matrix_general_mumps(k, pos);
    }
    inline virtual void add_to_matrix(size_t i, size_t j, double val) override{
        this->K.add(i, j, val);
    }
};

void EigenSparseAsymmetricTriplets::generate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache, const FiniteElement::MatrixType type){
    (void) u_size;
    (void) l_num;
    this->generate_base(mesh, u_size, l_num, node_positions, topopt, D_cache, type);
    this->triplets = K.get_eigen_triplets();
    this->K.clear();
}

}

void EigenSparseAsymmetric::generate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache, const FiniteElement::MatrixType type){
    logger::quick_log("Generating stiffness matrix...");
    if(this->first_time){
        const size_t matrix_width = u_size + 2*l_num;
        this->K = Mat(matrix_width, matrix_width);
        internal::EigenSparseAsymmetricTriplets trigen;
        trigen.generate(mesh, u_size, l_num, node_positions, topopt, D_cache, type);
        this->K.setFromTriplets(trigen.triplets.begin(), trigen.triplets.end());
        this->first_time = false;
    } else {
        std::fill(this->K.valuePtr(), this->K.valuePtr() + this->K.nonZeros(), 0);
        this->generate_base(mesh, u_size, l_num, node_positions, topopt, D_cache, type);
    }
    logger::quick_log("Done.");
}

}
