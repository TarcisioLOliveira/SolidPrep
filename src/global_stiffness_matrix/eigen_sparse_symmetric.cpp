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
#include "math/matrix.hpp"
#include "utils/sparse_matrix.hpp"
#include <Eigen/src/SparseCore/SparseMatrix.h>

namespace global_stiffness_matrix{

namespace internal{

class EigenSparseSymmetricTriplets : public GlobalStiffnessMatrix{
    public:
    typedef Eigen::Triplet<double, std::ptrdiff_t> T;

    EigenSparseSymmetricTriplets():GlobalStiffnessMatrix(1){}

    virtual ~EigenSparseSymmetricTriplets() = default;

    virtual void generate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type) override;

    inline virtual void dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const override{
        (void)v;
        (void)v_out;
        // Unused
    }
    inline virtual void reset_hessian() override{
        // Unused 
    };

    std::vector<T> triplets;

    protected:
    bool first_time = true;
    utils::SparseMatrix K;

    inline virtual void insert_block_symmetric(const math::Matrix& k, const std::vector<long>& posi, const std::vector<long>& posj) override{
        // UNTESTED HEURISTIC!!!!!!!
        if(posi[0] > posj[0]){
            this->K.insert_block(k, posi, posj, false);
        } else {
            this->K.insert_block(k, posi, posj, true);
        }
    }
    virtual inline void insert_element_matrix(const math::Matrix& k, const std::vector<long>& pos) override{
        this->K.insert_matrix_symmetric_mumps(k, pos);
    }
    inline virtual void reserve_block_symmetric(const math::Matrix& k, const std::vector<long>& posi, const std::vector<long>& posj) override{
        // UNTESTED HEURISTIC!!!!!!!
        if(posi[0] > posj[0]){
            this->K.reserve_block(k, posi, posj, false);
        } else {
            this->K.reserve_block(k, posi, posj, true);
        }
    }
    virtual inline void reserve_element_matrix(const math::Matrix& k, const std::vector<long>& pos) override{
        this->K.reserve_matrix_symmetric_mumps(k, pos);
    }
    inline virtual void add_to_matrix(size_t i, size_t j, double val) override{
        this->K.add(i, j, val);
    }
};

void EigenSparseSymmetricTriplets::generate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type){
    (void) u_size;
    (void) l_num;
    this->generate_base(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, lambda, type);
    this->triplets = K.get_eigen_triplets();
    this->K.clear();
}

}

void EigenSparseSymmetric::generate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type){
    logger::quick_log("Generating stiffness matrix...");
    this->u_size = u_size;
    this->l_num = l_num;
    if(this->first_time){
        size_t matrix_width = u_size;
        if(type != FiniteElement::ContactType::RIGID){
            matrix_width += l_num;
        }
        this->K = Mat(matrix_width, matrix_width);
        internal::EigenSparseSymmetricTriplets trigen;
        trigen.generate(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, lambda, type);
        this->K.setFromTriplets(trigen.triplets.begin(), trigen.triplets.end());
        this->first_time = false;
    } else {
        std::fill(this->K.valuePtr(), this->K.valuePtr() + this->K.nonZeros(), 0);
        this->generate_base(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, lambda, type);
    }
    if(type > FiniteElement::ContactType::RIGID){
        this->K_bkp = K;
    }
    logger::quick_log("Done.");
}

}
