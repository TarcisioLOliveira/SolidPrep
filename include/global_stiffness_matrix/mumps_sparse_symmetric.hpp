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

#ifndef MUMPS_SPARSE_SYMMETRIC_HPP
#define MUMPS_SPARSE_SYMMETRIC_HPP

#include "utils/coo.hpp"
#include "meshing.hpp"
#include "global_stiffness_matrix.hpp"

namespace global_stiffness_matrix{

class MUMPSSparseSymmetric : public GlobalStiffnessMatrix{
    public:
    virtual ~MUMPSSparseSymmetric() = default;

    virtual void generate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache, const FiniteElement::MatrixType type) override;

    inline virtual void dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const override{
        this->sK.dot_vector(v, v_out, false);
    }

    inline std::vector<int>& get_rows(){
        return this->sK.rows;
    }
    inline std::vector<int>& get_cols(){
        return this->sK.cols;
    }
    inline std::vector<double>& get_vals(){
        return this->sK.vals;
    }

    protected:
    bool first_time = true;
    utils::COO<int> sK = utils::COO<int>(1);

    inline virtual void insert_block_symmetric(const std::vector<double>& k, const std::vector<long>& posi, const std::vector<long>& posj) override{
        if(posi[0] > posj[0]){
            this->sK.insert_block(k, posi, posj, false);
        } else {
            this->sK.insert_block(k, posi, posj, true);
        }
    }
    inline virtual void insert_element_matrix(const std::vector<double>& k, const std::vector<long>& pos) override{
        this->sK.insert_matrix_symmetric(k, pos);
    }
    inline virtual void add_to_matrix(size_t i, size_t j, double val) override{
        this->sK.add(i, j, val);
    }
};

}

#endif
