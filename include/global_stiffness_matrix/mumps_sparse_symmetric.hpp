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

    virtual void generate(const Meshing * const mesh, const std::vector<long>& node_positions, const size_t matrix_width, const std::vector<double>& density, const double pc, const double psi) override;

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

    virtual inline void insert_element_matrix(const std::vector<double>& k, const std::vector<long>& pos) override{
        this->sK.insert_matrix_symmetric(k, pos);
    }
    inline virtual void add_to_matrix(size_t i, size_t j, double val) override{
        this->sK.add(i, j, val);
    }
};

}

#endif
