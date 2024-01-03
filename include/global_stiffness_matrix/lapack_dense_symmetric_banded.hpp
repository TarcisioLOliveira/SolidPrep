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

#ifndef LAPACK_DENSE_SYMMETRIC_BANDED_HPP
#define LAPACK_DENSE_SYMMETRIC_BANDED_HPP

#include "meshing.hpp"
#include "global_stiffness_matrix.hpp"

namespace global_stiffness_matrix{

class LAPACKDenseSymmetricBanded : public GlobalStiffnessMatrix{
    public:
    virtual ~LAPACKDenseSymmetricBanded() = default;

    virtual void generate(const Meshing * const mesh, const std::vector<long>& node_positions, const size_t matrix_width, bool topopt, const std::vector<std::vector<double>>& D_cache) override;

    inline std::vector<double>& get_K(){
        return K;
    }
    inline const std::vector<double>& get_K() const{
        return K;
    }

    inline const size_t& get_W() const{
        return W;
    }
    inline const size_t& get_N() const{
        return N;
    }

    protected:
    std::vector<double> K = std::vector<double>();
    bool first_time = true;

    inline virtual void insert_element_matrix(const std::vector<double>& k, const std::vector<long>& pos) override{
        const size_t w = pos.size();
        for(size_t i = 0; i < w; ++i){
            for(size_t j = i; j < w; ++j){
                if(pos[i] > -1 && pos[j] > -1){
                    this->K[utils::to_band(pos[i], pos[j], this->N)] += k[w*i + j];
                }
            }
        }
    }
    inline virtual void add_to_matrix(size_t i, size_t j, double val) override{
        if(j <= i){
            this->K[utils::to_band(i, j, this->N)] += val;
        } else {
            this->K[utils::to_band(j, i, this->N)] += val;
        }
    }
};


}

#endif
