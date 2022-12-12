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

#ifndef LAPACK_DENSE_SYMMETRIC_BANDED
#define LAPACK_DENSE_SYMMETRIC_BANDED

#include "meshing.hpp"

namespace global_stiffness_matrix{

class LAPACKDenseSymmetricBanded{
    public:
    const double K_MIN = 1e-6;

    virtual ~LAPACKDenseSymmetricBanded() = default;

    virtual void calculate_dimensions(const Meshing * const mesh, const std::vector<double>& load);

    virtual void generate(const Meshing * const mesh, const std::vector<double>& density, const double pc);

    inline std::vector<double>& get_K(){
        return K;
    }
    inline const std::vector<double>& get_K() const{
        return K;
    }

    inline size_t get_W() const{
        return W;
    }
    inline size_t get_N() const{
        return N;
    }
    inline bool calculated_dimensions() const{
        return recalculated_dimensions;
    }

    protected:
    std::vector<double> K = std::vector<double>();
    size_t W = 0;
    size_t N = 0;
    bool recalculated_dimensions = false;

    virtual void add_geometry(const Meshing * const mesh, const Geometry * const g);

    virtual void add_geometry(const Meshing * const mesh, const Geometry * const g, std::vector<double>::const_iterator& rho, const double pc);

    inline void insert_element_matrix(const std::vector<double>& k, const std::vector<long>& pos, const size_t n){
        const size_t w = pos.size();
        for(size_t i = 0; i < w; ++i){
            for(size_t j = i; j < w; ++j){
                if(pos[i] > -1 && pos[j] > -1){
                    this->K[utils::to_band(pos[i], pos[j], n)] += k[w*i + j];
                }
            }
        }
    }
};


}

#endif
