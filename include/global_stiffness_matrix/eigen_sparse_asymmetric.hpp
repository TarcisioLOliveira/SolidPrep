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

#ifndef EIGEN_SPARSE_ASYMMETRIC_HPP
#define EIGEN_SPARSE_ASYMMETRIC_HPP

#include <Eigen/SparseCore>
#include "global_stiffness_matrix.hpp"
#include "math/matrix.hpp"
#include "meshing.hpp"

namespace global_stiffness_matrix{

class EigenSparseAsymmetric : public GlobalStiffnessMatrix{
    public:
    typedef Eigen::SparseMatrix<double, Eigen::RowMajor, std::ptrdiff_t> Mat;

    EigenSparseAsymmetric(double EPS_DISPL):GlobalStiffnessMatrix(EPS_DISPL){}
    virtual ~EigenSparseAsymmetric() = default;

    virtual void generate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type) override;

    inline virtual void dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const override{
        Eigen::VectorXd u = Eigen::Map<const Eigen::VectorXd>(v.data(), v.size());

        Eigen::VectorXd u_out = this->K*u;
        std::copy(u_out.begin(), u_out.end(), v_out.begin());
    }
    inline virtual void reset_hessian() override{
       this->K = this->K_bkp; 
    };

    Mat& get_K() {
        return K;
    }

    protected:
    bool first_time = true;
    Mat K;
    Mat K_bkp;
    size_t u_size, l_num;

    inline virtual void insert_block_symmetric(const math::Matrix& k, const std::vector<long>& posi, const std::vector<long>& posj) override{
        const size_t w = posj.size();
        const size_t h = posi.size();
        for(size_t i = 0; i < h; ++i){
            if(posi[i] < 0){
                continue;
            }
            for(size_t j = 0; j < w; ++j){
                if(posj[j] < 0){
                    continue;
                }
                if(posi[i] > -1 && posj[j] > -1){
                    K.coeffRef(posi[i], posj[j]) += k(i, j);
                    K.coeffRef(posj[j], posi[i]) += k(i, j);
                }
            }
        }
    }

    inline virtual void insert_element_matrix(const math::Matrix& k, const std::vector<long>& pos) override{
        const size_t w = pos.size();
        for(size_t i = 0; i < w; ++i){
            if(pos[i] < 0){
                continue;
            }
            for(size_t j = 0; j < w; ++j){
                if(pos[j] < 0){
                    continue;
                }
                K.coeffRef(pos[i], pos[j]) += k(i, j);
            }
        }
    }
    inline virtual void add_to_matrix(size_t i, size_t j, double val) override{
        K.coeffRef(i, j) += val;
    }
};

}

#endif
