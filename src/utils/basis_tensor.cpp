/*
 *   Copyright (C) 2023 Tarcísio Ladeia de Oliveira.
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

#include "utils/basis_tensor.hpp"

namespace utils{

std::vector<double> basis_tensor_2D(const Eigen::Matrix<double, 2, 2>& R){
    return std::vector<double>
        {R(0,0)*R(0,0), R(0,1)*R(0,1), 2*R(0,0)*R(0,1),
         R(1,0)*R(1,0), R(1,1)*R(1,1), 2*R(1,0)*R(1,1),
         R(0,0)*R(1,0), R(0,1)*R(1,1),   R(0,0)*R(1,1)+R(0,1)*R(1,0)};
}

std::vector<double> basis_tensor_2D_inv_T(const Eigen::Matrix<double, 2, 2>& R){
    return std::vector<double>
        {  R(0,0)*R(0,0),   R(0,1)*R(0,1), R(0,0)*R(0,1),
           R(1,0)*R(1,0),   R(1,1)*R(1,1), R(1,0)*R(1,1),
         2*R(0,0)*R(1,0), 2*R(0,1)*R(1,1), R(0,0)*R(1,1)+R(0,1)*R(1,0)};
}

std::vector<double> basis_tensor_3D(const Eigen::Matrix<double, 3, 3>& R){
    std::array<double, 9> K1 =
        {R(0,0)*R(0,0), R(0,1)*R(0,1), R(0,2)*R(0,2),
         R(1,0)*R(1,0), R(1,1)*R(1,1), R(1,2)*R(1,2),
         R(2,0)*R(2,0), R(2,1)*R(2,1), R(2,2)*R(2,2)};

    std::array<double, 9> K2 =
        {R(0,0)*R(0,1), R(0,2)*R(0,0), R(0,1)*R(0,2),
         R(1,0)*R(1,1), R(1,2)*R(1,0), R(1,1)*R(1,2),
         R(2,0)*R(2,1), R(2,2)*R(2,0), R(2,1)*R(2,2)};

    std::array<double, 9> K3 =
        {R(0,0)*R(1,0), R(0,1)*R(1,1), R(0,2)*R(1,2),
         R(2,0)*R(0,0), R(2,1)*R(0,1), R(2,2)*R(0,2),
         R(1,0)*R(2,0), R(1,1)*R(2,1), R(1,2)*R(2,2)};

    std::array<double, 9> K4 =
        {R(0,0)*R(1,1)+R(0,1)*R(1,0), R(0,2)*R(1,0)+R(0,0)*R(1,2), R(0,1)*R(1,2)+R(0,2)*R(1,1),
         R(2,0)*R(0,1)+R(2,1)*R(0,0), R(2,2)*R(0,0)+R(2,0)*R(0,2), R(2,1)*R(0,2)+R(2,2)*R(0,1),
         R(1,0)*R(2,1)+R(1,1)*R(2,0), R(1,2)*R(2,0)+R(1,0)*R(2,2), R(1,1)*R(2,2)+R(1,2)*R(2,1)};

    std::vector<double> K(36, 0);
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            K[i*6 + j] = K1[i*3 + j];
            K[i*6 + (j+3)] = 2*K2[i*3 + j];
            K[(i+3)*6 + j] = K3[i*3 + j];
            K[(i+3)*6 + (j+3)] = K4[i*3 + j];
        }
    }

    return K;
}

std::vector<double> basis_tensor_3D_inv_T(const Eigen::Matrix<double, 3, 3>& R){
    std::array<double, 9> K1 =
        {R(0,0)*R(0,0), R(0,1)*R(0,1), R(0,2)*R(0,2),
         R(1,0)*R(1,0), R(1,1)*R(1,1), R(1,2)*R(1,2),
         R(2,0)*R(2,0), R(2,1)*R(2,1), R(2,2)*R(2,2)};

    std::array<double, 9> K2 =
        {R(0,0)*R(0,1), R(0,2)*R(0,0), R(0,1)*R(0,2),
         R(1,0)*R(1,1), R(1,2)*R(1,0), R(1,1)*R(1,2),
         R(2,0)*R(2,1), R(2,2)*R(2,0), R(2,1)*R(2,2)};

    std::array<double, 9> K3 =
        {R(0,0)*R(1,0), R(0,1)*R(1,1), R(0,2)*R(1,2),
         R(2,0)*R(0,0), R(2,1)*R(0,1), R(2,2)*R(0,2),
         R(1,0)*R(2,0), R(1,1)*R(2,1), R(1,2)*R(2,2)};

    std::array<double, 9> K4 =
        {R(0,0)*R(1,1)+R(0,1)*R(1,0), R(0,2)*R(1,0)+R(0,0)*R(1,2), R(0,1)*R(1,2)+R(0,2)*R(1,1),
         R(2,0)*R(0,1)+R(2,1)*R(0,0), R(2,2)*R(0,0)+R(2,0)*R(0,2), R(2,1)*R(0,2)+R(2,2)*R(0,1),
         R(1,0)*R(2,1)+R(1,1)*R(2,0), R(1,2)*R(2,0)+R(1,0)*R(2,2), R(1,1)*R(2,2)+R(1,2)*R(2,1)};

    std::vector<double> K(36, 0);
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            K[i*6 + j] = K1[i*3 + j];
            K[i*6 + (j+3)] = K2[i*3 + j];
            K[(i+3)*6 + j] = 2*K3[i*3 + j];
            K[(i+3)*6 + (j+3)] = K4[i*3 + j];
        }
    }

    return K;
}


}
