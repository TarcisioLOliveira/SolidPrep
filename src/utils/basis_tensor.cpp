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

math::Matrix basis_tensor_2D(const math::Matrix& R){
    return math::Matrix(
        {R(0,0)*R(0,0), R(0,1)*R(0,1), 2*R(0,0)*R(0,1),
         R(1,0)*R(1,0), R(1,1)*R(1,1), 2*R(1,0)*R(1,1),
         R(0,0)*R(1,0), R(0,1)*R(1,1),   R(0,0)*R(1,1)+R(0,1)*R(1,0)}, 3, 3);
}

math::Matrix basis_tensor_2D_inv_T(const math::Matrix& R){
    return math::Matrix(
        {  R(0,0)*R(0,0),   R(0,1)*R(0,1), R(0,0)*R(0,1),
           R(1,0)*R(1,0),   R(1,1)*R(1,1), R(1,0)*R(1,1),
         2*R(0,0)*R(1,0), 2*R(0,1)*R(1,1), R(0,0)*R(1,1)+R(0,1)*R(1,0)}, 3, 3);
}

math::Matrix basis_tensor_3D(const math::Matrix& R){
    math::Matrix K1(
        {R(0,0)*R(0,0), R(0,1)*R(0,1), R(0,2)*R(0,2),
         R(1,0)*R(1,0), R(1,1)*R(1,1), R(1,2)*R(1,2),
         R(2,0)*R(2,0), R(2,1)*R(2,1), R(2,2)*R(2,2)}, 3, 3);

    math::Matrix K2(
        {R(0,0)*R(0,1), R(0,2)*R(0,0), R(0,1)*R(0,2),
         R(1,0)*R(1,1), R(1,2)*R(1,0), R(1,1)*R(1,2),
         R(2,0)*R(2,1), R(2,2)*R(2,0), R(2,1)*R(2,2)}, 3, 3);

    math::Matrix K3(
        {R(0,0)*R(1,0), R(0,1)*R(1,1), R(0,2)*R(1,2),
         R(2,0)*R(0,0), R(2,1)*R(0,1), R(2,2)*R(0,2),
         R(1,0)*R(2,0), R(1,1)*R(2,1), R(1,2)*R(2,2)}, 3, 3);

    math::Matrix K4(
        {R(0,0)*R(1,1)+R(0,1)*R(1,0), R(0,2)*R(1,0)+R(0,0)*R(1,2), R(0,1)*R(1,2)+R(0,2)*R(1,1),
         R(2,0)*R(0,1)+R(2,1)*R(0,0), R(2,2)*R(0,0)+R(2,0)*R(0,2), R(2,1)*R(0,2)+R(2,2)*R(0,1),
         R(1,0)*R(2,1)+R(1,1)*R(2,0), R(1,2)*R(2,0)+R(1,0)*R(2,2), R(1,1)*R(2,2)+R(1,2)*R(2,1)}, 3, 3);

    math::Matrix K(6, 6);
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            K(i, j) = K1(i, j);
            K(i, (j+3)) = 2*K2(i, j);
            K((i+3), j) = K3(i, j);
            K((i+3), (j+3)) = K4(i, j);
        }
    }

    return K;
}

math::Matrix basis_tensor_3D_inv_T(const math::Matrix& R){
    math::Matrix K1(
        {R(0,0)*R(0,0), R(0,1)*R(0,1), R(0,2)*R(0,2),
         R(1,0)*R(1,0), R(1,1)*R(1,1), R(1,2)*R(1,2),
         R(2,0)*R(2,0), R(2,1)*R(2,1), R(2,2)*R(2,2)}, 3, 3);

    math::Matrix K2(
        {R(0,0)*R(0,1), R(0,2)*R(0,0), R(0,1)*R(0,2),
         R(1,0)*R(1,1), R(1,2)*R(1,0), R(1,1)*R(1,2),
         R(2,0)*R(2,1), R(2,2)*R(2,0), R(2,1)*R(2,2)}, 3, 3);

    math::Matrix K3(
        {R(0,0)*R(1,0), R(0,1)*R(1,1), R(0,2)*R(1,2),
         R(2,0)*R(0,0), R(2,1)*R(0,1), R(2,2)*R(0,2),
         R(1,0)*R(2,0), R(1,1)*R(2,1), R(1,2)*R(2,2)}, 3, 3);

    math::Matrix K4(
        {R(0,0)*R(1,1)+R(0,1)*R(1,0), R(0,2)*R(1,0)+R(0,0)*R(1,2), R(0,1)*R(1,2)+R(0,2)*R(1,1),
         R(2,0)*R(0,1)+R(2,1)*R(0,0), R(2,2)*R(0,0)+R(2,0)*R(0,2), R(2,1)*R(0,2)+R(2,2)*R(0,1),
         R(1,0)*R(2,1)+R(1,1)*R(2,0), R(1,2)*R(2,0)+R(1,0)*R(2,2), R(1,1)*R(2,2)+R(1,2)*R(2,1)}, 3, 3);

    math::Matrix K(36, 0);
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            K(i, j) = K1(i, j);
            K(i, (j+3)) = K2(i, j);
            K((i+3), j) = 2*K3(i, j);
            K((i+3), (j+3)) = K4(i, j);
        }
    }

    return K;
}


}
