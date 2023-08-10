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

#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "utils/D_operations.hpp"

namespace utils::D_op{

std::vector<double> invert_2D(const std::vector<double>& d){
    std::vector<double> D = {
    d[4]/(d[0]*d[4] - d[1]*d[3])
    ,
    -d[1]/(d[0]*d[4] - d[1]*d[3])
    ,
    0
    ,
    -d[3]/(d[0]*d[4] - d[1]*d[3])
    ,
    d[0]/(d[0]*d[4] - d[1]*d[3])
    ,
    0
    ,
    0
    ,
    0
    ,
    1/d[8]
    };
    return D;
}
std::vector<double> invert_3D(const std::vector<double>& d){
    std::vector<double> D = {
    (d[13]*d[8] - d[14]*d[7])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    (d[1]*d[14] - d[13]*d[2])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    (-d[1]*d[8] + d[2]*d[7])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    0
    ,
    0
    ,
    0
    ,
    (-d[12]*d[8] + d[14]*d[6])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    (-d[0]*d[14] + d[12]*d[2])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    (d[0]*d[8] - d[2]*d[6])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    0
    ,
    0
    ,
    0
    ,
    (d[12]*d[7] - d[13]*d[6])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    (d[0]*d[13] - d[1]*d[12])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    (-d[0]*d[7] + d[1]*d[6])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    1/d[21]
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    1/d[28]
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    1/d[35]
    };

    return D;
}
std::vector<double> square_2D(const std::vector<double>& d){
    std::vector<double> D = {
    d[0]*d[0] + d[1]*d[3]
    ,
    d[1]*(d[0] + d[4])
    ,
    0
    ,
    d[3]*(d[0] + d[4])
    ,
    d[1]*d[3] + d[4]*d[4]
    ,
    0
    ,
    0
    ,
    0
    ,
    d[8]*d[8]
    };
    return D;
}
std::vector<double> square_3D(const std::vector<double>& d){
    std::vector<double> D = {
    d[0]*d[0] + d[1]*d[6] + d[12]*d[2]
    ,
    d[0]*d[1] + d[1]*d[7] + d[13]*d[2]
    ,
    d[0]*d[2] + d[1]*d[8] + d[14]*d[2]
    ,
    0
    ,
    0
    ,
    0
    ,
    d[0]*d[6] + d[12]*d[8] + d[6]*d[7]
    ,
    d[1]*d[6] + d[13]*d[8] + d[7]*d[7]
    ,
    d[14]*d[8] + d[2]*d[6] + d[7]*d[8]
    ,
    0
    ,
    0
    ,
    0
    ,
    d[0]*d[12] + d[12]*d[14] + d[13]*d[6]
    ,
    d[1]*d[12] + d[13]*d[14] + d[13]*d[7]
    ,
    d[12]*d[2] + d[13]*d[8] + d[14]*d[14]
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    d[21]*d[21]
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    d[28]*d[28]
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    d[35]*d[35]
    };

    return D;
}
std::vector<double> mult_2D(const std::vector<double>& d, const std::vector<double>& e){
    std::vector<double> D = {
    d[0]*e[0] + d[1]*e[3]
    ,
    d[0]*e[1] + d[1]*e[4]
    ,
    0
    ,
    d[3]*e[0] + d[4]*e[3]
    ,
    d[3]*e[1] + d[4]*e[4]
    ,
    0
    ,
    0
    ,
    0
    ,
    d[8]*e[8]
    };

    return D;
}
std::vector<double> mult_3D(const std::vector<double>& d, const std::vector<double>& e){
    std::vector<double> D = {
    d[0]*e[0] + d[1]*e[6] + d[2]*e[12]
    ,
    d[0]*e[1] + d[1]*e[7] + d[2]*e[13]
    ,
    d[0]*e[2] + d[1]*e[8] + d[2]*e[14]
    ,
    0
    ,
    0
    ,
    0
    ,
    d[6]*e[0] + d[7]*e[6] + d[8]*e[12]
    ,
    d[6]*e[1] + d[7]*e[7] + d[8]*e[13]
    ,
    d[6]*e[2] + d[7]*e[8] + d[8]*e[14]
    ,
    0
    ,
    0
    ,
    0
    ,
    d[12]*e[0] + d[13]*e[6] + d[14]*e[12]
    ,
    d[12]*e[1] + d[13]*e[7] + d[14]*e[13]
    ,
    d[12]*e[2] + d[13]*e[8] + d[14]*e[14]
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    d[21]*e[21]
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    d[28]*e[28]
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    d[35]*e[35]
    };

    return D;
}
std::vector<double> square_root_2D(const std::vector<double>& d){
    std::vector<double> D = {
    std::sqrt(2)*d[1]*d[3]*((-d[0] + d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*std::sqrt(d[0] + d[4] - std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4])) + (d[0] - d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*std::sqrt(d[0] + d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4])))/((-d[0] + d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*(d[0] - d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))
    ,
    std::sqrt(2)*d[1]*(2*d[1]*d[3]*(d[0] - d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*std::sqrt(d[0] + d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4])) - (-d[0] + d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*std::sqrt(d[0] + d[4] - std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*(d[0]*d[0] - 2*d[0]*d[4] + d[0]*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]) + 2*d[1]*d[3] + d[4]*d[4] - d[4]*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4])))/((-d[0] + d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*(d[0] - d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*(d[0]*d[0] - 2*d[0]*d[4] + d[0]*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]) + 4*d[1]*d[3] + d[4]*d[4] - d[4]*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4])))
    ,
    0
    ,
    std::sqrt(2)*d[3]*(-std::sqrt(d[0] + d[4] - std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4])) + std::sqrt(d[0] + d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4])))/(2*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))
    ,
    std::sqrt(2)*(d[1]*d[3]*std::sqrt(d[0] + d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4])) + std::sqrt(d[0] + d[4] - std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*(d[0]*d[0] - 2*d[0]*d[4] + d[0]*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]) + 2*d[1]*d[3] + d[4]*d[4] - d[4]*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))/2)/(d[0]*d[0] - 2*d[0]*d[4] + d[0]*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]) + 4*d[1]*d[3] + d[4]*d[4] - d[4]*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] +d[4]*d[4]))
    ,
    0
    ,
    0
    ,
    0
    ,
    std::sqrt(d[8])
    };

    return D;
}
std::vector<double> square_root_3D(const std::vector<double>& d){
    typedef Eigen::Matrix<double, 6, 6> MatType;
    MatType M = Eigen::Map<const MatType>(d.data(), 6, 6);

    Eigen::SelfAdjointEigenSolver<MatType> eigen(M);

    M = eigen.operatorSqrt();

    std::vector<double> result(36);
    std::copy(M.data(), M.data()+36, result.begin());

    return result;
}

}
