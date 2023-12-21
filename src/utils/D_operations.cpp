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
#include <lapacke.h>
#include <cblas.h>

namespace utils::D_op{

std::vector<double> invert_2D(const std::vector<double>& d){
    constexpr size_t N = 3;
    std::vector<int> ipiv(N*N);
    std::vector<double> D = d;
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N, N, D.data(), N, ipiv.data());
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, N, D.data(), N, ipiv.data());
    return D;
}
std::vector<double> invert_3D(const std::vector<double>& d){
    constexpr size_t N = 6;
    std::vector<int> ipiv(N*N);
    std::vector<double> D = d;
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N, N, D.data(), N, ipiv.data());
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, N, D.data(), N, ipiv.data());
    return D;
}
std::vector<double> square_2D(const std::vector<double>& d){
    constexpr size_t N = 3;
    std::vector<double> D = d;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, d.data(), N, d.data(), N, 0, D.data(), N);
    return D;
}
std::vector<double> square_3D(const std::vector<double>& d){
    constexpr size_t N = 6;
    std::vector<double> D = d;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, d.data(), N, d.data(), N, 0, D.data(), N);
    return D;
}
std::vector<double> mult_2D(const std::vector<double>& d, const std::vector<double>& e){
    constexpr size_t N = 3;
    std::vector<double> D = d;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, d.data(), N, e.data(), N, 0, D.data(), N);
    return D;
}
std::vector<double> mult_3D(const std::vector<double>& d, const std::vector<double>& e){
    constexpr size_t N = 6;
    std::vector<double> D = d;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, d.data(), N, e.data(), N, 0, D.data(), N);
    return D;
}
std::vector<double> square_root_2D(const std::vector<double>& d){
    constexpr size_t N = 3;
    typedef Eigen::Matrix<double, N, N> MatType;
    MatType M = Eigen::Map<const MatType>(d.data(), N, N);

    Eigen::SelfAdjointEigenSolver<MatType> eigen(M);

    M = eigen.operatorSqrt();

    std::vector<double> result(N*N);
    std::copy(M.data(), M.data()+N*N, result.begin());

    return result;
}
std::vector<double> square_root_3D(const std::vector<double>& d){
    constexpr size_t N = 6;
    typedef Eigen::Matrix<double, N, N> MatType;
    MatType M = Eigen::Map<const MatType>(d.data(), N, N);

    Eigen::SelfAdjointEigenSolver<MatType> eigen(M);

    M = eigen.operatorSqrt();

    std::vector<double> result(N*N);
    std::copy(M.data(), M.data()+N*N, result.begin());

    return result;
}

}
