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

#include <lapacke.h>
#include "utils/LIN3B.hpp"
#include "utils/gauss_legendre.hpp"

namespace utils{


LIN3B::LIN3B(std::array<double, 3> Y):Y(std::move(Y)){
    this->len = std::abs(Y[0] - Y[2]);

    constexpr size_t N = 3;
    std::array<double, N*N> M = 
        {1, Y[0], Y[0]*Y[0],
         1, Y[1], Y[1]*Y[1],
         1, Y[2], Y[2]*Y[2]};

    std::array<int, N> ipiv;

    // M*C = I -> C=M^-1
    // C = {a[0], a[1], a[2], a[3],
    //      b[0], b[1], b[2], b[3],
    //      c[0], c[1], c[2], c[3],
    //      d[0], d[1], d[2], d[3]}
    int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating LU in LIN3B.", info);
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating computing inverse from LU in LIN3B.", info);

    for(size_t i = 0; i < N; ++i){
        this->a[i] = M[i];
        this->b[i] = M[i+N];
        this->c[i] = M[i+2*N];
    }
}

Eigen::Matrix<double, 3, 3> LIN3B::absorption() const{
    Eigen::Matrix<double, 3, 3> result{{0,0,0},{0,0,0},{0,0,0}};
    const auto& gsi = GaussLegendre<3>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const double y = it->x*(Y[2] - Y[0]) + Y[0];
        auto NN = this->NN(y);

        result += it->w*NN*NN.transpose();
    }

    return this->len*result;
}

Eigen::Vector<double, 3> LIN3B::source(std::function<double(double)> fn) const{
    Eigen::Vector<double, 3> result{0,0,0};
    const auto& gsi = GaussLegendre<6>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const double y = it->x*(Y[2] - Y[0]) + Y[0];
        auto NN = this->NN(y);

        result += it->w*NN*fn(y);
    }

    return this->len*result;
}


}
