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
#include "boundary_element/BTRI3.hpp"
#include "utils/gauss_legendre.hpp"

namespace boundary_element{

BTRI3::BTRI3(ElementShape s):
    BoundaryMeshElement(s.nodes){
    const size_t N = BTRI3::NODES_PER_ELEM;
    
    std::array<double, N> x, y, z;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
        z[i] = this->nodes[i]->point.Z();
    }
    std::array<double, N*N> M = 
        {x[0], y[0], z[0],
         x[1], y[1], z[1],
         x[2], y[2], z[2]};

    std::array<int, N> ipiv;

    // M*C = I -> C=M^-1
    // C = {a[0], a[1], a[2],
    //      b[0], b[1], b[2],
    //      c[0], c[1], c[2]}
    int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating LU in Q4.", info);
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating computing inverse from LU in Q4.", info);

    for(size_t i = 0; i < N; ++i){
        this->a[i] = M[i];
        this->b[i] = M[i+3];
        this->c[i] = M[i+6];
    }

    gp_Vec v1(this->nodes[0]->point, this->nodes[1]->point);
    gp_Vec v2(this->nodes[0]->point, this->nodes[2]->point);


    this->delta = 0.5*(v1.Crossed(v2)).Magnitude();
}

Eigen::MatrixXd BTRI3::diffusion_1dof(const Eigen::MatrixXd& A) const{
    const auto B = this->dN_mat_1dof();
    return delta*B.transpose()*A*B;
}
Eigen::MatrixXd BTRI3::advection_1dof(const Eigen::VectorXd& v) const{
    const auto NN = this->N_mat_1dof(this->GS_point(1.0/3.0, 1.0/3.0, 1.0/3.0));
    const auto B = this->dN_mat_1dof();
    return this->delta*NN*(v.transpose()*B);
}
Eigen::MatrixXd BTRI3::absorption_1dof() const{
    const auto& gsi = utils::GaussLegendreTri<2>::get();
    Eigen::MatrixXd result{{0, 0, 0},{0, 0, 0},{0, 0, 0}};
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const auto NN = this->N_mat_1dof(this->GS_point(it->a, it->b, it->c));
        result += NN.transpose()*NN;
    }

    return delta*result;
}
Eigen::VectorXd BTRI3::source_1dof() const{
    return this->delta*this->N_mat_1dof(this->GS_point(1.0/3.0, 1.0/3.0, 1.0/3.0));
}

}