/*
 *   Copyright (C) 2024 Tarcísio Ladeia de Oliveira.
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

#include <algorithm>
#include <lapacke.h>
#include "math/matrix.hpp"
#include "logger.hpp"

namespace math{

Matrix::Matrix(size_t W, size_t H, Scalar s):
    W(W), H(H), M(new Scalar[W*H]){

    this->fill(s);
}

Matrix::Matrix(Scalar* m, size_t W, size_t H):
    W(W), H(H), M(new Scalar[W*H]){

    std::copy(m, m + W*H, this->M);

}
Matrix::Matrix(std::vector<Scalar> m, size_t W, size_t H):
    W(W), H(H), M(new Scalar[W*H]){

    std::copy(m.begin(), m.end(), this->M);
}
Matrix::Matrix(std::initializer_list<Scalar> m, size_t W, size_t H):
    W(W), H(H), M(new Scalar[W*H]){

    std::copy(m.begin(), m.end(), this->M);
}

void Matrix::fill(Scalar s){
    std::fill(M, M + W*H, s);
}

Matrix::~Matrix(){
    delete[] M;
}

Matrix::Matrix(const Matrix& m):
    W(m.W), H(m.H), M(new Scalar[W*H]){
    
    std::copy(m.M, m.M+W*H, this->M);
}

Matrix::Matrix(const MatrixTransposeView& m):
    W(m.H), H(m.W), M(new Scalar[W*H]){
    
    for(size_t i = 0; i < W; ++i){
        for(size_t j = 0; j < H; ++j){
            this->at(i,j) = m(i, j);
        }
    }
}

Matrix::Matrix(Matrix&& m):
    W(m.W), H(m.H), M(m.M){

    m.M = nullptr;
}

Matrix& Matrix::operator=(const Matrix& m){
    W = m.W;
    H = m.H;
    delete[] this->M;
    this->M = new Scalar[W*H];

    std::copy(m.M, m.M+W*H, this->M);

    return *this;
}
Matrix& Matrix::operator=(Matrix&& m){
    W = m.W;
    H = m.H;
    M = m.M;
    m.M = nullptr;

    return *this;
}

Matrix& Matrix::operator=(const MatrixTransposeView& m){
    W = m.H;
    H = m.W;
    delete[] this->M;
    this->M = new Scalar[W*H];

    for(size_t i = 0; i < W; ++i){
        for(size_t j = 0; j < H; ++j){
            this->at(i,j) = m(i, j);
        }
    }

    return *this;
}

Matrix Matrix::get_inverted() const{
    Matrix m(*this);
    m.invert();

    return m;
}
void Matrix::invert(){
    logger::log_assert(W == H,
                       logger::ERROR,
                       "unable to invert matrix, it is not square");

    const size_t N = W;

    // This works correctly even without the matrix being column major and
    // symmetric
    std::vector<int> ipiv(W);
    int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, this->M, N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating LU.", info);
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR, N, this->M, N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating computing inverse from LU.", info);
}

Matrix& Matrix::operator+=(const Matrix& m){
    logger::log_assert(W == m.W && H == m.H,
                       logger::ERROR,
                       "incompatible dimensions in matrix sum: ({}, {}) and ({}, {})",
                       H, W, m.H, m.W);

    for(size_t i = 0; i < W*H; ++i){
        this->M[i] += m.M[i];
    }

    return *this;
}
Matrix& Matrix::operator-=(const Matrix& m){
    logger::log_assert(W == m.W && H == m.H,
                       logger::ERROR,
                       "incompatible dimensions in matrix sum: ({}, {}) and ({}, {})",
                       H, W, m.H, m.W);

    for(size_t i = 0; i < W*H; ++i){
        this->M[i] -= m.M[i];
    }

    return *this;
}
Matrix& Matrix::operator+=(const MatrixTransposeView& m){
    logger::log_assert(W == m.W && H == m.H,
                       logger::ERROR,
                       "incompatible dimensions in matrix sum: ({}, {}) and ({}, {})",
                       H, W, m.H, m.W);

    for(size_t i = 0; i < W; ++i){
        for(size_t j = 0; j < H; ++j){
            this->at(i,j) += m(i, j);
        }
    }

    return *this;
}
Matrix& Matrix::operator-=(const MatrixTransposeView& m){
    logger::log_assert(W == m.W && H == m.H,
                       logger::ERROR,
                       "incompatible dimensions in matrix sum: ({}, {}) and ({}, {})",
                       H, W, m.H, m.W);

    for(size_t i = 0; i < W; ++i){
        for(size_t j = 0; j < H; ++j){
            this->at(i,j) -= m(i, j);
        }
    }

    return *this;
}

Matrix& Matrix::operator*=(Scalar s){
    for(Scalar* mi = M; mi < M + W*H; ++mi){
        *mi *= s;
    }
    return *this;
}
Matrix& Matrix::operator/=(Scalar s){
    for(Scalar* mi = M; mi < M + W*H; ++mi){
        *mi *= s;
    }
    return *this;
}


Matrix Matrix::operator+(const Matrix& m) const{
    Matrix m2(*this);
    m2 += m;

    return m2;
}
Matrix Matrix::operator-(const Matrix& m) const{
    Matrix m2(*this);
    m2 -= m;

    return m2;
}

Matrix Matrix::operator*(const Matrix& m) const{
    logger::log_assert(W == m.H,
                       logger::ERROR,
                       "incompatible dimensions in matrix sum: ({}, {}) and ({}, {})",
                       H, W, m.H, m.W);

    Matrix r(H, m.W);

    for(size_t i = 0; i < r.H; ++i){
        for(size_t k = 0; k < this->W; ++k){
            const double tmp = this->at(i,k);
            for(size_t j = 0; j < r.W; ++j){
                r(i,j) += tmp*m(k,j);
            }
        }
    }

    return r;
}

Matrix Matrix::operator*(Scalar s) const{
    Matrix m(*this);
    m *= s;

    return m;
}
Matrix Matrix::operator/(Scalar s) const{
    Matrix m(*this);
    m /= s;

    return m;
}

bool Matrix::operator==(const Matrix& m) const{
    if(W != m.W) return false;
    if(H != m.H) return false;
    for(size_t i = 0; i < W*H; ++i){
        if(M[i] != m.M[i]) return false;
    }
    return true;
}

bool Matrix::operator!=(const Matrix& m) const{
    return !(*this == m);
}

Matrix Matrix::operator-() const{
    return -1*(*this);
}

Matrix operator*(Scalar s, const Matrix& m){
    return m*s;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////// MATRIX TRANSPOSE VIEW ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

MatrixTransposeView::MatrixTransposeView(size_t non_T_W, size_t non_T_H, Scalar* M):
    W(non_T_H), H(non_T_W), M(M){

}

Matrix MatrixTransposeView::operator+(const Matrix& m) const{
    Matrix m2(m);
    m2 += *this;

    return m2;
}
Matrix MatrixTransposeView::operator-(const Matrix& m) const{
    const size_t mW = m.get_W();
    const size_t mH = m.get_H();
    logger::log_assert(W == mW && H == mH,
                       logger::ERROR,
                       "incompatible dimensions in matrix sum: ({}, {}) and ({}, {})",
                       H, W, mH, mW);
    Matrix m2(W, H);
    for(size_t i = 0; i < W; ++i){
        for(size_t j = 0; j < H; ++j){
            m2(i,j) = this->at(i,j) -  m(i, j);
        }
    }

    return m2;
}
Matrix MatrixTransposeView::operator*(const Matrix& m) const{
    const size_t mW = m.get_W();
    const size_t mH = m.get_H();
    logger::log_assert(W == mH,
                       logger::ERROR,
                       "incompatible dimensions in matrix sum: ({}, {}) and ({}, {})",
                       H, W, mH, mW);

    Matrix r(H, mW);

    for(size_t i = 0; i < H; ++i){
        for(size_t k = 0; k < this->W; ++k){
            const double tmp = this->at(i,k);
            for(size_t j = 0; j < mW; ++j){
                r(i,j) += tmp*m(k,j);
            }
        }
    }

    return r;
}

Matrix MatrixTransposeView::operator*(Scalar s) const{
    Matrix m(*this);
    m *= s;

    return m;
}
Matrix MatrixTransposeView::operator/(Scalar s) const{
    Matrix m(*this);
    m /= s;

    return m;
}

}
