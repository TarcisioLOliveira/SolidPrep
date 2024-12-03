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
#include <utility>
#include "math/matrix.hpp"
#include "logger.hpp"

namespace math{

///////////////////////////////////////////////////////////////////////////////
//////////////////////////////// MATRIX ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

Matrix::Matrix(size_t H, size_t W, Scalar s):
    H(H), W(W), M(new Scalar[W*H]){

    this->fill(s);
}

Matrix::Matrix(Scalar* m, size_t H, size_t W):
    H(H), W(W), M(new Scalar[W*H]){

    std::copy(m, m + W*H, this->M);

}
Matrix::Matrix(std::vector<Scalar> m, size_t H, size_t W):
    H(H), W(W), M(new Scalar[W*H]){

    logger::log_assert(W*H == m.size(),
            logger::ERROR,
            "incorrect dimensions for Matrix");

    std::copy(m.begin(), m.end(), this->M);
}
Matrix::Matrix(std::initializer_list<Scalar> m, size_t H, size_t W):
    H(H), W(W), M(new Scalar[W*H]){

    logger::log_assert(W*H == m.size(),
            logger::ERROR,
            "incorrect dimensions for Matrix");

    std::copy(m.begin(), m.end(), this->M);
}

void Matrix::fill(Scalar s){
    std::fill(M, M + W*H, s);
}

Matrix::~Matrix(){
    delete[] M;
}

Matrix::Matrix(const Matrix& m):
    H(m.H), W(m.W), M(new Scalar[W*H]){
    
    std::copy(m.M, m.M+W*H, this->M);
}

Matrix::Matrix(const MatrixTransposeView& m):
    H(m.H), W(m.W), M(new Scalar[W*H]){
    
    for(size_t i = 0; i < H; ++i){
        for(size_t j = 0; j < W; ++j){
            this->at(i,j) = m(i, j);
        }
    }
}

Matrix::Matrix(Matrix&& m):
    H(m.H), W(m.W), M(m.M){

    m.M = nullptr;
}


Matrix Matrix::identity(size_t N){
    Matrix I(N, N);
    for(size_t i = 0; i < N; ++i){
        I(i, i) = 1;
    }
    return I;
}

Matrix Matrix::diag(const Vector& v){
    const size_t N = v.get_N();
    Matrix d(N, N);
    for(size_t i = 0; i < N; ++i){
        d(i, i) = v[i];
    }

    return d;
}

bool Matrix::is_equal(const Matrix& m, Scalar eps) const{
    if(W != m.W) return false;
    if(H != m.H) return false;
    for(size_t i = 0; i < W*H; ++i){
        if(std::abs(M[i] - m.M[i]) >= eps) return false;
    }
    return true;
}
bool Matrix::is_equal(const MatrixTransposeView& m, Scalar eps) const{
    if(W != m.W) return false;
    if(H != m.H) return false;
    for(size_t i = 0; i < H; ++i){
        for(size_t j = 0; j < W; ++j){
            if(std::abs(this->at(i,j) - m(i,j)) >= eps) return false;
        }
    }
    return true;
}

Scalar Matrix::determinant() const{
    logger::log_assert(W == H,
                       logger::ERROR,
                       "unable to calculate determinant, matrix is not square");
    std::vector<double> tmp(M, M + W*H);
    const size_t N = W;

    std::vector<int> ipiv(W);
    int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, tmp.data(), N, ipiv.data());
    logger::log_assert(info >= 0, logger::ERROR, "LAPACKE returned {} while calculating LU.", info);

    if(info > 0) return 0;

    Scalar det = 1;
    for(size_t i = 0; i < N; ++i){
        det *= tmp[i*N + i];
    }
    for(int i = 0; static_cast<size_t>(i) < N; i++){
        if(i+1 != ipiv[i]){
            det *= -1;
        }
    }
    
    return det;
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
    W = m.W;
    H = m.H;
    delete[] this->M;
    this->M = new Scalar[W*H];

    for(size_t i = 0; i < H; ++i){
        for(size_t j = 0; j < W; ++j){
            this->at(i,j) = m(i, j);
        }
    }

    return *this;
}

Matrix Matrix::get_inverted_LU() const{
    Matrix m(*this);
    m.invert_LU();

    return m;
}
void Matrix::invert_LU(){
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

Matrix Matrix::get_inverted_cholesky() const{
    Matrix m(*this);
    m.invert_cholesky();

    return m;
}
void Matrix::invert_cholesky(){
    logger::log_assert(W == H,
                       logger::ERROR,
                       "unable to invert matrix, it is not square");

    const size_t N = W;

    // This works correctly even without the matrix being column major and
    // symmetric
    int info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', N, this->M, N);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating Cholesky.", info);
    info = LAPACKE_dpotri(LAPACK_COL_MAJOR, 'L', N, this->M, N);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating computing inverse from Cholesky.", info);

    for(size_t i = 0; i < N; ++i){
        for(size_t j = i+1; j < N; ++j){
            M[j*N + i] = M[i*N + j];
        }
    }
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
                       "incompatible dimensions in matrix difference: ({}, {}) and ({}, {})",
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

    for(size_t i = 0; i < H; ++i){
        for(size_t j = 0; j < W; ++j){
            this->at(i,j) += m(i, j);
        }
    }

    return *this;
}
Matrix& Matrix::operator-=(const MatrixTransposeView& m){
    logger::log_assert(W == m.W && H == m.H,
                       logger::ERROR,
                       "incompatible dimensions in matrix difference: ({}, {}) and ({}, {})",
                       H, W, m.H, m.W);

    for(size_t i = 0; i < H; ++i){
        for(size_t j = 0; j < W; ++j){
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
        *mi /= s;
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
                       "incompatible dimensions in matrix multiplication: ({}, {}) and ({}, {})",
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
Matrix Matrix::operator*(const MatrixTransposeView& m) const{
    logger::log_assert(W == m.H,
                       logger::ERROR,
                       "incompatible dimensions in matrix multiplication: ({}, {}) and ({}, {})",
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

Vector Matrix::operator*(const VectorNonTransposed& v) const{
    logger::log_assert(W == v.get_N(), logger::ERROR,
            "matrix and vector have incompatible dimensions: ({}, {}), ({}, 1)",
            H, W, v.get_N());
    Vector v2(H);
    for(size_t j = 0; j < W; ++j){
        const double tmp = v[j];
        for(size_t i = 0; i < H; ++i){
            v2[i] += this->at(i,j)*tmp;
        }
    }

    return v2;
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

std::ostream& operator<<(std::ostream& output, const Matrix& m){
    for(size_t i = 0; i < m.get_H(); ++i){
        for(size_t j = 0; j < m.get_W(); ++j){
            output << m(i,j) << " ";
        }
        output << std::endl;
    }
    return output;
}


MatrixTransposeView Matrix::T() const{
    return MatrixTransposeView(H, W, M);
}

bool Matrix::operator==(const MatrixTransposeView&& m) const{
    if(W != m.W) return false;
    if(H != m.H) return false;

    for(size_t i = 0; i < H; ++i){
        for(size_t j = 0; j < W; ++j){
            if(this->at(i,j) != m(i, j)){
                return false;
            }
        }
    }
    return true;
}
bool Matrix::operator!=(const MatrixTransposeView&& m) const{
    return !(*this == m);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////// MATRIX TRANSPOSE VIEW ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

MatrixTransposeView::MatrixTransposeView(const size_t non_T_H, const size_t non_T_W, const Scalar* const M):
    H(non_T_W), W(non_T_H), M(M){

}

bool MatrixTransposeView::is_equal(const MatrixTransposeView& m, Scalar eps) const{
    if(W != m.W) return false;
    if(H != m.H) return false;
    for(size_t i = 0; i < W*H; ++i){
        if(std::abs(M[i] - m.M[i]) >= eps) return false;
    }
    return true;
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
                       "incompatible dimensions in matrix difference: ({}, {}) and ({}, {})",
                       H, W, mH, mW);
    Matrix m2(H, W);
    for(size_t i = 0; i < H; ++i){
        for(size_t j = 0; j < W; ++j){
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
                       "incompatible dimensions in matrix multiplication: ({}, {}) and ({}, {})",
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

Vector MatrixTransposeView::operator*(const VectorNonTransposed& v) const{
    logger::log_assert(W == v.get_N(), logger::ERROR,
            "matrix and vector have incompatible dimensions: ({}, {}), ({}, 1)",
            H, W, v.get_N());
    Vector v2(H);
    for(size_t j = 0; j < W; ++j){
        const double tmp = v[j];
        for(size_t i = 0; i < H; ++i){
            v2[i] += this->at(i,j)*tmp;
        }
    }

    return v2;
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

bool MatrixTransposeView::operator==(const MatrixTransposeView&& m) const{
    if(W != m.W) return false;
    if(H != m.H) return false;
    for(size_t i = 0; i < W*H; ++i){
        if(M[i] != m.M[i]) return false;
    }
    return true;
}
bool MatrixTransposeView::operator!=(const MatrixTransposeView&& m) const{
    return !(*this == m);
}

std::ostream& operator<<(std::ostream& output, const MatrixTransposeView& m){
    for(size_t i = 0; i < m.get_H(); ++i){
        for(size_t j = 0; j < m.get_W(); ++j){
            output << m(i,j) << " ";
        }
        output << std::endl;
    }
    return output;
}

///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// LU ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

LU::LU(Matrix& m, bool copy):
    N((m.get_W() == m.get_H()) ? m.get_W() : 0),
    M((copy) ? ((N > 0) ? new Scalar[N] : nullptr) : m.data()),
    ipiv(N, 0),
    copy(copy){

    logger::log_assert(N > 0 && M != nullptr, logger::ERROR,
        "unable to factorize matrix, not square");

    int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N, N, M, N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating LU.", info);
}

LU::~LU(){
    if(copy) delete[] M;
}

Scalar LU::determinant() const{
    Scalar det = 1;
    for(size_t i = 0; i < N; ++i){
        det *= M[i*N + i];
    }
    for(int i = 0; static_cast<size_t>(i) < N; i++){
        if(i+1 != ipiv[i]){
            det *= -1;
        }
    }
    
    return det;
}

void LU::solve(Vector& v) const{
    logger::log_assert(v.get_N() == N,
            logger::ERROR,
            "unable to solve problem, different dimensions: {} (matrix) and {} (vector)",
            N, v.get_N());
    LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', N, 1, M, N, ipiv.data(), v.data(), 1);
}
void LU::solve(Matrix& m) const{
    logger::log_assert(m.get_H() == N,
            logger::ERROR,
            "unable to solve problem, different dimensions: {} (LU matrix) and ({}, {})  (RHS)",
            N, m.get_H(), m.get_W());
    const size_t W = m.get_W();
    LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', N, m.get_W(), M, N, ipiv.data(), m.data(), W);
}
void LU::solve_transposed(Vector& v) const{
    logger::log_assert(v.get_N() == N,
            logger::ERROR,
            "unable to solve_transposed problem, different dimensions: {} (matrix) and {} (vector)",
            N, v.get_N());
    LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'T', N, 1, M, N, ipiv.data(), v.data(), 1);
}
void LU::solve_transposed(Matrix& m) const{
    logger::log_assert(m.get_H() == N,
            logger::ERROR,
            "unable to solve_transposed problem, different dimensions: {} (LU matrix) and ({}, {})  (RHS)",
            N, m.get_H(), m.get_W());
    const size_t W = m.get_W();
    LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'T', N, m.get_W(), M, N, ipiv.data(), m.data(), W);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Cholesky ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

Cholesky::Cholesky(Matrix& m, bool copy):
    N((m.get_W() == m.get_H()) ? m.get_W() : 0),
    M((copy) ? ((N > 0) ? new Scalar[N] : nullptr) : m.data()),
    ipiv(N, 0),
    copy(copy){

    logger::log_assert(N > 0 && M != nullptr, logger::ERROR,
        "unable to factorize matrix, not square");

    int info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', N, this->M, N);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating Cholesky.", info);
}

Cholesky::~Cholesky(){
    if(copy) delete[] M;
}

Scalar Cholesky::determinant() const{
    Scalar det = 1;
    for(size_t i = 0; i < N; ++i){
        det *= M[i*N + i];
    }
    for(int i = 0; static_cast<size_t>(i) < N; i++){
        if(i+1 != ipiv[i]){
            det *= -1;
        }
    }
    
    return det*det;
}

void Cholesky::solve(Vector& v) const{
    logger::log_assert(v.get_N() == N,
            logger::ERROR,
            "unable to solve problem, different dimensions: {} (matrix) and {} (vector)",
            N, v.get_N());
    LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'L', N, 1, M, N, v.data(), 1);
}
void Cholesky::solve(Matrix& m) const{
    logger::log_assert(m.get_H() == N,
            logger::ERROR,
            "unable to solve problem, different dimensions: {} (Cholesky matrix) and ({}, {})  (RHS)",
            N, m.get_H(), m.get_W());
    const size_t W = m.get_W();
    LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'L', N, W, M, N, m.data(), 1);
}

///////////////////////////////////////////////////////////////////////////////
////////////////////////////// VECTOR AGNOSTIC ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

VectorAgnostic::VectorAgnostic(const Vector& v):
        N(v.get_N()), V(v.data()){}
VectorAgnostic::VectorAgnostic(const VectorView& v):
        N(v.get_N()), V(v.data()){}
VectorAgnostic::VectorAgnostic(const VectorTranspose& v):
        N(v.get_N()), V(v.data()){}
VectorAgnostic::VectorAgnostic(const VectorTransposeView& v):
        N(v.get_N()), V(v.data()){}
VectorAgnostic::VectorAgnostic(const VectorTransposed& v):
        N(v.get_N()), V(v.V){}
VectorAgnostic::VectorAgnostic(const VectorNonTransposed& v):
        N(v.get_N()), V(v.V){}

VectorNonTransposed::VectorNonTransposed(const Vector& v):
        N(v.get_N()), V(v.data()){}
VectorNonTransposed::VectorNonTransposed(const VectorView& v):
        N(v.get_N()), V(v.data()){}

VectorTransposed::VectorTransposed(const VectorTranspose& v):
        N(v.get_N()), V(v.data()){}
VectorTransposed::VectorTransposed(const VectorTransposeView& v):
        N(v.get_N()), V(v.data()){}

///////////////////////////////////////////////////////////////////////////////
//////////////////////////////// VECTOR ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

Vector::Vector(Scalar* v, size_t N):
    VectorMut(N, new Scalar[N]){

    std::copy(v, v + N, this->V);
}
Vector::Vector(std::vector<Scalar> v):
    VectorMut(v.size(), new Scalar[v.size()]){

    std::copy(v.begin(), v.begin() + N, this->V);
}
Vector::Vector(std::initializer_list<Scalar> v):
    VectorMut(v.size(), new Scalar[v.size()]){

    std::copy(v.begin(), v.begin() + N, this->V);
}
Vector::Vector(size_t N, Scalar s):
    VectorMut(N, new Scalar[N]){

    this->fill(s);
}
Vector::~Vector(){
    delete[] this->V;
}
Vector::Vector(const VectorAgnostic& v):
    VectorMut(v.N, new Scalar[v.N]){

    std::copy(v.V, v.V + N, this->V);
}
Vector::Vector(const Vector& v):
    VectorMut(v.N, new Scalar[v.N]){

    std::copy(v.V, v.V + N, this->V);
}
Vector::Vector(Vector&& v):
    VectorMut(v.N, v.V){
    
    v.V = nullptr;
}
Vector::Vector(const VectorMut<Vector>& v):
    VectorMut(v){
}
Vector::Vector(VectorMut<Vector>&& v):
    VectorMut(std::forward<VectorMut<Vector>>(v)){
}


VectorTransposeView Vector::T() const{
    return VectorTransposeView(N, V);
}

Vector Vector::operator*(Scalar s) const{
    Vector v(*this);
    v *= s;

    return v;
}
Vector Vector::operator/(Scalar s) const{
    Vector v(*this);
    v /= s;

    return v;
}

Vector Vector::operator+(const VectorNonTransposed& v) const{
    Vector v2(*this);
    v2 += v;

    return v2;
}
Vector Vector::operator-(const VectorNonTransposed& v) const{
    Vector v2(*this);
    v2 -= v;

    return v2;
}

Vector& Vector::operator*=(Scalar s){
    for(size_t i = 0; i < N; ++i){
        V[i] *= s;
    }
    return *this;
}
Vector& Vector::operator/=(Scalar s){
    for(size_t i = 0; i < N; ++i){
        V[i] /= s;
    }
    return *this;
}

Vector& Vector::operator+=(const VectorNonTransposed& v){
    logger::log_assert(N == v.N, logger::ERROR,
            "vectors have different sizes: {}, {}",
            N, v.N);

    for(size_t i = 0; i < N; ++i){
        V[i] += v[i];
    }
    return *this;
}
Vector& Vector::operator-=(const VectorNonTransposed& v){
    logger::log_assert(N == v.N, logger::ERROR,
            "vectors have different sizes: {}, {}",
            N, v.N);
    for(size_t i = 0; i < N; ++i){
        V[i] -= v[i];
    }
    return *this;
}

Matrix Vector::operator*(const VectorTransposed& v) const{
    Matrix m(N, v.N);
    for(size_t i = 0; i < N; ++i){
        double tmp = V[i];
        for(size_t j = 0; j < v.N; ++j){
            m(i, j) = tmp*v[j];
        }
    }

    return m;
}

///////////////////////////////////////////////////////////////////////////////
//////////////////////////// VECTOR TRANSPOSE /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

VectorTranspose::VectorTranspose(const VectorMut<VectorTranspose>& v):
    VectorMut(v){}
VectorTranspose::VectorTranspose(VectorMut<VectorTranspose>&& v):
    VectorMut(std::forward<VectorMut<VectorTranspose>>(v)){}

VectorTranspose::VectorTranspose(const VectorAgnostic& v):
    VectorMut(v.N, new Scalar[v.N]){

    std::copy(v.V, v.V + N, this->V);
}
VectorTranspose::VectorTranspose(const VectorTranspose& v):
    VectorMut(v.N, new Scalar[v.N]){

    std::copy(v.V, v.V + N, this->V);
}
VectorTranspose::VectorTranspose(VectorTranspose&& v):
    VectorMut(v.N, v.V){
    
    v.V = nullptr;
}
VectorTranspose::VectorTranspose(size_t N, Scalar s):
    VectorMut(N, new Scalar[N]){

    this->fill(s);
}
VectorTranspose::~VectorTranspose(){
    delete[] this->V;
}
VectorTranspose& VectorTranspose::operator=(const VectorTranspose& v){
    this->N = v.N;
    delete[] this->V;
    this->V = new Scalar[N];

    std::copy(v.V, v.V + N, this->V);

    return *this;
}
VectorTranspose& VectorTranspose::operator=(VectorTranspose&& v){
    this->N = v.N;
    delete[] this->V;
    this->V = v.V;

    v.V = nullptr;

    return *this;
}
VectorTranspose VectorTranspose::operator*(Scalar s) const{
    VectorTranspose v(*this);
    v *= s;

    return v;
}
VectorTranspose VectorTranspose::operator/(Scalar s) const{
    VectorTranspose v(*this);
    v /= s;

    return v;
}

VectorTranspose VectorTranspose::operator+(const VectorTransposed& v) const{
    VectorTranspose v2(*this);
    v2 += v;

    return v2;
}
VectorTranspose VectorTranspose::operator-(const VectorTransposed& v) const{
    VectorTranspose v2(*this);
    v2 -= v;

    return v2;
}

VectorTranspose& VectorTranspose::operator+=(const VectorTransposed& v){
    logger::log_assert(N == v.N, logger::ERROR,
            "vectors have different sizes: {}, {}",
            N, v.N);

    for(size_t i = 0; i < N; ++i){
        V[i] += v[i];
    }
    return *this;
}
VectorTranspose& VectorTranspose::operator-=(const VectorTransposed& v){
    logger::log_assert(N == v.N, logger::ERROR,
            "vectors have different sizes: {}, {}",
            N, v.N);
    for(size_t i = 0; i < N; ++i){
        V[i] -= v[i];
    }
    return *this;
}
VectorTranspose& VectorTranspose::operator*=(Scalar s){
    for(size_t i = 0; i < N; ++i){
        V[i] *= s;
    }
    return *this;
}
VectorTranspose& VectorTranspose::operator/=(Scalar s){
    for(size_t i = 0; i < N; ++i){
        V[i] /= s;
    }
    return *this;
}
VectorTranspose VectorTranspose::operator*(const Matrix& m) const{
    const size_t H = m.get_H();
    const size_t W = m.get_W();
    logger::log_assert(N == m.get_H(), logger::ERROR,
            "vector and matrix have incompatible dimensions: (1, {}), ({}, {})",
            N, H, W);
    VectorTranspose v(W);
    for(size_t i = 0; i < H; ++i){
        const double tmp = V[i];
        for(size_t j = 0; j < W; ++j){
            v[j] += tmp*m(i,j);
        }
    }

    return v;
}
VectorTranspose VectorTranspose::operator*(const MatrixTransposeView& m) const{
    const size_t H = m.get_H();
    const size_t W = m.get_W();
    logger::log_assert(N == m.get_H(), logger::ERROR,
            "vector and matrix have incompatible dimensions: (1, {}), ({}, {})",
            N, H, W);
    VectorTranspose v(W);
    for(size_t i = 0; i < H; ++i){
        const double tmp = V[i];
        for(size_t j = 0; j < W; ++j){
            v[j] += tmp*m(i,j);
        }
    }

    return v;
}
double VectorTranspose::operator*(const VectorNonTransposed& v) const{
    return this->dot(v);
}

VectorView VectorTranspose::T() const{
    return VectorView(N, V);
}

///////////////////////////////////////////////////////////////////////////////
////////////////////////////// VECTOR VIEW ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


VectorView::VectorView(const VectorConst<VectorView>& v):
    VectorConst(v){}
VectorView::VectorView(VectorConst<VectorView>&& v):
    VectorConst(std::forward<VectorConst<VectorView>>(v)){}

VectorView::VectorView(const size_t N, const Scalar* const V):
    VectorConst(N, V){}

Vector VectorView::operator*(Scalar s) const{
    Vector v(*this);
    v *= s;

    return v;
}
Vector VectorView::operator/(Scalar s) const{
    Vector v(*this);
    v /= s;

    return v;
}

Vector VectorView::operator+(const VectorNonTransposed& v) const{
    Vector v2(*this);
    v2 += v;

    return v2;
}
Vector VectorView::operator-(const VectorNonTransposed& v) const{
    Vector v2(*this);
    v2 -= v;

    return v2;
}

Matrix VectorView::operator*(const VectorTransposed& v) const{
    Matrix m(N, v.N);
    for(size_t i = 0; i < N; ++i){
        double tmp = V[i];
        for(size_t j = 0; j < v.N; ++j){
            m(i, j) = tmp*v[j];
        }
    }

    return m;
}

VectorTransposeView VectorView::T() const{
    return VectorTransposeView(N, V);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////// VECTOR TRANSPOSE VIEW ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////


VectorTransposeView::VectorTransposeView(const VectorConst<VectorTransposeView>& v):
    VectorConst(v){}
VectorTransposeView::VectorTransposeView(VectorConst<VectorTransposeView>&& v):
    VectorConst(std::forward<VectorConst<VectorTransposeView>>(v)){}

VectorTransposeView::VectorTransposeView(const size_t N, const Scalar* const V):
    VectorConst(N, V){}

VectorTranspose VectorTransposeView::operator*(Scalar s) const{
    VectorTranspose v(*this);
    v *= s;

    return v;
}
VectorTranspose VectorTransposeView::operator/(Scalar s) const{
    VectorTranspose v(*this);
    v /= s;

    return v;
}

VectorTranspose VectorTransposeView::operator+(const VectorTransposed& v) const{
    VectorTranspose v2(*this);
    v2 += v;

    return v2;
}
VectorTranspose VectorTransposeView::operator-(const VectorTransposed& v) const{
    VectorTranspose v2(*this);
    v2 -= v;

    return v2;
}

VectorTranspose VectorTransposeView::operator*(const Matrix& m) const{
    const size_t H = m.get_H();
    const size_t W = m.get_W();
    logger::log_assert(N == m.get_H(), logger::ERROR,
            "vector and matrix have incompatible dimensions: (1, {}), ({}, {})",
            N, H, W);
    VectorTranspose v(W);
    for(size_t i = 0; i < H; ++i){
        const double tmp = V[i];
        for(size_t j = 0; j < W; ++j){
            v[j] += tmp*m(i,j);
        }
    }

    return v;
}

VectorTranspose VectorTransposeView::operator*(const MatrixTransposeView& m) const{
    const size_t H = m.get_H();
    const size_t W = m.get_W();
    logger::log_assert(N == m.get_H(), logger::ERROR,
            "vector and matrix have incompatible dimensions: (1, {}), ({}, {})",
            N, H, W);
    VectorTranspose v(W);
    for(size_t i = 0; i < H; ++i){
        const double tmp = V[i];
        for(size_t j = 0; j < W; ++j){
            v[j] += tmp*m(i,j);
        }
    }

    return v;
}

double VectorTransposeView::operator*(const VectorNonTransposed& v) const{
    return this->dot(v);
}

VectorView VectorTransposeView::T() const{
    return VectorView(N, V);
}

std::ostream& operator<<(std::ostream& output, const VectorAgnostic& v){
    const size_t N = v.get_N();
    for(size_t i = 0; i < N; ++i){
        output << v[i] << " ";
    }
    output << std::endl;

    return output;
}

///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// EIGEN /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

Eigen::Eigen(Matrix A):
    eigenvectors(A.get_H(), A.get_W()), eigenvalues(A.get_W()){

    logger::log_assert(A.get_H() == A.get_W(), logger::ERROR,
            "unable to calculate eigenvectors and eigenvalues: matrix is not square");

    const size_t N = A.get_W();

    std::vector<int> ISUPPZ(2*N);
    // std::vector<double> WORK(26*N);
    // std::vector<int> IWORK(10*N);
    int M;

    int info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'A', 'L', N, A.data(), N, 0, 0, 0, 0, 1e-7, &M,
         eigenvalues.data(), eigenvectors.data(), N, ISUPPZ.data());

    logger::log_assert(info == 0, logger::ERROR,
            "LAPACK returned {} while calculating eigenvalues and eigenvectors", info);
}

Matrix Eigen::square_root() const{
    const size_t N = this->eigenvalues.get_N();
    Matrix D(Matrix::diag(this->eigenvalues));

    for(size_t i = 0; i < N; ++i){
        D(i,i) = std::sqrt(D(i,i));
    }

    const Matrix& Z = this->eigenvectors;

    return Z*D*Z.T();
}

}
