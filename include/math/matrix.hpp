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

#ifndef MATH_MATRIX_HPP
#define MATH_MATRIX_HPP

#include "logger.hpp"
#include <cstddef>
#include <initializer_list>
#include <ostream>
#include <utility>
#include <vector>

namespace math{

typedef double Scalar;
class MatrixTransposeView;

class Matrix{
    public:
    Matrix() = default;
    Matrix(Scalar* m, size_t H, size_t W);
    Matrix(std::vector<Scalar> m, size_t H, size_t W);
    Matrix(std::initializer_list<Scalar> m, size_t H, size_t W);
    Matrix(size_t H, size_t W, Scalar s = 0);
    ~Matrix();
    Matrix(const Matrix& m);
    Matrix(Matrix&& m);
    Matrix(const MatrixTransposeView& m);

    static Matrix identity(size_t N);

    inline size_t get_H() const{
        return this->H;
    }
    inline size_t get_W() const{
        return this->W;
    }

    inline Scalar at(const size_t i, const size_t j) const{
        return this->M[i*W + j];
    }
    inline Scalar& at(const size_t i, const size_t j){
        return this->M[i*W + j];
    }
    inline Scalar operator()(const size_t i, const size_t j) const{
        return this->M[i*W + j];
    }
    inline Scalar& operator()(const size_t i, const size_t j){
        return this->M[i*W + j];
    }
    inline Scalar* data(){
        return this->M;
    }
    inline const Scalar* data() const{
        return this->M;
    }

    void fill(Scalar s);
    Matrix get_inverted() const;
    void invert();
    bool is_equal(const Matrix& m, Scalar eps = 1e-7) const;
    bool is_equal(const MatrixTransposeView& m, Scalar eps = 1e-7) const;
    Scalar determinant() const;

    MatrixTransposeView T() const;
    
    Matrix& operator=(const Matrix& m);
    Matrix& operator=(Matrix&& m);
    Matrix& operator=(const MatrixTransposeView& m);

    inline const Matrix& operator+() const{
        return *this;
    }
    Matrix operator-() const;

    Matrix& operator+=(const Matrix& m);
    Matrix& operator-=(const Matrix& m);
    Matrix& operator+=(const MatrixTransposeView& m);
    Matrix& operator-=(const MatrixTransposeView& m);

    Matrix operator+(const Matrix& m) const;
    Matrix operator-(const Matrix& m) const;
    Matrix operator*(const Matrix& m) const;
    Matrix operator*(const MatrixTransposeView& m) const;

    Matrix& operator*=(Scalar s);
    Matrix& operator/=(Scalar s);
    Matrix operator*(Scalar s) const;
    Matrix operator/(Scalar s) const;

    bool operator==(const Matrix& m) const;
    bool operator!=(const Matrix& m) const;
    bool operator==(const MatrixTransposeView&& m) const;
    bool operator!=(const MatrixTransposeView&& m) const;


    private:
    size_t H = 0, W = 0;
    Scalar* M = nullptr;
};

std::ostream& operator<<(std::ostream& output, const Matrix& m);
std::ostream& operator<<(std::ostream& output, const MatrixTransposeView& m);

class MatrixTransposeView{
    friend class Matrix;
    public:
    ~MatrixTransposeView() = default;
    MatrixTransposeView(const MatrixTransposeView&) = delete;
    MatrixTransposeView(MatrixTransposeView&&) = delete;
    MatrixTransposeView& operator=(const MatrixTransposeView& m) = delete;
    MatrixTransposeView& operator=(MatrixTransposeView&& m) = delete;

    inline size_t get_H() const{
        return this->H;
    }
    inline size_t get_W() const{
        return this->W;
    }

    inline Scalar at(const size_t i, const size_t j) const{
        return this->M[j*H + i];
    }
    inline Scalar operator()(const size_t i, const size_t j) const{
        return this->M[j*H + i];
    }
    inline const Scalar* data() const{
        return this->M;
    }
    inline bool is_equal(const Matrix& m, Scalar eps = 1e-7) const{
        return m.is_equal(*this, eps);
    }
    bool is_equal(const MatrixTransposeView& m, Scalar eps = 1e-7) const;

    Matrix operator+(const Matrix& m) const;
    Matrix operator-(const Matrix& m) const;
    Matrix operator*(const Matrix& m) const;

    Matrix operator*(Scalar s) const;
    Matrix operator/(Scalar s) const;

    inline bool operator==(const Matrix& m) const{
        return m == *this;
    }
    inline bool operator!=(const Matrix& m) const{
        return m != *this;
    }
    bool operator==(const MatrixTransposeView&& m) const;
    bool operator!=(const MatrixTransposeView&& m) const;

    private:
    MatrixTransposeView() = default;
    MatrixTransposeView(const size_t non_T_H, const size_t non_T_W, const Scalar* const M);

    const size_t H = 0, W = 0;
    const Scalar* M = nullptr;
};

inline Matrix operator*(Scalar s, const Matrix& m){
    return m*s;
}
inline Matrix operator*(Scalar s, const MatrixTransposeView& m){
    return m*s;
}

inline Matrix&& operator*(Scalar s, Matrix&& m){
    m *= s;
    return std::forward<Matrix>(m);
}
inline Matrix operator*(Scalar s, MatrixTransposeView&& m){
    return m*s;
}

}

#endif
