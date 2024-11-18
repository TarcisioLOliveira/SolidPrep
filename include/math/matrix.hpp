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

#include <cstddef>
#include <initializer_list>
#include <ostream>
#include <vector>

namespace math{

typedef double Scalar;
class MatrixTransposeView;

class Matrix{
    public:
    Matrix() = default;
    Matrix(Scalar* m, size_t W, size_t H);
    Matrix(std::vector<Scalar> m, size_t W, size_t H);
    Matrix(std::initializer_list<Scalar> m, size_t W, size_t H);
    Matrix(size_t W, size_t H, Scalar s = 0);
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

    void fill(Scalar s);
    Matrix get_inverted() const;
    void invert();
    bool is_equal(const Matrix& m, Scalar eps = 1e-7) const;
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

    Matrix& operator*=(Scalar s);
    Matrix& operator/=(Scalar s);
    Matrix operator*(Scalar s) const;
    Matrix operator/(Scalar s) const;

    bool operator==(const Matrix& m) const;
    bool operator!=(const Matrix& m) const;

    friend std::ostream& operator<<(std::ostream& output, const Matrix& m);

    private:
    size_t W = 0, H = 0;
    Scalar* M = nullptr;
};

Matrix operator*(Scalar s, const Matrix& m);

class MatrixTransposeView{
    friend class Matrix;
    public:
    ~MatrixTransposeView() = default;
    MatrixTransposeView(const MatrixTransposeView&) = delete;
    MatrixTransposeView(MatrixTransposeView&&) = delete;
    MatrixTransposeView& operator=(const MatrixTransposeView& m);
    MatrixTransposeView& operator=(MatrixTransposeView&& m);

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

    Matrix operator+(const Matrix& m) const;
    Matrix operator-(const Matrix& m) const;
    Matrix operator*(const Matrix& m) const;

    Matrix operator*(Scalar s) const;
    Matrix operator/(Scalar s) const;

    private:
    MatrixTransposeView() = default;
    MatrixTransposeView(const size_t non_T_W, const size_t non_T_H, const Scalar* const M);

    const size_t W = 0, H = 0;
    const Scalar* M = nullptr;
};


}

#endif
