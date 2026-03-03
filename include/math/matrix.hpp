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
#include "math/vector_view.hpp"
#include <cstddef>
#include <initializer_list>
#include <numeric>
#include <ostream>
#include <utility>
#include <vector>

namespace math{

typedef double Scalar;
class MatrixTransposeView;
class Vector;
class VectorView;
class VectorTranspose;
class VectorTransposeView;
class VectorAgnostic;
class VectorTransposed;
class VectorNonTransposed;

class MatrixSlice{
    public:
    friend class Matrix;
    MatrixSlice() = delete;
    ~MatrixSlice() = default;

    inline size_t get_H() const{
        return this->pos_i.size();
    }
    inline size_t get_W() const{
        return this->pos_j.size();
    }

    inline Scalar at(const size_t i, const size_t j) const{
        return this->M[pos_i[i]*W + pos_j[j]];
    }
    inline Scalar& at(const size_t i, const size_t j){
        return this->M[pos_i[i]*W + pos_j[j]];
    }
    inline Scalar operator()(const size_t i, const size_t j) const{
        return this->M[pos_i[i]*W + pos_j[j]];
    }
    inline Scalar& operator()(const size_t i, const size_t j){
        return this->M[pos_i[i]*W + pos_j[j]];
    }

    MatrixSlice& operator+=(const Matrix& m);
    MatrixSlice& operator-=(const Matrix& m);
    MatrixSlice& operator+=(const MatrixTransposeView& m);
    MatrixSlice& operator-=(const MatrixTransposeView& m);
    MatrixSlice& operator+=(const MatrixSlice& m);
    MatrixSlice& operator-=(const MatrixSlice& m);

    private:
    MatrixSlice(Scalar* m, const size_t W, std::vector<size_t> pos_i, std::vector<size_t> pos_j):
        pos_i(std::move(pos_i)),
        pos_j(std::move(pos_j)),
        W(W),
        M(m)
    {}
    const std::vector<size_t> pos_i, pos_j;
    const size_t W;
    Scalar* M;
};

class Matrix{
    public:
    Matrix() = default;
    Matrix(Scalar* m, size_t H, size_t W);
    Matrix(std::vector<Scalar> m, size_t H, size_t W);
    Matrix(std::initializer_list<Scalar> m, size_t H, size_t W);
    Matrix(size_t H, size_t W, Scalar s = 0);
    Matrix(const MatrixSlice& m);
    ~Matrix();
    Matrix(const Matrix& m);
    Matrix(Matrix&& m);
    Matrix(const MatrixTransposeView& m);

    static Matrix identity(size_t N);
    static Matrix diag(const Vector& v);

    inline std::vector<size_t> all_i() const{
        std::vector<size_t> pos_i(this->H);
        std::iota(pos_i.begin(), pos_i.end(), 0);
        return pos_i;
    }
    inline std::vector<size_t> all_j() const{
        std::vector<size_t> pos_j(this->W);
        std::iota(pos_j.begin(), pos_j.end(), 0);
        return pos_j;
    }

    inline MatrixSlice slice(std::vector<size_t> pos_i, std::vector<size_t> pos_j){
        return MatrixSlice(this->M, this->W, pos_i, pos_j);
    }
    inline MatrixSlice slice(std::vector<size_t> pos_i, std::vector<size_t> pos_j) const{
        return MatrixSlice(this->M, this->W, pos_i, pos_j);
    }

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
    Matrix get_inverted_LU() const;
    void invert_LU();
    Matrix get_inverted_cholesky() const;
    void invert_cholesky();
    bool is_equal(const Matrix& m, Scalar eps = 1e-7) const;
    bool is_equal(const MatrixTransposeView& m, Scalar eps = 1e-7) const;
    Scalar determinant() const;

    MatrixTransposeView T() const;
    
    Matrix& operator=(const Matrix& m);
    Matrix& operator=(Matrix&& m);
    Matrix& operator=(const MatrixTransposeView& m);
    Matrix& operator=(const MatrixSlice& m);

    inline const Matrix& operator+() const{
        return *this;
    }
    Matrix operator-() const;

    Matrix& operator+=(const Matrix& m);
    Matrix& operator-=(const Matrix& m);
    Matrix& operator+=(const MatrixTransposeView& m);
    Matrix& operator-=(const MatrixTransposeView& m);
    Matrix& operator+=(const MatrixSlice& m);
    Matrix& operator-=(const MatrixSlice& m);

    Matrix operator+(const Matrix& m) const;
    Matrix operator-(const Matrix& m) const;
    Matrix operator*(const Matrix& m) const;
    Matrix operator+(const MatrixSlice& m) const;
    Matrix operator-(const MatrixSlice& m) const;
    Matrix operator*(const MatrixSlice& m) const;
    Matrix operator*(const MatrixTransposeView& m) const;
    Vector operator*(const VectorNonTransposed& v) const;
    Vector operator*(const VectorSlice& v) const;
    Vector operator*(const VectorSliceView& v) const;
    template<typename P>
    Vector operator*(const VectorSliceGeneralView<Scalar, P>& v) const;

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
    Matrix operator+(const MatrixSlice& m) const;
    Matrix operator-(const MatrixSlice& m) const;
    Matrix operator*(const MatrixSlice& m) const;
    Vector operator*(const VectorNonTransposed& v) const;
    Vector operator*(const VectorSlice& v) const;
    Vector operator*(const VectorSliceView& v) const;
    template<typename P>
    Vector operator*(const VectorSliceGeneralView<Scalar, P>& v) const;

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

class LU{
    public:
    // MODIFIES VALUE OF MATRIX IF copy == false !!!!!!!
    LU(Matrix& m, bool copy = false);
    ~LU();

    Scalar determinant() const;
    void solve(Vector& v) const;
    void solve(Matrix& m) const;
    void solve_transposed(Vector& v) const;
    void solve_transposed(Matrix& m) const;

    private:
    size_t N;
    Scalar* M;
    std::vector<int> ipiv;
    bool copy;
};

class Cholesky{
    public:
    // MODIFIES VALUE OF MATRIX IF copy == false !!!!!!!
    Cholesky(Matrix& m, bool copy = false);
    ~Cholesky();

    Scalar determinant() const;
    void solve(Vector& v) const;
    void solve(Matrix& m) const;

    private:
    size_t N;
    Scalar* M;
    std::vector<int> ipiv;
    bool copy;
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

class VectorAgnostic{
    public:
    VectorAgnostic() = delete;
    VectorAgnostic(const Vector& v);
    VectorAgnostic(const VectorView& v);
    VectorAgnostic(const VectorTranspose& v);
    VectorAgnostic(const VectorTransposeView& v);
    VectorAgnostic(const VectorTransposed& v);
    VectorAgnostic(const VectorNonTransposed& v);

    inline size_t get_N() const{
        return this->N;
    }

    inline Scalar at(const size_t i) const{
        return this->V[i];
    }
    inline Scalar operator()(const size_t i) const{
        return this->V[i];
    }
    inline Scalar operator[](const size_t i) const{
        return this->V[i];
    }

    const size_t N;
    const Scalar* V;
};

class VectorNonTransposed{
    public:
    VectorNonTransposed() = delete;
    VectorNonTransposed(const Vector& v);
    VectorNonTransposed(const VectorView& v);

    inline size_t get_N() const{
        return this->N;
    }

    inline Scalar at(const size_t i) const{
        return this->V[i];
    }
    inline Scalar operator()(const size_t i) const{
        return this->V[i];
    }
    inline Scalar operator[](const size_t i) const{
        return this->V[i];
    }

    const size_t N;
    const Scalar* V;
};

class VectorTransposed{
    public:
    VectorTransposed() = delete;
    VectorTransposed(const VectorTranspose& v);
    VectorTransposed(const VectorTransposeView& v);

    inline size_t get_N() const{
        return this->N;
    }

    inline Scalar at(const size_t i) const{
        return this->V[i];
    }
    inline Scalar operator()(const size_t i) const{
        return this->V[i];
    }
    inline Scalar operator[](const size_t i) const{
        return this->V[i];
    }

    const size_t N;
    const Scalar* V;
};

template<typename T1, typename T2>
class VectorBase{
    public:
    virtual ~VectorBase() = default;

    inline T1 get_N() const{
        return this->N;
    }

    double norm() const{
        double sum = 0;
        for(size_t i = 0; i < N; ++i){
            sum += V[i] * V[i];
        }
        return std::sqrt(sum);
    }
    bool is_equal(const VectorAgnostic& v, Scalar eps = 1e-7) const{
        if(N != v.N) return false;
        for(size_t i = 0; i < N; ++i){
            if(std::abs(V[i] - v.V[i]) >= eps){
                return false;
            }
        }
        return true;
    }
    double dot(const VectorAgnostic& v) const{
        logger::log_assert(N == v.N, logger::ERROR,
                "vectors have different sizes: {}, {}",
                N, v.N);
        double sum = 0;
        for(size_t i = 0; i < N; ++i){
            sum += V[i] * v.V[i];
        }

        return sum;
    }
    double dot(const VectorSlice& v) const{
        logger::log_assert(N == v.size(), logger::ERROR,
                "vectors have different sizes: {}, {}",
                N, v.size());
        double sum = 0;
        for(size_t i = 0; i < N; ++i){
            sum += V[i] * v[i];
        }

        return sum;
    }
    double dot(const VectorSliceView& v) const{
        logger::log_assert(N == v.size(), logger::ERROR,
                "vectors have different sizes: {}, {}",
                N, v.size());
        double sum = 0;
        for(size_t i = 0; i < N; ++i){
            sum += V[i] * v[i];
        }

        return sum;
    }
    template<typename P>
    double dot(const VectorSliceGeneralView<Scalar, P>& v) const{
        logger::log_assert(N == v.size(), logger::ERROR,
                "vectors have different sizes: {}, {}",
                N, v.size());
        double sum = 0;
        for(size_t i = 0; i < N; ++i){
            sum += V[i] * v[i];
        }

        return sum;
    }

    bool operator==(const VectorAgnostic& v) const{
        if(N != v.N) return false;
        for(size_t i = 0; i < N; ++i){
            if(V[i] != v[i]) return false;
        }
        return true;
    }
    bool operator!=(const VectorAgnostic& v) const{
        return !(*this == v);
    }

    protected:
    VectorBase() = default;
    VectorBase(T1 N, T2* V):
        N(N), V(V){}
    T1 N = 0;
    T2* V = nullptr;
};

template<class T>
class VectorConst : public VectorBase<const size_t, const Scalar>{
    public:
    virtual ~VectorConst() = default;

    inline Scalar at(const size_t i) const{
        return this->V[i];
    }
    inline Scalar operator()(const size_t i) const{
        return this->V[i];
    }
    inline Scalar operator[](const size_t i) const{
        return this->V[i];
    }
    inline const Scalar* data() const{
        return this->V;
    }
    
    inline const T& operator+() const{
        return *this;
    }
    T operator-() const{
        T v(*this);

        return -1.0*v;
    }

    protected:
    VectorConst() = default;
    VectorConst(size_t N, const Scalar* V):
        VectorBase(N, V){}

};

template<class T>
class VectorMut : public VectorBase<size_t, Scalar>{
    public:
    virtual ~VectorMut(){
        delete[] this->V;
    }

    inline Scalar& at(const size_t i){
        return this->V[i];
    }
    inline Scalar& operator()(const size_t i){
        return this->V[i];
    }
    inline Scalar& operator[](const size_t i){
        return this->V[i];
    }
    inline Scalar* data(){
        return this->V;
    }
    inline Scalar at(const size_t i) const{
        return this->V[i];
    }
    inline Scalar operator()(const size_t i) const{
        return this->V[i];
    }
    inline Scalar operator[](const size_t i) const{
        return this->V[i];
    }
    inline const Scalar* data() const{
        return this->V;
    }

    void fill(Scalar s){
        std::fill(V, V + N, s);
    }

    void normalize(){
        const double norm = this->norm();
        for(size_t i = 0; i < this->N; ++i){
            this->V[i] /= norm;
        }
    }

    T normalized() const{
        T nv(*this);
        nv.normalize();

        return nv;
    }

    T& operator=(const VectorAgnostic& v){
        this->N = v.N;
        delete[] this->V;
        this->V = new Scalar[N];

        std::copy(v.V, v.V + N, this->V);

        return *this;
    }
    T& operator=(const T& v){
        this->N = v.N;
        delete[] this->V;
        this->V = new Scalar[N];

        std::copy(v.V, v.V + N, this->V);

        return *this;
    }
    T& operator=(T&& v){
        this->N = v.N;
        delete[] this->V;
        this->V = v.V;
        v.V = nullptr;

        return *this;
    }
    T& operator=(const VectorMut<T>& v){
        this->N = v.N;
        delete[] this->V;
        this->V = new Scalar[N];

        std::copy(v.V, v.V + N, this->V);

        return *this;
    }
    T& operator=(VectorMut<T>&& v){
        this->N = v.N;
        delete[] this->V;
        this->V = v.V;
        v.V = nullptr;

        return *this;
    }
    
    inline const T& operator+() const{
        return *this;
    }
    T operator-() const{
        T v(*this);

        return -1.0*v;
    }

    protected:
    VectorMut() = default;
    VectorMut(size_t N):
        VectorBase(N, new Scalar[N]){}
    VectorMut(size_t N, Scalar*&& v):
        VectorBase(N, v){
        v = nullptr;     
    }
    VectorMut(const VectorMut& v):
        VectorBase(v.N, new Scalar[v.N]){
        std::copy(v.V, v.V + N, this->V);
    }
    VectorMut(VectorMut&& v):
        VectorBase(v.N, v.V){

        v.V = nullptr;
    }
};

class Vector : public VectorMut<Vector>{
    public:
    Vector() = default;
    Vector(Scalar* v, size_t N);
    Vector(std::vector<Scalar> v);
    Vector(std::initializer_list<Scalar> v);
    Vector(size_t N, Scalar s = 0);
    virtual ~Vector() = default;
    Vector(const VectorAgnostic& v);
    Vector(const Vector& v);
    Vector(Vector&& v);
    Vector(const VectorMut<Vector>& v);
    Vector(VectorMut<Vector>&& v);

    VectorTransposeView T() const;

    Vector cross(const VectorAgnostic& v) const;

    Vector operator*(Scalar s) const;
    Vector operator/(Scalar s) const;

    Vector operator+(const VectorNonTransposed& v) const;
    Vector operator-(const VectorNonTransposed& v) const;

    Vector& operator*=(Scalar s);
    Vector& operator/=(Scalar s);

    Vector& operator+=(const VectorNonTransposed& v);
    Vector& operator-=(const VectorNonTransposed& v);

    Vector& operator=(const Vector& v);
    Vector& operator=(Vector&& v);

    Matrix operator*(const VectorTransposed& v) const;
};

class VectorTranspose : public VectorMut<VectorTranspose>{
    friend class VectorTransposeView;
    public:
    ~VectorTranspose() = default;
    VectorTranspose(const VectorMut<VectorTranspose>& v);
    VectorTranspose(VectorMut<VectorTranspose>&& v);
    VectorTranspose(const VectorAgnostic& v);

    VectorView T() const;

    VectorTranspose operator+(const VectorTransposed& v) const;
    VectorTranspose operator-(const VectorTransposed& v) const;
    VectorTranspose operator*(const Matrix& m) const;
    VectorTranspose operator*(const MatrixTransposeView& m) const;

    VectorTranspose operator*(Scalar s) const;
    VectorTranspose operator/(Scalar s) const;

    VectorTranspose& operator*=(Scalar s);
    VectorTranspose& operator/=(Scalar s);

    VectorTranspose& operator+=(const VectorTransposed& v);
    VectorTranspose& operator-=(const VectorTransposed& v);

    double operator*(const VectorNonTransposed& v) const;

    private:
    VectorTranspose& operator=(const VectorTranspose& v);
    VectorTranspose& operator=(VectorTranspose&& v);
    VectorTranspose() = default;
    VectorTranspose(const VectorTranspose&);
    VectorTranspose(VectorTranspose&&);
    VectorTranspose(size_t N, Scalar s = 0);
};

class VectorView : public VectorConst<VectorView>{
    friend class Vector;
    friend class VectorTranspose;
    friend class VectorTransposeView;
    public:
    ~VectorView() = default;
    VectorView(const VectorView&) = delete;
    VectorView(VectorView&&) = delete;
    VectorView& operator=(const VectorView& v) = delete;
    VectorView& operator=(VectorView&& v) = delete;

    VectorTransposeView T() const;

    Vector operator+(const VectorNonTransposed& v) const;
    Vector operator-(const VectorNonTransposed& v) const;

    Vector operator*(Scalar s) const;
    Vector operator/(Scalar s) const;

    Matrix operator*(const VectorTransposed& v) const;

    private:
    VectorView() = default;
    VectorView(const size_t N, const Scalar* const V);
    VectorView(const VectorConst<VectorView>& v);
    VectorView(VectorConst<VectorView>&& v);
};

class VectorTransposeView : public VectorConst<VectorTransposeView>{
    friend class Vector;
    friend class VectorView;
    public:
    ~VectorTransposeView() = default;
    VectorTransposeView(const VectorTransposeView&) = delete;
    VectorTransposeView(VectorTransposeView&&) = delete;
    VectorTransposeView& operator=(const VectorTransposeView& v) = delete;
    VectorTransposeView& operator=(VectorTransposeView&& v) = delete;

    VectorTranspose operator+(const VectorTransposed& v) const;
    VectorTranspose operator-(const VectorTransposed& v) const;
    VectorTranspose operator*(const Matrix& m) const;
    VectorTranspose operator*(const MatrixTransposeView& m) const;

    VectorView T() const;

    VectorTranspose operator*(Scalar s) const;
    VectorTranspose operator/(Scalar s) const;

    double operator*(const VectorNonTransposed& v) const;

    private:
    VectorTransposeView(const VectorConst<VectorTransposeView>& v);
    VectorTransposeView(VectorConst<VectorTransposeView>&& v);
    VectorTransposeView() = default;
    VectorTransposeView(const size_t N, const Scalar* const V);
};

std::ostream& operator<<(std::ostream& output, const VectorAgnostic& v);

inline Vector operator*(Scalar s, const Vector& v){
    return v*s;
}
inline Vector operator*(Scalar s, const VectorView& v){
    return v*s;
}
inline VectorTranspose operator*(Scalar s, const VectorTranspose& v){
    return v*s;
}
inline VectorTranspose operator*(Scalar s, const VectorTransposeView& v){
    return v*s;
}

inline Vector&& operator*(Scalar s, Vector&& v){
    v *= s;
    return std::forward<Vector>(v);
}
inline Vector operator*(Scalar s, VectorView&& v){
    return v*s;
}
inline VectorTranspose&& operator*(Scalar s, VectorTranspose&& v){
    v *= s;
    return std::forward<VectorTranspose>(v);
}
inline VectorTranspose operator*(Scalar s, VectorTransposeView&& v){
    return v*s;
}


class Eigen{
    public:
    Eigen(Matrix A);

    Matrix square_root() const;

    Matrix eigenvectors;
    Vector eigenvalues;
};


template<typename P>
Vector Matrix::operator*(const VectorSliceGeneralView<Scalar, P>& v) const{
    logger::log_assert(W == v.size(), logger::ERROR,
            "matrix and vector have incompatible dimensions: ({}, {}), ({}, 1)",
            H, W, v.size());
    Vector v2(H);
    for(size_t j = 0; j < W; ++j){
        const double tmp = v[j];
        for(size_t i = 0; i < H; ++i){
            v2[i] += this->at(i,j)*tmp;
        }
    }

    return v2;
}
template<typename P>
Vector MatrixTransposeView::operator*(const VectorSliceGeneralView<Scalar, P>& v) const{
    logger::log_assert(W == v.size(), logger::ERROR,
            "matrix and vector have incompatible dimensions: ({}, {}), ({}, 1)",
            H, W, v.size());
    Vector v2(H);
    for(size_t j = 0; j < W; ++j){
        const double tmp = v[j];
        for(size_t i = 0; i < H; ++i){
            v2[i] += this->at(i,j)*tmp;
        }
    }

    return v2;
}

}

#endif
