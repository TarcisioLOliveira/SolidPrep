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
class Vector;
class VectorView;
class VectorTranspose;
class VectorTransposeView;
class VectorAgnostic;
class VectorTransposed;
class VectorNonTransposed;

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
    Vector operator*(const VectorNonTransposed& v) const;

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
    Vector operator*(const VectorNonTransposed& v) const;

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
    virtual ~VectorMut() = default;

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
    VectorMut(size_t N, Scalar* V):
        VectorBase(N, V){}
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
    virtual ~Vector();
    Vector(const VectorAgnostic& v);
    Vector(const Vector& v);
    Vector(Vector&& v);
    Vector(const VectorMut<Vector>& v);
    Vector(VectorMut<Vector>&& v);

    VectorTransposeView T() const;

    Vector operator*(Scalar s) const;
    Vector operator/(Scalar s) const;

    Vector operator+(const VectorNonTransposed& v) const;
    Vector operator-(const VectorNonTransposed& v) const;

    Vector& operator*=(Scalar s);
    Vector& operator/=(Scalar s);

    Vector& operator+=(const VectorNonTransposed& v);
    Vector& operator-=(const VectorNonTransposed& v);

    Matrix operator*(const VectorTransposed& v) const;
};

class VectorTranspose : public VectorMut<VectorTranspose>{
    friend class VectorTransposeView;
    public:
    ~VectorTranspose();
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

}

#endif
