/*
 *   Copyright (C) 2025 Tarcísio Ladeia de Oliveira.
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

#ifndef MATH_VECTOR_VIEW_HPP
#define MATH_VECTOR_VIEW_HPP

#include <vector>

namespace math{

class Vector;

class VectorSlice{
    public:
    typedef double Scalar;
    VectorSlice() = default;
    VectorSlice(Scalar* data, std::vector<size_t> points);
    VectorSlice(std::vector<Scalar>& data, std::vector<size_t> points);
    VectorSlice(math::Vector& data, std::vector<size_t> points);
    ~VectorSlice() = default;

    math::Vector as_math_vector() const;
    std::vector<Scalar> as_std_vector() const;

    void assign(const VectorSlice& v);
    void assign(const math::Vector& v);
    void assign(const std::vector<Scalar>& v);

    inline size_t size() const{
        return this->points.size();
    }

    VectorSlice& operator+=(const VectorSlice& v);
    VectorSlice& operator+=(const math::Vector& v);
    VectorSlice& operator+=(const std::vector<Scalar>& v);
    VectorSlice& operator-=(const VectorSlice& v);
    VectorSlice& operator-=(const math::Vector& v);
    VectorSlice& operator-=(const std::vector<Scalar>& v);

    inline Scalar& operator[](const size_t i){
        return this->data[this->points[i]];
    }
    inline Scalar operator[](const size_t i) const{
        return this->data[this->points[i]];
    }

    private:
    std::vector<size_t> points;
    Scalar* data;
};

class VectorSliceView{
    public:
    typedef double Scalar;
    VectorSliceView() = default;
    VectorSliceView(const Scalar* data, std::vector<size_t> points);
    VectorSliceView(const std::vector<Scalar>& data, std::vector<size_t> points);
    VectorSliceView(const math::Vector& data, std::vector<size_t> points);
    ~VectorSliceView() = default;

    math::Vector as_math_vector() const;
    std::vector<Scalar> as_std_vector() const;

    inline size_t size() const{
        return this->points.size();
    }

    inline Scalar operator[](const size_t i) const{
        return this->data[this->points[i]];
    }

    private:
    std::vector<size_t> points;
    const Scalar* data;
};

template<typename T, typename P>
class VectorSliceGeneralView{
    public:
    VectorSliceGeneralView() = default;
    VectorSliceGeneralView(const T* const data, std::vector<P> points):
        points(std::move(points)), data(data){}
    VectorSliceGeneralView(const std::vector<T>& data, std::vector<P> points):
        points(std::move(points)), data(data.data()){}
    ~VectorSliceGeneralView() = default;

    inline size_t size() const{
        return this->points.size();
    }

    inline T operator[](const size_t i) const{
        return this->get(i);
    }

    private:
    inline T get(const size_t i) const{
        if constexpr (std::is_signed<T>::value){
            const P p(this->points[i]);
            if(p >= 0){
                return this->data[p];
            } else {
                return 0;
            }
        } else {
            const P p(this->points[i]);
            return this->data[p];
        }
    }

    const std::vector<P> points;
    const T* const data;
};

}

#endif
