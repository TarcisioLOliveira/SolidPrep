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


#include "math/vector_view.hpp"
#include "logger.hpp"
#include "math/matrix.hpp"
#include "element.hpp"

namespace math{

VectorSlice::VectorSlice(Scalar* data, std::vector<size_t> points):
    points(std::move(points)), data(data){}
VectorSlice::VectorSlice(std::vector<Scalar>& data, std::vector<size_t> points):
    points(std::move(points)), data(data.data()){}
VectorSlice::VectorSlice(math::Vector& data, std::vector<size_t> points):
    points(std::move(points)), data(data.data()){}

math::Vector VectorSlice::as_math_vector() const{
    math::Vector v(this->size());
    for(size_t i = 0; i < this->size(); ++i){
        v[i] = this->data[this->points[i]];
    }

    return v;
}
std::vector<Scalar> VectorSlice::as_std_vector() const{
    std::vector<Scalar> v(this->size());
    for(size_t i = 0; i < this->size(); ++i){
        v[i] = this->data[this->points[i]];
    }

    return v;
}

void VectorSlice::assign(const VectorSlice& v){
    const size_t N = std::min(this->size(), v.size());
    for(size_t i = 0; i < N; ++i){
        this->data[this->points[i]] = v[i];
    }
}
void VectorSlice::assign(const math::Vector& v){
    const size_t N = std::min(this->size(), v.get_N());
    for(size_t i = 0; i < N; ++i){
        this->data[this->points[i]] = v[i];
    }
}
void VectorSlice::assign(const std::vector<Scalar>& v){
    const size_t N = std::min(this->size(), v.size());
    for(size_t i = 0; i < N; ++i){
        this->data[this->points[i]] = v[i];
    }
}

VectorSlice& VectorSlice::operator+=(const VectorSlice& v){
    const size_t N = std::min(this->size(), v.size());
    for(size_t i = 0; i < N; ++i){
        this->data[this->points[i]] += v[i];
    }
    return *this;
}
VectorSlice& VectorSlice::operator+=(const math::Vector& v){
    const size_t N = std::min(this->size(), v.get_N());
    for(size_t i = 0; i < N; ++i){
        this->data[this->points[i]] += v[i];
    }
    return *this;
}
VectorSlice& VectorSlice::operator+=(const std::vector<Scalar>& v){
    const size_t N = std::min(this->size(), v.size());
    for(size_t i = 0; i < N; ++i){
        this->data[this->points[i]] += v[i];
    }
    return *this;
}
VectorSlice& VectorSlice::operator-=(const VectorSlice& v){
    const size_t N = std::min(this->size(), v.size());
    for(size_t i = 0; i < N; ++i){
        this->data[this->points[i]] -= v[i];
    }
    return *this;
}
VectorSlice& VectorSlice::operator-=(const math::Vector& v){
    const size_t N = std::min(this->size(), v.get_N());
    for(size_t i = 0; i < N; ++i){
        this->data[this->points[i]] -= v[i];
    }
    return *this;
}
VectorSlice& VectorSlice::operator-=(const std::vector<Scalar>& v){
    const size_t N = std::min(this->size(), v.size());
    for(size_t i = 0; i < N; ++i){
        this->data[this->points[i]] -= v[i];
    }
    return *this;
}


VectorSliceView::VectorSliceView(const Scalar* data, std::vector<size_t> points):
    points(std::move(points)), data(data){}
VectorSliceView::VectorSliceView(const std::vector<Scalar>& data, std::vector<size_t> points):
    points(std::move(points)), data(data.data()){}
VectorSliceView::VectorSliceView(const math::Vector& data, std::vector<size_t> points):
    points(std::move(points)), data(data.data()){}

math::Vector VectorSliceView::as_math_vector() const{
    math::Vector v(this->size());
    for(size_t i = 0; i < this->size(); ++i){
        v[i] = this->data[this->points[i]];
    }

    return v;
}
std::vector<Scalar> VectorSliceView::as_std_vector() const{
    std::vector<Scalar> v(this->size());
    for(size_t i = 0; i < this->size(); ++i){
        v[i] = this->data[this->points[i]];
    }

    return v;
}

}
