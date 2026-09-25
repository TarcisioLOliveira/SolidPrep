/*
 *   Copyright (C) 2026 Tarcísio Ladeia de Oliveira.
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

#ifndef VEC_WRAP_HPP
#define VEC_WRAP_HPP

#include <vector>

namespace utils{

template<typename T>
class VectorWrapper{
    public:
    std::vector<T>& vec;
    VectorWrapper() = default;
    VectorWrapper(std::vector<T>& v):vec(v){}

    inline VectorWrapper& operator+=(const std::vector<T>& v){
        for(size_t i = 0; i < vec.size(); ++i){
            this->vec[i] += v[i];
        }

        return *this;
    }
    inline VectorWrapper& operator+=(const VectorWrapper& v){
        for(size_t i = 0; i < vec.size(); ++i){
            this->vec[i] += v.vec[i];
        }

        return *this;
    }
    inline VectorWrapper& operator*=(const double s){
        for(size_t i = 0; i < vec.size(); ++i){
            this->vec[i] *= s;
        }

        return *this;
    }
};

template<typename T>
inline std::vector<T> operator*(const double s, const VectorWrapper<T>& v){
    std::vector<T> v2 = v.vec;
    for(size_t i = 0; i < v2.size(); ++i){
        v2[i] *= s;
    }

    return v2;
}

}

#endif
