/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
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

#ifndef SPARSE_MATRIX_HPP
#define SPARSE_MATRIX_HPP

#include "math/matrix.hpp"
#include <string>
#include <unordered_map>
#include <vector>
#include <cstddef>
#include <Eigen/Sparse>

namespace utils{

// from boost (functional/hash):
// see http://www.boost.org/doc/libs/1_35_0/doc/html/hash/combine.html template
template <class T> inline void hash_combine(size_t &seed, T const &v) {
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

class SparseMatrix{
    public:
    class Point{
        public:
        Point(size_t i, size_t j): i(i), j(j){}

        size_t i, j;
        bool operator<(const Point &other) const {
            if (i < other.i) return true;
            if (other.i < i) return false;
            return j < other.j;
        }
        bool operator==(const Point &other) const{
            return this->i == other.i && this->j == other.j;
        }
    };

    class HashPoint{
        public:
        HashPoint() = default;
        size_t operator()(const SparseMatrix::Point& p) const{
            size_t seed = 0;
            hash_combine(seed, p.i);
            hash_combine(seed, p.j);
            return seed;
        }
    };

    void set(size_t i, size_t j, double val);
    void add(size_t i, size_t j, double val);
    double get(size_t i, size_t j) const;
    void insert_matrix(std::vector<double> M, std::vector<long> pos);
    void merge(SparseMatrix& M);
    std::vector<double> multiply(std::vector<double> vec) const;
    std::vector<double> to_general_band(size_t diag_size, size_t& ku, size_t& kl) const;
    std::vector<size_t> affected_ids(const std::vector<size_t>& ids) const;
    // One-indexed
    void to_mumps_format(std::vector<int>& rows, std::vector<int>& cols, std::vector<double>& vals);
    // Assumes you'll only use to_mumps_format(), so ku/kl are not calculated
    inline void insert_matrix_symmetric_mumps(const math::Matrix& M, const std::vector<long>& pos){
        size_t W = pos.size();
        for(size_t i = 0; i < W; ++i){
            if(pos[i] < 0){
                continue;
            }
            for(size_t j = 0; j <= i; ++j){
                if(pos[j] < 0){
                    continue;
                }
                //if(M[i*W + j] != 0){
                    if(pos[i] >= pos[j]){
                        this->data[Point(pos[i], pos[j])] += M(i, j);
                    } else {
                        this->data[Point(pos[j], pos[i])] += M(i, j);
                    }
                //}
            }
        }
    }
    inline void insert_matrix_general_mumps(const math::Matrix& M, const std::vector<long>& pos){
        size_t W = pos.size();
        for(size_t i = 0; i < W; ++i){
            if(pos[i] < 0){
                continue;
            }
            for(size_t j = 0; j < W; ++j){
                if(pos[j] < 0){
                    continue;
                }
                //if(M[i*W + j] != 0){
                    this->data[Point(pos[i], pos[j])] += M(i, j);
                //}
            }
        }
    }
    inline void insert_block(const math::Matrix& M, const std::vector<long>& pos_i, const std::vector<long>& pos_j, const bool transpose){
        size_t Wi = pos_i.size();
        size_t Wj = pos_j.size();
        if(!transpose){
            for(size_t i = 0; i < Wi; ++i){
                if(pos_i[i] < 0){
                    continue;
                }
                for(size_t j = 0; j < Wj; ++j){
                    if(pos_j[j] < 0){
                        continue;
                    }
                    //if(std::abs(M[i*Wj + j]) > 0){
                        this->data[Point(pos_i[i], pos_j[j])] += M(i, j);
                    //}
                }
            }
        } else {
            for(size_t i = 0; i < Wi; ++i){
                if(pos_i[i] < 0){
                    continue;
                }
                for(size_t j = 0; j < Wj; ++j){
                    if(pos_j[j] < 0){
                        continue;
                    }
                    //if(std::abs(M[i*Wj + j]) > 0){
                        this->data[Point(pos_j[j], pos_i[i])] += M(i, j);
                    //}
                }
            }
        }
    }
    inline auto get_eigen_triplets() const{
        typedef Eigen::Triplet<double, std::ptrdiff_t> T;
        std::vector<T> t;
        t.reserve(this->data.size());
        for(const auto& v:this->data){
            t.emplace_back(v.first.i, v.first.j, v.second);
        }
        return t;
    }
    void zero();
    inline void clear(){this->data.clear();}

    private:
    std::unordered_map<Point, double, HashPoint> data;
    size_t ku = 0;
    size_t kl = 0;
    bool mumps_first_time = true;

    Point point_to_general_band(Point p) const;
};

}

#endif
