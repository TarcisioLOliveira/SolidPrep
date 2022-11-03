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

#include <map>
#include <vector>
#include <cstddef>

namespace utils{

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
    };

    void set(size_t i, size_t j, double val);
    void add(size_t i, size_t j, double val);
    double get(size_t i, size_t j) const;
    void insert_matrix(std::vector<double> M, std::vector<long> pos);
    std::vector<double> multiply(std::vector<double> vec) const;
    std::vector<double> to_general_band(size_t diag_size, size_t& ku, size_t& kl) const;
    std::vector<size_t> affected_ids(const std::vector<size_t>& ids) const;
    void to_mumps_format(std::vector<double>& rows, std::vector<double>& cols, std::vector<double>& vals) const;
    inline void clear(){this->data.clear();}

    private:
    std::map<Point, double> data;
    size_t ku = 0;
    size_t kl = 0;

    Point point_to_general_band(Point p) const;
};

}

#endif
