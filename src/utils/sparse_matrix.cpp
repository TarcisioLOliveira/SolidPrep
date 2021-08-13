/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   Foobar is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Foobar is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with SolidPrep.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "utils/sparse_matrix.hpp"
#include "logger.hpp"

namespace utils{

void SparseMatrix::set(size_t i, size_t j, double val){
    this->data[Point(i, j)] = val;
    size_t k = 0;
    if(j > i){
        k = j - i;
        if(k > this->ku){
            ku = k;
        }
    } else if(i > j){
        k = i - j;
        if(k > this->kl){
            kl = k;
        }
    }
}

double SparseMatrix::get(size_t i, size_t j) const{
    auto pos = this->data.find(Point(i, j));
    if(pos != this->data.end()){
        return pos->second;
    }
    return 0;
}

void SparseMatrix::insert_matrix(std::vector<double> M, std::vector<long> pos){
    size_t W = pos.size();
    for(size_t i = 0; i < pos.size(); ++i){
        for(size_t j = 0; j < pos.size(); ++j){
            if(pos[i] > -1 && pos[j] > -1){
                if(std::abs(M[i*W + j]) > 1e-7){
                    this->set(pos[i], pos[j], M[i*W + j]);
                }
            }
        }
    }
}

std::vector<double> SparseMatrix::multiply(std::vector<double> vec) const{
    std::vector<double> result(vec.size(), 0);
    for(auto& v:this->data){
        size_t i = v.first.i;
        size_t j = v.first.j;
        result[i] = v.second*vec[j];
    }

    return result;
}

std::vector<double> SparseMatrix::to_general_band(size_t diag_size, size_t& ku, size_t& kl) const{
    ku = this->ku;
    kl = this->kl;
    size_t H = 2*kl + ku + 1;
    std::vector<double> band(H*diag_size, 0);
    for(auto& v:this->data){
        Point place = this->point_to_general_band(v.first);
        band[place.i*diag_size + place.j] = v.second;
    }

    return band;
}

SparseMatrix::Point SparseMatrix::point_to_general_band(Point p) const{
    long center = kl + ku;
    long i = 0;
    long j = p.j;
    if(p.j > p.i){
        i = center - (p.j - p.i);
    } else {
        i = center + (p.i - p.j);
    }
    return Point(i, j);
}

}
