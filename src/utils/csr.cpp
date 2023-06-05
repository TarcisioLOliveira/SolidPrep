/*
 *   Copyright (C) 2023 Tarcísio Ladeia de Oliveira.
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

#include "utils/csr.hpp"
#include "logger.hpp"
#include <algorithm>
#include <limits>
#include <set>

namespace utils{

void CSR::set(size_t i, size_t j, double val){
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

void CSR::add(size_t i, size_t j, double val){
    this->data[Point(i, j)] += val;
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

double CSR::get(size_t i, size_t j) const{
    auto pos = this->data.find(Point(i, j));
    if(pos != this->data.end()){
        return pos->second;
    }
    return 0;
}

void CSR::insert_matrix(std::vector<double> M, std::vector<long> pos){
    size_t W = pos.size();
    for(size_t i = 0; i < pos.size(); ++i){
        for(size_t j = 0; j < pos.size(); ++j){
            if(pos[i] > -1 && pos[j] > -1){
                if(std::abs(M[i*W + j]) > 1e-7){
                    this->add(pos[i], pos[j], M[i*W + j]);
                }
            }
        }
    }
}

std::vector<double> CSR::multiply(std::vector<double> vec) const{
    std::vector<double> result(vec.size(), 0);
    for(auto& v:this->data){
        size_t i = v.first.i;
        size_t j = v.first.j;
        result[i] = v.second*vec[j];
    }

    return result;
}

std::vector<double> CSR::to_general_band(size_t diag_size, size_t& ku, size_t& kl) const{
    ku = this->ku;
    kl = this->kl;
    size_t H = kl + ku + 1;
    std::vector<double> band(H*diag_size, 0);
    for(auto& v:this->data){
        Point place = this->point_to_general_band(v.first);
        band[place.i*diag_size + place.j] = v.second;
    }

    return band;
}

void CSR::zero(){
    for(auto& v:this->data){
        v.second = 0;    
    }
}

void CSR::merge(CSR& M){
    //this->data.reserve(this->data.size() + M.data.size());
    this->data.merge(M.data);
    for(const auto& v:M.data){
        this->data[v.first] += v.second;
    }
}
void CSR::to_mumps_format(std::vector<int>& rows, std::vector<int>& cols, std::vector<double>& vals) const{
    size_t N = this->data.size();
    if(N > rows.size()){
        rows.resize(N, 0);
        cols.resize(N, 0);
        vals.resize(N, 0);
    }
    auto v = this->data.cbegin();
    for(size_t i = 0; i < N; ++i){
        rows[i] = v->first.i+1;
        cols[i] = v->first.j+1;
        vals[i] = v->second;
        ++v;
    }
}

CSR::Point CSR::point_to_general_band(Point p) const{
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

std::vector<size_t> CSR::affected_ids(const std::vector<size_t>& ids) const{
    std::set<size_t> affected;
    size_t cur_i = this->data.begin()->first.i;
    size_t cur_id = 0;
    bool skip_line = false;
    for(auto& v:this->data){
        if(cur_i != v.first.i){
            cur_i = v.first.i;
            cur_id = 0;
            skip_line = false;
        }
        if(skip_line){
            continue;
        }
        if(v.first.j < ids[cur_id]){
            continue;
        }
        if(v.first.j > ids[cur_id]){
            ++cur_id;
        }
        if(v.first.j == ids[cur_id]){
            if(v.second != 0){
                affected.insert(v.first.i);
                skip_line = true;
            }
            ++cur_id;
        }
    }
    std::vector<size_t> aff(affected.begin(), affected.end());

    return aff;
}

void CSR::generate_csr(int n){
    if(this->csr_first_time){
        this->csr_first_time = false;
        this->n = n;
        logger::log_assert(this->data.size() < std::numeric_limits<int>::max(), logger::ERROR, "number of non-zero elements is greater than INT_MAX: {}", this->data.size());
        this->nnz = this->data.size();
        this->csrRowPtr.resize(n+1);
        this->csrColInd.resize(nnz);
        this->csrVal.resize(nnz);

        int count = 0;
        auto r = this->csrRowPtr.begin();
        auto c = this->csrColInd.begin();
        auto v = this->csrVal.begin();
        size_t cur_row = 0;
        *r = 0;
        ++r;
        for(const auto& p:this->data){
            if(p.first.i != cur_row){
                *r = count;
                ++cur_row;
                ++r;
            }
            *c = p.first.j;
            *v = p.second;

            ++c;
            ++v;
            ++count;
        }
        *r = nnz;
    } else {
        auto v = this->csrVal.begin();
        for(auto& p:this->data){
            *v = p.second;
            ++v;
        }
    }
}

}
