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

#ifndef COO_HPP
#define COO_HPP

#include <string>
#include <map>
#include <vector>
#include <cstddef>
#include <limits>
#include "logger.hpp"
#include "math/matrix.hpp"

namespace utils{

template<typename INT>
class COO{
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

    COO(INT offset = 0):
        offset(offset){}

    inline void dot_vector(const std::vector<double>& v, std::vector<double>& v_out, const bool full_matrix) const{
        if(this->coo_first_time){
            if(full_matrix){
                this->map_full_dot_vector(v, v_out);
            } else {
                this->map_sym_dot_vector(v, v_out);
            }
        } else {
            if(full_matrix){
                this->coo_full_dot_vector(v, v_out);
            } else {
                this->coo_sym_dot_vector(v, v_out);
            }
        }
    }

    inline void backup_matrix(){
        this->cooVal_bkp = this->cooVal;
    }

    inline void restore_matrix(){
        std::copy(cooVal_bkp.begin(), cooVal_bkp.end(), cooVal.begin());
    }
    inline void reserve_matrix_symmetric(const math::Matrix& M, const std::vector<long>& pos){
        size_t W = pos.size();
        logger::log_assert(this->coo_first_time, logger::ERROR, "cannot reserve coo space after coo structure is defined.");
        const auto P = this->index_sort(pos);
        size_t start_i = 0;
        while(start_i < W && pos[P[start_i]] < 0){
            ++start_i;
        }
        for(size_t i = start_i; i < W; ++i){
            const auto Pi = pos[P[i]];
            for(size_t j = start_i; j <= i; ++j){
                const auto Pj = pos[P[j]];

                const auto Mij = M(P[i], P[j]);

                // Pi >= Pj
                if(std::abs(Mij) > ZERO_TOL){
                    this->data[Point(Pi, Pj)] += 0;
                }
            }
        }
    }
    inline void reserve_matrix_general(const math::Matrix& M, const std::vector<long>& pos){
        size_t W = pos.size();
        logger::log_assert(this->coo_first_time, logger::ERROR, "cannot reserve coo space after coo structure is defined.");
        const auto P = this->index_sort(pos);
        size_t start_i = 0;
        while(start_i < W && pos[P[start_i]] < 0){
            ++start_i;
        }
        for(size_t i = start_i; i < W; ++i){
            const auto Pi = pos[P[i]];
            for(size_t j = start_i; j < W; ++j){
                const auto Pj = pos[P[j]];

                const auto Mij = M(P[i], P[j]);

                // Pi >= Pj
                if(std::abs(Mij) > ZERO_TOL){
                    this->data[Point(Pi, Pj)] += 0;
                }
            }
        }
    }
    inline void reserve_block(const math::Matrix& M, const std::vector<long>& pos_i, const std::vector<long>& pos_j, const bool transpose){
        size_t Wi = pos_i.size();
        size_t Wj = pos_j.size();
        logger::log_assert(this->coo_first_time, logger::ERROR, "cannot reserve coo space after coo structure is defined.");
        const auto Pi = this->index_sort(pos_i);
        size_t start_i = 0;
        while(start_i < Wi && pos_i[Pi[start_i]] < 0){
            ++start_i;
        }
        const auto Pj = this->index_sort(pos_j);
        size_t start_j = 0;
        while(start_j < Wj && pos_j[Pj[start_j]] < 0){
            ++start_j;
        }
        if(!transpose){
            for(size_t i = start_i; i < Wi; ++i){
                const auto pi = pos_i[Pi[i]];
                for(size_t j = start_i; j < Wj; ++j){
                    const auto pj = pos_j[Pj[j]];

                    const auto Mij = M(Pi[i], Pj[j]);

                    if(std::abs(Mij) > ZERO_TOL){
                        this->data[Point(pi, pj)] += 0;
                    }
                }
            }
        } else {
            for(size_t i = start_i; i < Wi; ++i){
                const auto pi = pos_i[Pi[i]];
                for(size_t j = start_i; j < Wj; ++j){
                    const auto pj = pos_j[Pj[j]];

                    const auto Mij = M(Pi[i], Pj[j]);

                    if(std::abs(Mij) > ZERO_TOL){
                        this->data[Point(pj, pi)] += 0;
                    }
                }
            }
        }
    }

    inline void insert_matrix_symmetric(const math::Matrix& M, const std::vector<long>& pos){
        size_t W = pos.size();
        if(this->coo_first_time){
            const auto P = this->index_sort(pos);
            size_t start_i = 0;
            while(start_i < W && pos[P[start_i]] < 0){
                ++start_i;
            }
            for(size_t i = start_i; i < W; ++i){
                const auto Pi = pos[P[i]];
                for(size_t j = start_i; j <= i; ++j){
                    const auto Pj = pos[P[j]];

                    const auto Mij = M(P[i], P[j]);

                    // Pi >= Pj
                    if(std::abs(Mij) > ZERO_TOL){
                        this->data[Point(Pi, Pj)] += Mij;
                    }
                }
            }
        } else {
            const auto P = this->index_sort(pos);
            for(size_t i = 0; i < W; ++i){
                if(pos[i] > -1){
                    this->coo_add_row(i, pos[i], P, pos, M);
                }
            }
        }
    }
    inline void insert_matrix_general(const math::Matrix& M, const std::vector<long>& pos){
        size_t W = pos.size();
        if(this->coo_first_time){
            const auto P = this->index_sort(pos);
            size_t start_i = 0;
            while(start_i < W && pos[P[start_i]] < 0){
                ++start_i;
            }
            for(size_t i = start_i; i < W; ++i){
                const auto Pi = pos[P[i]];
                for(size_t j = start_i; j < W; ++j){
                    const auto Pj = pos[P[j]];

                    const auto Mij = M(P[i], P[j]);

                    // Pi >= Pj
                    if(std::abs(Mij) > ZERO_TOL){
                        this->data[Point(Pi, Pj)] += Mij;
                    }
                }
            }
        } else {
            const auto P = this->index_sort(pos);
            for(size_t i = 0; i < W; ++i){
                if(pos[i] > -1){
                    this->coo_add_row(i, pos[i], P, pos, M);
                }
            }
        }
    }
    inline void insert_block(const math::Matrix& M, const std::vector<long>& pos_i, const std::vector<long>& pos_j, const bool transpose){
        size_t Wi = pos_i.size();
        size_t Wj = pos_j.size();
        if(this->coo_first_time){
            const auto Pi = this->index_sort(pos_i);
            size_t start_i = 0;
            while(start_i < Wi && pos_i[Pi[start_i]] < 0){
                ++start_i;
            }
            const auto Pj = this->index_sort(pos_j);
            size_t start_j = 0;
            while(start_j < Wj && pos_j[Pj[start_j]] < 0){
                ++start_j;
            }
            if(!transpose){
                for(size_t i = start_i; i < Wi; ++i){
                    const auto pi = pos_i[Pi[i]];
                    for(size_t j = start_i; j < Wj; ++j){
                        const auto pj = pos_j[Pj[j]];

                        const auto Mij = M(Pi[i], Pj[j]);

                        if(std::abs(Mij) > ZERO_TOL){
                            this->data[Point(pi, pj)] += Mij;
                        }
                    }
                }
            } else {
                for(size_t i = start_i; i < Wi; ++i){
                    const auto pi = pos_i[Pi[i]];
                    for(size_t j = start_i; j < Wj; ++j){
                        const auto pj = pos_j[Pj[j]];

                        const auto Mij = M(Pi[i], Pj[j]);

                        if(std::abs(Mij) > ZERO_TOL){
                            this->data[Point(pj, pi)] += Mij;
                        }
                    }
                }
            }
        } else {
            if(!transpose){
                const auto P(this->index_sort(pos_j));
                for(size_t i = 0; i < Wi; ++i){
                    if(pos_i[i] > -1){
                        this->coo_add_row(i, pos_i[i], P, pos_j, M);
                    }
                }
            } else {
                const auto P(this->index_sort(pos_i));
                const math::Matrix MT = M.T();
                for(size_t j = 0; j < Wj; ++j){
                    if(pos_j[j] > -1){
                        this->coo_add_row(j, pos_j[j], P, pos_i, MT);
                    }
                }
            }
        }
    }
    inline void add(const INT i, const INT j, const double v){
        if(coo_first_time){
            this->data[Point(i, j)] += v;
        } else {
            add_coo(i, j, v);
        }
    }
    inline double get(const INT i, const INT j) const{
        if(coo_first_time){
            auto it = this->data.find(Point(i,j));
            if(it != this->data.end()){
                return it->second;
            }
        } else {
            return get_coo(i, j);
        }
        return 0;
    }
    void zero(){
        if(coo_first_time){
            for(auto &v:this->data){
                v.second = 0.0;
            }
        } else {
            std::fill(this->vals.begin(), this->vals.end(), 0.0);
        }
    }

    void generate_coo(INT n);

    const INT& nvals = this->nnz;
    std::vector<INT>& rows = this->cooRowInd;
    std::vector<INT>& cols = this->cooColInd;
    std::vector<double>& vals = this->cooVal;

    inline void clear(){
        this->data.clear();
        this->cooRowPtr.clear();
        this->cooRowInd.clear();
        this->cooColInd.clear();
        this->cooVal.clear();
        this->n = 0;
        this->nnz = 0;
        this->coo_first_time = true;
    }

    void dump_matrix(){
        for(size_t i = 0; i < cooRowPtr.size() - 1; ++i){
            for(INT c = cooRowPtr[i]; c < cooRowPtr[i+1]; ++c){
                std::cout << cooVal[c] << " ";
            }
            std::cout << std::endl;
        }
    }

    private:
    std::map<Point, double> data;
    INT offset = 0;
    INT n = 0;
    INT nnz = 0;
    bool coo_first_time = true;
    std::vector<INT> cooRowPtr;
    std::vector<INT> cooRowInd;
    std::vector<INT> cooColInd;
    std::vector<double> cooVal;
    std::vector<double> cooVal_bkp;
    const double ZERO_TOL = 1e-6;

    inline std::vector<size_t> index_sort(const std::vector<long>& pos) const{
        // https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes

        std::vector<size_t> P(pos.size());
        std::iota(P.begin(), P.end(), 0);
        sort(P.begin(), P.end(), 
            [&](size_t i, size_t j){
                return pos[i] < pos[j];
            }
        );

        return P;
    };

    inline void coo_add_row(const long i, const size_t row, const std::vector<size_t>& P, const std::vector<long>& pos, const math::Matrix& M){
        size_t j = 0;
        while(pos[P[j]] < 0){
            ++j;
        }
        if(j >= pos.size()){
            return;
        }
        for(int c = cooRowPtr[row]; c < cooRowPtr[row+1]; ++c){
            const auto curr = pos[P[j]];
            if(cooColInd[c] == curr + offset){
                cooVal[c] += M(i, P[j]);
                ++j;
                if(j >= pos.size()){
                    break;
                }
            } else if(cooColInd[c] > curr + offset){
                ++j;
                if(j >= pos.size()){
                    break;
                }
                --c;
            }
        }
    } 


    void map_full_dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const;
    void map_sym_dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const;
    void coo_full_dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const;
    void coo_sym_dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const;

    inline void add_coo(const INT i, const INT j, const double v){
        for(INT c = cooRowPtr[i]; c < cooRowPtr[i+1]; ++c){
            if(cooColInd[c] == j + offset){
                cooVal[c] += v;
                break;
            }
        }
    }
    inline double get_coo(const INT i, const INT j) const{
        for(INT c = cooRowPtr[i]; c < cooRowPtr[i+1]; ++c){
            if(j + offset <= cooColInd[c]){
                if(cooColInd[c] == j + offset){
                    return cooVal[c];
                } else {
                    return 0;
                }
            }
        }
        return 0;
    }
};

template<typename INT>
void COO<INT>::map_full_dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const{
    for(auto &p:this->data){
        const size_t i = p.first.i;
        const size_t j = p.first.j;
        const double val = p.second;
        v_out[i] += val*v[j];
    }
}
template<typename INT>
void COO<INT>::map_sym_dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const{
    for(auto &p:this->data){
        const size_t i = p.first.i;
        const size_t j = p.first.j;
        const double val = p.second;
        if(i == j){
            v_out[i] += val*v[j];
        } else {
            v_out[i] += val*v[j];
            v_out[j] += val*v[i];
        }
    }
}
template<typename INT>
void COO<INT>::coo_full_dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const{
    #pragma omp parallel for
    for(size_t i = 0; i < cooRowPtr.size() - 1; ++i){
        for(INT c = cooRowPtr[i]; c < cooRowPtr[i+1]; ++c){
            v_out[i] += cooVal[c]*v[cooColInd[c] - offset];
        }
    }
}
template<typename INT>
void COO<INT>::coo_sym_dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const{
    for(size_t i = 0; i < cooRowPtr.size() - 1; ++i){
        INT c = 0;
        for(c = cooRowPtr[i]; c < cooRowPtr[i+1]-1; ++c){
            const size_t j = cooColInd[c] - offset;
            v_out[i] += cooVal[c]*v[j];
            v_out[j] += cooVal[c]*v[i];
        }
        const size_t j = cooColInd[c] - offset;
        v_out[i] += cooVal[c]*v[j];
    }
}

template<typename INT>
void COO<INT>::generate_coo(INT n){
    if(this->coo_first_time){
        this->coo_first_time = false;
        this->n = n;
        logger::log_assert(this->data.size() < std::numeric_limits<int>::max(), logger::ERROR, "number of non-zero elements is greater than INT_MAX: {}", this->data.size());
        this->nnz = this->data.size();

        this->cooRowPtr.resize(n+1);
        this->cooColInd.resize(nnz);
        this->cooRowInd.resize(nnz);
        this->cooVal.resize(nnz);

        int count = 0;
        auto r = this->cooRowPtr.begin();
        auto ri = this->cooRowInd.begin();
        auto c = this->cooColInd.begin();
        auto v = this->cooVal.begin();
        size_t cur_row = 0;
        *r = 0;
        ++r;
        for(const auto& p:this->data){
            if(p.first.i != cur_row){
                *r = count;
                ++cur_row;
                ++r;
                //if(r >= this->cooRowPtr.end()){
                //    break;
                //}
            }
            *ri = p.first.i+offset;
            *c = p.first.j+offset;
            *v = p.second;
            ++ri;
            ++c;
            ++v;
            ++count;
        }
        *r = nnz;
        this->data.clear();
    }
}

}

#endif
