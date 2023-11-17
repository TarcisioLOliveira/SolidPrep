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

    inline void insert_matrix_symmetric(const std::vector<double>& M, const std::vector<long>& pos){
        size_t W = pos.size();
        if(this->coo_first_time){
            for(size_t i = 0; i < W; ++i){
                for(size_t j = 0; j <= i; ++j){
                    if(std::abs(M[i*W + j]) > 0){
                        if(pos[i] > -1 && pos[j] > -1){
                            if(pos[i] >= pos[j]){
                                this->data[Point(pos[i], pos[j])] += M[i*W + j];
                            } else {
                                this->data[Point(pos[j], pos[i])] += M[i*W + j];
                            }
                        }
                    }
                }
            }
        } else {
            for(size_t i = 0; i < W; ++i){
                for(size_t j = 0; j <= i; ++j){
                    if(pos[i] > -1 && pos[j] > -1){
                        if(pos[i] >= pos[j]){
                            this->add_coo(pos[i], pos[j], M[i*W + j]);
                        } else {
                            this->add_coo(pos[j], pos[i], M[i*W + j]);
                        }
                    }
                }
            }
            //for(size_t i = 0; i < W; ++i){
            //    if(pos[i] < 0){
            //        continue;
            //    }
            //    size_t j = 0;
            //    while(j <= i && (pos[j] < 0 || M[i*W+j] == 0)){
            //        ++j;
            //    }
            //    // pos is not ordered!
            //    while(j <= i){
            //        bool found_one = false;
            //        for(int c = cooRowPtr[pos[i]]; c < cooRowPtr[pos[i]+1]; ++c){
            //            if(cooColInd[c] == pos[j] + offset){
            //                found_one = true;
            //                cooVal[c] += M[i*W+j];
            //                do{
            //                    ++j;
            //                } while(j <= i && (pos[j] < 0 || M[i*W+j] == 0));
            //                if(j > i){
            //                    break;
            //                }
            //            }
            //        }
            //        if(!found_one){
            //            do{
            //                ++j;
            //            } while(j <= i && (pos[j] < 0 || M[i*W+j] == 0));
            //        }
            //    }
            //}
        }
    }
    inline void insert_matrix_general(const std::vector<double>& M, const std::vector<long>& pos){
        size_t W = pos.size();
        if(this->coo_first_time){
            for(size_t i = 0; i < W; ++i){
                for(size_t j = 0; j < W; ++j){
                    if(std::abs(M[i*W + j]) > 0){
                        if(pos[i] > -1 && pos[j] > -1){
                            this->data[Point(pos[i], pos[j])] += M[i*W + j];
                        }
                    }
                }
            }
        } else {
            for(size_t i = 0; i < W; ++i){
                if(pos[i] < 0){
                    continue;
                }
                size_t j = 0;
                while(j < W && pos[j] < 0){
                    ++j;
                }
                // pos is not ordered!
                while(j < W){
                    bool found_one = false;
                    for(int c = cooRowPtr[pos[i]]; c < cooRowPtr[pos[i]+1]; ++c){
                        if(cooColInd[c] == pos[j] + offset){
                            found_one = true;
                            cooVal[c] += M[i*W+j];
                            size_t j_prev = j;
                            do{
                                ++j;
                            } while(j < W && pos[j] < 0);
                            if(j >= W || pos[j] < pos[j_prev]){
                                break;
                            }
                        }
                    }
                    if(!found_one){
                        do{
                            ++j;
                        } while(j < W && pos[j] < 0);
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
    
    inline void add_coo(const INT i, const INT j, const double v){
        for(int c = cooRowPtr[i]; c < cooRowPtr[i+1]; ++c){
            if(cooColInd[c] == j + offset){
                cooVal[c] += v;
                break;
            }
        }
    }
};

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
