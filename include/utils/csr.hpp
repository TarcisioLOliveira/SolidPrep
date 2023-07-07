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

#ifndef CSR_HPP
#define CSR_HPP

#include <string>
#include <map>
#include <vector>
#include <cstddef>

namespace utils{

class CSR{
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

    inline void insert_matrix_symmetric(const std::vector<double>& M, const std::vector<long>& pos){
        size_t W = pos.size();
        if(this->csr_first_time){
            for(size_t i = 0; i < W; ++i){
                for(size_t j = 0; j <= i; ++j){
                    if(M[i*W + j] != 0){
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
                if(pos[i] < 0){
                    continue;
                }
                size_t j = 0;
                while(j <= i && (pos[j] < 0 || M[i*W+j] == 0)){
                    ++j;
                }
                // pos is not ordered!
                while(j <= i){
                    bool found_one = false;
                    for(int c = csrRowPtr[pos[i]]; c < csrRowPtr[pos[i]+1]; ++c){
                        if(csrColInd[c] == pos[j]){
                            found_one = true;
                            csrVal[c] += M[i*W+j];
                            do{
                                ++j;
                            } while(j <= i && (pos[j] < 0 || M[i*W+j] == 0));
                            if(j > i){
                                break;
                            }
                        }
                    }
                    if(!found_one){
                        do{
                            ++j;
                        } while(j <= i && (pos[j] < 0 || M[i*W+j] == 0));
                    }
                }
            }
        }
    }
    void zero();

    void generate_csr(int n);

    const int& nvals = this->nnz;
    const std::vector<int>& rows = this->csrRowPtr;
    const std::vector<int>& cols = this->csrColInd;
    const std::vector<double>& vals = this->csrVal;

    inline void clear(){this->data.clear();}

    private:
    std::map<Point, double> data;
    int n = 0;
    int nnz = 0;
    bool csr_first_time = true;
    std::vector<int> csrRowPtr;
    std::vector<int> csrColInd;
    std::vector<double> csrVal;
    
    inline void add_csr(const int i, const int j, const double v){
        for(int c = csrRowPtr[i]; c < csrRowPtr[i+1]; ++c){
            if(csrColInd[c] == j){
                csrVal[c] += v;
                break;
            }
        }
    }
};

}

#endif
