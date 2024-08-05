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

    inline bool generate_hessian(const size_t hoffset, const size_t l_len, std::vector<double>& l, const std::vector<double>& Ku, const bool full_matrix){
        if(full_matrix){
            return this->coo_full_hessian(hoffset, l_len, l, Ku);
        } else {
            return this->coo_sym_hessian(hoffset, l_len, l, Ku);
        }
    }

    inline void insert_matrix_symmetric(const std::vector<double>& M, const std::vector<long>& pos){
        size_t W = pos.size();
        if(this->coo_first_time){
            for(size_t i = 0; i < W; ++i){
                if(pos[i] < 0){
                    continue;
                }
                for(size_t j = 0; j <= i; ++j){
                    if(pos[j] < 0){
                        continue;
                    }
                    if(std::abs(M[i*W + j]) > 0){
                        if(pos[i] >= pos[j]){
                            this->data[Point(pos[i], pos[j])] += M[i*W + j];
                        } else {
                            this->data[Point(pos[j], pos[i])] += M[i*W + j];
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
                if(pos[i] < 0){
                    continue;
                }
                for(size_t j = 0; j < W; ++j){
                    if(pos[j] < 0){
                        continue;
                    }
                    if(std::abs(M[i*W + j]) > 0){
                        this->data[Point(pos[i], pos[j])] += M[i*W + j];
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
    inline void insert_block(const std::vector<double>& M, const std::vector<long>& pos_i, const std::vector<long>& pos_j, const bool transpose){
        size_t Wi = pos_i.size();
        size_t Wj = pos_j.size();
        if(this->coo_first_time){
            if(!transpose){
                for(size_t i = 0; i < Wi; ++i){
                    if(pos_i[i] < 0){
                        continue;
                    }
                    for(size_t j = 0; j < Wj; ++j){
                        if(pos_j[j] < 0){
                            continue;
                        }
                        if(std::abs(M[i*Wj + j]) > 0){
                            this->data[Point(pos_i[i], pos_j[j])] += M[i*Wj + j];
                        }
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
                        if(std::abs(M[i*Wj + j]) > 0){
                            this->data[Point(pos_j[j], pos_i[i])] += M[i*Wj + j];
                        }
                    }
                }
            }
        } else {
            if(!transpose){
                for(size_t i = 0; i < Wi; ++i){
                    if(pos_i[i] < 0){
                        continue;
                    }
                    size_t j = 0;
                    while(j < Wj && pos_j[j] < 0){
                        ++j;
                    }
                    // pos is not ordered!
                    while(j < Wj){
                        bool found_one = false;
                        for(int c = cooRowPtr[pos_i[i]]; c < cooRowPtr[pos_i[i]+1]; ++c){
                            if(cooColInd[c] == pos_j[j] + offset){
                                found_one = true;
                                cooVal[c] += M[i*Wj+j];
                                size_t j_prev = j;
                                do{
                                    ++j;
                                } while(j < Wj && pos_j[j] < 0);
                                if(j >= Wj || pos_j[j] < pos_j[j_prev]){
                                    break;
                                }
                            }
                        }
                        if(!found_one){
                            do{
                                ++j;
                            } while(j < Wj && pos_j[j] < 0);
                        }
                    }
                }
            } else {
                for(size_t j = 0; j < Wj; ++j){
                    if(pos_j[j] < 0){
                        continue;
                    }
                    size_t i = 0;
                    while(i < Wi && pos_i[i] < 0){
                        ++i;
                    }
                    // pos is not ordered!
                    while(i < Wi){
                        bool found_one = false;
                        for(int c = cooRowPtr[pos_j[j]]; c < cooRowPtr[pos_j[j]+1]; ++c){
                            if(cooColInd[c] == pos_i[i] + offset){
                                found_one = true;
                                cooVal[c] += M[i*Wj+j];
                                size_t i_prev = i;
                                do{
                                    ++i;
                                } while(i < Wi && pos_i[i] < 0);
                                if(i >= Wi || pos_i[i] < pos_i[i_prev]){
                                    break;
                                }
                            }
                        }
                        if(!found_one){
                            do{
                                ++i;
                            } while(i < Wi && pos_i[i] < 0);
                        }
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


    void map_full_dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const;
    void map_sym_dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const;
    void coo_full_dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const;
    void coo_sym_dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const;

    bool coo_full_hessian(const size_t hoffset, const size_t l_len, std::vector<double>& l, const std::vector<double>& Ku);
    bool coo_sym_hessian(const size_t hoffset, const size_t l_len, std::vector<double>& l, const std::vector<double>& Ku);
    
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
bool COO<INT>::coo_full_hessian(const size_t hoffset, const size_t l_len, std::vector<double>& l, const std::vector<double>& Ku){
    // hoffset = u_size + 2*l_len
    //bool modified = false;
    //for(size_t i = hoffset; i < cooRowPtr.size() - 1; ++i){
    //    INT c = cooRowPtr[i+1];
    //    double& li = l[2*l_len + i - hoffset];
    //    const double dli = 2*li;
    //    const double ddli = 2;
    //    while(c < cooRowPtr[i+1]){
    //        const size_t j = cooColInd[c] - offset;
    //        if(i == j){
    //            break;
    //        }
    //        ++c;
    //    }
    //    const double diag = cooVal[c]*dli*dli + ddli*Ku[i];
    //    if(diag <= 1e-9 && li != 0){
    //        li = 0;
    //        modified = true;
    //    }
    //}
    //if(modified){
    //    return true;
    //}

    // Shift modification
    double min_diag = 0;
    for(size_t i = hoffset; i < cooRowPtr.size() - 1; ++i){
        INT c = cooRowPtr[i];
        double& li = l[2*l_len + i - hoffset];
        const double dli = 2*li;
        const double ddli = 2;
        while(c < cooRowPtr[i+1]){
            const size_t j = cooColInd[c] - offset;
            if(i == j){
                break;
            }
            ++c;
        }
        const double diag = cooVal[c]*dli*dli + ddli*Ku[i];
        if(diag < 0 && li != 0 && diag < min_diag){
            min_diag = diag;
        }
    }
    min_diag *= 1.01;
    logger::quick_log("min diag", min_diag);

    // "Vertical" step
    for(size_t i = 0; i < hoffset; ++i){
        INT c = cooRowPtr[i];
        while(cooColInd[c] - offset < hoffset && c < cooRowPtr[i+1]){
            //const size_t j = cooColInd[c] - offset;
            //if(i == j){
            //    cooVal[c] -= min_diag;
            //}
            ++c;
        }
        while(c < cooRowPtr[i+1]){
            const size_t j = cooColInd[c] - offset;
            const double lj = l[2*l_len + j - hoffset];
            const double dlj = -1.0/(lj*lj);//2*lj;
            cooVal[c] *= dlj;
            ++c;
        }
    }
    // "Horizontal" step
    for(size_t i = hoffset; i < cooRowPtr.size() - 1; ++i){
        INT c = cooRowPtr[i];
        const double li = l[2*l_len + i - hoffset];
        const double dli = -1.0/(li*li);//2*li;
        while(cooColInd[c] - offset < hoffset && c < cooRowPtr[i+1]){
            cooVal[c] *= dli;
            ++c;
        }
        // "Square" step
        while(c < cooRowPtr[i+1]){
            const size_t j = cooColInd[c] - offset;
            const double lj = l[2*l_len + j - hoffset];
            const double dlj = -1.0/(lj*lj);//2*lj;
            cooVal[c] *= dlj*dli;

            // "Sum" step
            if(i == j){
                const double ddli = 2;
                cooVal[c] += ddli*Ku[i];// + 1e-9;
                cooVal[c] -= min_diag;
                cooVal[c] += 1e-9;
                //if(cooVal[c] <= 1e-9){// || std::isnan(cooVal[c]) || std::isinf(cooVal[c])){
                //    cooVal[c] = 1e-9;
                //    //logger::quick_log(cooVal[c], ddli, li, Ku[i], i, j);
                //}
            }
            ++c;
        }
    }
    return false;
}

template<typename INT>
bool COO<INT>::coo_sym_hessian(const size_t hoffset, const size_t l_len, std::vector<double>& l, const std::vector<double>& Ku){
    //bool modified = false;
    //for(size_t i = hoffset; i < cooRowPtr.size() - 1; ++i){
    //    INT c = cooRowPtr[i+1];
    //    double& li = l[2*l_len + i - hoffset];
    //    const double dli = 2*li;
    //    const double ddli = 2;
    //    while(c < cooRowPtr[i+1]){
    //        const size_t j = cooColInd[c] - offset;
    //        if(i == j){
    //            break;
    //        }
    //        ++c;
    //    }
    //    const double diag = cooVal[c]*dli*dli + ddli*Ku[i];
    //    if(diag <= 1e-9 && li != 0){
    //        li = 0;
    //        modified = true;
    //    }
    //}
    //if(modified){
    //    return true;
    //}
    // Shift modification
    double min_diag = 0;
    for(size_t i = hoffset; i < cooRowPtr.size() - 1; ++i){
        INT c = cooRowPtr[i];
        double& li = l[2*l_len + i - hoffset];
        const double dli = 2*li;
        const double ddli = 2;
        while(c < cooRowPtr[i+1]){
            const size_t j = cooColInd[c] - offset;
            if(i == j){
                break;
            }
            ++c;
        }
        const double diag = cooVal[c]*dli*dli + ddli*Ku[i];
        if(diag < 0 && li != 0 && diag < min_diag){
            min_diag = diag;
        }
    }
    min_diag *= 1.01;
    logger::quick_log("min diag", min_diag);
    // "Horizontal" step
    for(size_t i = hoffset; i < cooRowPtr.size() - 1; ++i){
        INT c = cooRowPtr[i];
        const double li = l[2*l_len + i - hoffset];
        const double dli = 2*li;//li/std::sqrt(li*li + 1e-14); // Approx li, except near 0, then 0
        while(cooColInd[c] - offset < hoffset && c < cooRowPtr[i+1]){
            cooVal[c] *= dli;
            ++c;
        }
        // "Square" step
        while(c < cooRowPtr[i+1]){
            const size_t j = cooColInd[c] - offset;
            const double lj = l[2*l_len + j - hoffset];
            const double dlj = 2*lj;//lj/std::sqrt(lj*lj + 1e-14);
            cooVal[c] *= dlj*dli;

            ++c;
        }
        const double ddli = 2;//1e-14/std::pow(li*li + 1e-14, 1.5); // Approx 0
        // "Sum" step
        cooVal[c - 1] += ddli*Ku[i];
        cooVal[c - 1] -= min_diag;
        cooVal[c - 1] += 1e-9;
        //std::cout << cooVal[c - 1] << " ";
    }
    //std::cout << std::endl;
    return false;
}

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
