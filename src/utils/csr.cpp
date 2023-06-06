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

void CSR::zero(){
    for(auto& v:this->data){
        v.second = 0;    
    }
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
    }
}

}
