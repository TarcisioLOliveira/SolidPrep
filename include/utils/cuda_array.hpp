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

#ifndef CUDA_ARRAY_HPP
#define CUDA_ARRAY_HPP

#include "logger.hpp"
#include <vector>
#include <cuda_runtime.h>


namespace utils::cuda{

template<typename T>
class CUDAArray{
    public:
    CUDAArray() = default;
    ~CUDAArray(){
        if(this->data != nullptr){
            cudaFree(this->data);
        }
    }
    CUDAArray(const CUDAArray&) = delete;
    CUDAArray(CUDAArray&& c):data(c.data), size(c.size){
        c.data = nullptr;
    }

    CUDAArray(size_t size):size(size){
        this->check_error(cudaMalloc(&this->data, sizeof(T)*size));
    }
    CUDAArray(size_t size, T init):size(size){
        this->check_error(cudaMalloc(&this->data, sizeof(T)*size));
        this->check_error(cudaMemset(this->data, init, sizeof(T)*size));
    }
    CUDAArray(const std::vector<T>& host_data):size(host_data.size()){
        this->check_error(cudaMalloc(&this->data, sizeof(T)*size));
        this->check_error(cudaMemcpy(static_cast<void*>(this->data), static_cast<const void*>(host_data.data()), sizeof(T)*this->size, cudaMemcpyHostToDevice));
    }
    CUDAArray(T* host_data, size_t size):size(size){
        this->check_error(cudaMalloc(&this->data, sizeof(T)*size));
        this->check_error(cudaMemcpy(static_cast<void*>(this->data), static_cast<void*>(host_data), sizeof(T)*this->size, cudaMemcpyHostToDevice));
    }

    inline void allocate(size_t size){
        if(this->data != nullptr){
            cudaFree(this->data);
            this->data = nullptr;
        }
        this->size = size;
        this->check_error(cudaMalloc(&this->data, sizeof(T)*size));
    }
    inline void allocate(size_t size, T init){
        if(this->data != nullptr){
            cudaFree(this->data);
            this->data = nullptr;
        }
        this->size = size;
        this->check_error(cudaMalloc(&this->data, sizeof(T)*size));
        this->check_error(cudaMemset(this->data, init, sizeof(T)*size));
    }
    inline void allocate(const std::vector<T>& host_data){
        if(this->data != nullptr){
            cudaFree(this->data);
            this->data = nullptr;
        }
        this->size = host_data.size();
        this->check_error(cudaMalloc(&this->data, sizeof(T)*this->size));
        this->check_error(cudaMemcpy(static_cast<void*>(this->data), static_cast<const void*>(host_data.data()), sizeof(T)*this->size, cudaMemcpyHostToDevice));
    }

    inline T* pointer(){
        return this->data;
    }
    inline T* pointer() const{
        return this->data;
    }
    inline size_t get_size() const{
        return this->size;
    }

    inline void set_data(const std::vector<T>& host_data, size_t range){
        this->check_error(cudaMemcpy(static_cast<void*>(this->data), static_cast<const void*>(host_data.data()), sizeof(T)*range, cudaMemcpyHostToDevice));
    }
    inline void set_data(const std::vector<T>& host_data){
        this->check_error(cudaMemcpy(static_cast<void*>(this->data), static_cast<const void*>(host_data.data()), sizeof(T)*this->size, cudaMemcpyHostToDevice));
    }
    inline void get_data(std::vector<T>& host_data) const{
        this->check_error(cudaMemcpy(static_cast<void*>(host_data.data()), static_cast<void*>(this->data), sizeof(T)*this->size, cudaMemcpyDeviceToHost));
    }
    inline void get_data(std::vector<T>& host_data, size_t range) const{
        this->check_error(cudaMemcpy(static_cast<void*>(host_data.data()), static_cast<void*>(this->data), sizeof(T)*range, cudaMemcpyDeviceToHost));
    }
    inline void set_zero(){
        if(this->data != nullptr){
            this->check_error(cudaMemset(this->data, 0.0, sizeof(T)*size));
        }
    }

    inline CUDAArray& operator=(const CUDAArray&) = delete;
    inline CUDAArray& operator=(CUDAArray&& c){
        if(this->data != nullptr){
            cudaFree(this->data);
        }
        this->data = c.data;
        c.data = nullptr;
        this->size = c.size;

        return *this;
    }

    private:
    T* data = nullptr;
    size_t size = 0;

    inline void check_error(cudaError_t error) const{
        logger::log_assert(error == cudaSuccess, logger::ERROR, "CUDA runtime call returned {} in CUDAArray instance.", cudaGetErrorString(error));
    }
};

}

#endif
