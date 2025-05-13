/*
 *   Copyright (C) 2025 Tarcísio Ladeia de Oliveira.
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

#ifndef PROJSPEC_DATA_MAP_HPP
#define PROJSPEC_DATA_MAP_HPP

#include <memory>
#include <string>
#include <unordered_map>
#include "cross_section.hpp"

class ProjectData;

namespace projspec{

class DataMap;

class DataArray{
    public:

    inline void set_string(size_t key, std::string value){
        this->string_map[key] = std::move(value);
    }
    inline void set_bool(size_t key, bool value){
        this->bool_map[key] = value;
    }
    inline void set_int(size_t key, int64_t value){
        this->int_map[key] = value;
    }
    inline void set_double(size_t key, double value){
        this->double_map[key] = value;
    }
    inline void set_object(size_t key, std::unique_ptr<DataMap> value){
        this->object_map[key] = std::move(value);
    }
    inline void set_array(size_t key, std::unique_ptr<DataArray> value){
        this->array_map[key] = std::move(value);
    }
    inline void set_cross_section(size_t key, CrossSection value){
        this->cross_section_map[key] = value;
    }

    inline std::string get_string(size_t key, const std::string& none = "") const{
        auto found = this->string_map.find(key);
        if(found != this->string_map.end()){
            return found->second;
        } else {
            return none;
        }
    }
    inline bool get_bool(size_t key, bool none = false) const{
        auto found = this->bool_map.find(key);
        if(found != this->bool_map.end()){
            return found->second;
        } else {
            return none;
        }
    }
    inline int64_t get_int(size_t key, int64_t none = 0) const{
        auto found = this->int_map.find(key);
        if(found != this->int_map.end()){
            return found->second;
        } else {
            return none;
        }
    }
    inline double get_double(size_t key, double none = 0) const{
        auto found = this->double_map.find(key);
        if(found != this->double_map.end()){
            return found->second;
        } else {
            return none;
        }
    }
    DataMap* get_object(size_t key, DataMap* none = nullptr) const;
    inline DataArray* get_array(size_t key, DataArray* none = nullptr) const{
        static DataArray empty_array;

        if(none == nullptr){
            none = &empty_array;
        }
        auto found = this->array_map.find(key);
        if(found != this->array_map.end()){
            return found->second.get();
        } else {
            return none;
        }
    }
    inline CrossSection get_cross_section(size_t key, CrossSection none = CrossSection()) const{
        auto found = this->cross_section_map.find(key);
        if(found != this->cross_section_map.end()){
            return found->second;
        } else {
            return none;
        }
    }

    inline std::vector<std::string> get_string_array() const{
        std::vector<std::string> array(this->string_map.size());
        for(size_t i = 0; i < array.size(); ++i){
            array[i] = this->string_map.at(i);
        }
        return array;
    }
    inline std::vector<bool> get_bool_array() const{
        std::vector<bool> array(this->bool_map.size());
        for(size_t i = 0; i < array.size(); ++i){
            array[i] = this->bool_map.at(i);
        }
        return array;
    }
    inline std::vector<int64_t> get_int_array() const{
        std::vector<int64_t> array(this->int_map.size());
        for(size_t i = 0; i < array.size(); ++i){
            array[i] = this->int_map.at(i);
        }
        return array;
    }
    inline std::vector<double> get_double_array() const{
        std::vector<double> array(this->double_map.size());
        for(size_t i = 0; i < array.size(); ++i){
            array[i] = this->double_map.at(i);
        }
        return array;
    }

    template<size_t N>
    inline std::array<std::string, N> get_string_array_fixed() const{
        std::array<std::string, N> array;
        for(size_t i = 0; i < array.size(); ++i){
            array[i] = this->string_map.at(i);
        }
        return array;
    }
    template<size_t N>
    inline std::array<bool, N> get_bool_array_fixed() const{
        std::array<bool, N> array;
        for(size_t i = 0; i < array.size(); ++i){
            array[i] = this->bool_map.at(i);
        }
        return array;
    }
    template<size_t N>
    inline std::array<int64_t, N> get_int_array_fixed() const{
        std::array<int64_t, N> array;
        for(size_t i = 0; i < array.size(); ++i){
            array[i] = this->int_map.at(i);
        }
        return array;
    }
    template<size_t N>
    inline std::array<double, N> get_double_array_fixed() const{
        std::array<double, N> array;
        for(size_t i = 0; i < array.size(); ++i){
            array[i] = this->double_map.at(i);
        }
        return array;
    }

    size_t size() const{
        return std::max({
            string_map.size(),
            bool_map.size(),
            int_map.size(),
            double_map.size(),
            object_map.size(),
            array_map.size()
        });
    }
    private:
    std::unordered_map<size_t, std::string> string_map;
    std::unordered_map<size_t, bool> bool_map;
    std::unordered_map<size_t, int64_t> int_map;
    std::unordered_map<size_t, double> double_map;
    std::unordered_map<size_t, std::unique_ptr<DataMap>> object_map;
    std::unordered_map<size_t, std::unique_ptr<DataArray>> array_map;
    std::unordered_map<size_t, CrossSection> cross_section_map;
};

class DataMap{
    public:
    DataMap(ProjectData* proj)
        : proj(proj)
    {}

    ProjectData* proj;
    
    inline void set_string(const std::string& key, std::string value){
        this->string_map[key] = std::move(value);
    }
    inline void set_bool(const std::string& key, bool value){
        this->bool_map[key] = value;
    }
    inline void set_int(const std::string& key, int64_t value){
        this->int_map[key] = value;
    }
    inline void set_double(const std::string& key, double value){
        this->double_map[key] = value;
    }
    inline void set_object(const std::string& key, std::unique_ptr<DataMap> value){
        this->object_map[key] = std::move(value);
    }
    inline void set_array(const std::string& key, std::unique_ptr<DataArray> value){
        this->array_map[key] = std::move(value);
    }
    inline void set_cross_section(const std::string& key, CrossSection value){
        this->cross_section_map[key] = value;
    }

    inline bool exists_string(const std::string& key) const{
        return this->string_map.contains(key);
    }
    inline bool exists_bool(const std::string& key) const{
        return this->bool_map.contains(key);
    }
    inline bool exists_int(const std::string& key) const{
        return this->int_map.contains(key);
    }
    inline bool exists_double(const std::string& key) const{
        return this->double_map.contains(key);
    }
    inline bool exists_object(const std::string& key) const{
        return this->object_map.contains(key);
    }
    inline bool exists_array(const std::string& key) const{
        return this->array_map.contains(key);
    }
    inline bool exists_cross_section(const std::string& key) const{
        return this->cross_section_map.contains(key);
    }

    inline std::string get_string(const std::string& key, const std::string& none = "") const{
        auto found = this->string_map.find(key);
        if(found != this->string_map.end()){
            return found->second;
        } else {
            return none;
        }
    }
    inline bool get_bool(const std::string& key, bool none = false) const{
        auto found = this->bool_map.find(key);
        if(found != this->bool_map.end()){
            return found->second;
        } else {
            return none;
        }
    }
    inline int64_t get_int(const std::string& key, int64_t none = 0) const{
        auto found = this->int_map.find(key);
        if(found != this->int_map.end()){
            return found->second;
        } else {
            return none;
        }
    }
    inline double get_double(const std::string& key, double none = 0) const{
        auto found = this->double_map.find(key);
        if(found != this->double_map.end()){
            return found->second;
        } else {
            return none;
        }
    }
    inline DataMap* get_object(const std::string& key, DataMap* none = nullptr) const{
        static DataMap empty_map(nullptr);

        if(none == nullptr){
            none = &empty_map;
        }
        auto found = this->object_map.find(key);
        if(found != this->object_map.end()){
            return found->second.get();
        } else {
            return none;
        }
    }
    inline DataArray* get_array(const std::string& key, DataArray* none = nullptr) const{
        static DataArray empty_array;

        if(none == nullptr){
            none = &empty_array;
        }
        auto found = this->array_map.find(key);
        if(found != this->array_map.end()){
            return found->second.get();
        } else {
            return none;
        }
    }
    inline CrossSection get_cross_section(const std::string& key, CrossSection none = CrossSection()) const{
        auto found = this->cross_section_map.find(key);
        if(found != this->cross_section_map.end()){
            return found->second;
        } else {
            return none;
        }
    }
    private:
    std::unordered_map<std::string, std::string> string_map;
    std::unordered_map<std::string, bool> bool_map;
    std::unordered_map<std::string, int64_t> int_map;
    std::unordered_map<std::string, double> double_map;
    std::unordered_map<std::string, std::unique_ptr<DataMap>> object_map;
    std::unordered_map<std::string, std::unique_ptr<DataArray>> array_map;
    std::unordered_map<std::string, CrossSection> cross_section_map;
};

}

#endif
