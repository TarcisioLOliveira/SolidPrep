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

#ifndef PROJSPEC_REGISTRY_HPP
#define PROJSPEC_REGISTRY_HPP

#include "element_factory.hpp"
#include "logger.hpp"
#include "project_specification/data_map.hpp"
#include "utils.hpp"
#include <functional>
#include <json/value.h>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <list>
#include <vector>

namespace projspec{

namespace object_type{

struct ObjectType{
    std::string name;
    bool array = false;
};

// const static ObjectType SOLID_TYPE{"solid_type"};
// const static ObjectType ANALYSIS{"analysis"};
// const static ObjectType MATERIAL{"material", true}; // ARRAY
// const static ObjectType GEOMETRY{"geometry", true}; // ARRAY
// const static ObjectType FIELDS{"fields", true}; // ARRAY
// const static ObjectType SIMULATION{"simulation"};
// const static ObjectType FINITE_ELEMENT{"finite_element"};
// const static ObjectType SIZING{"sizing"};
// const static ObjectType MESHER{"mesher"};
// const static ObjectType SHAPE_OPT{"shape_opt"};
// const static ObjectType TOPOPT{"topopt"};
// const static ObjectType LOADS{"loads", true}; // ARRAY
// const static ObjectType SUPPORTS{"supports", true}; // ARRAY
// const static ObjectType SPRINGS{"springs", true}; // ARRAY
// const static ObjectType INTERNAL_LOADS{"internal_loads", true}; // ARRAY

}

struct RequirementConditions{
    bool do_meshing = false;
    bool generate_beams = false;
    bool do_topopt = false;
    bool do_fea = false;
    bool do_shape_opt = false;
    bool do_simulation = false;
    utils::ProblemType problem_type = utils::PROBLEM_TYPE_3D;
};

enum DataType{
    TYPE_NULL,
    TYPE_BOOL,
    TYPE_INT,
    TYPE_DOUBLE,
    TYPE_STRING,
    TYPE_ARRAY,
    TYPE_OBJECT,
    TYPE_RELATIVE_PATH,
    TYPE_CROSS_SECTION,
    TYPE_POINTER,
    TYPE_MATRIX,
    TYPE_SHAPE_OP
};

struct ArrayRequirements;

struct ObjectRequirements;

struct MatrixRequirements;

struct DataEntry{
    std::string name;
    DataType type;
    bool required = true;
    std::vector<DataEntry> object_data = std::vector<DataEntry>();
    std::shared_ptr<ArrayRequirements> array_data = nullptr;
    std::function<bool(const RequirementConditions&)> required_if = nullptr;
    std::shared_ptr<MatrixRequirements> matrix_data = nullptr;
};

struct ArrayRequirements{
    size_t size = 0;
    DataType type = TYPE_DOUBLE;
    std::vector<DataEntry> object_data = std::vector<DataEntry>();
    std::shared_ptr<ArrayRequirements> array_data = nullptr;
    std::shared_ptr<MatrixRequirements> matrix_data = nullptr;
};

struct ObjectRequirements{
    std::string type_name;
    std::vector<DataEntry> object_entries;
};

struct MatrixRequirements{
    size_t W;
    size_t H;
};

class Registry{
    public:
    Registry(const Registry&) = delete;
    Registry(Registry&&) = delete;
    Registry& operator=(const Registry&) = delete;
    Registry& operator=(Registry&&) = delete;

    inline static void add(const std::string& type, ObjectRequirements specs){
        get().reqs[type].push_back(std::move(specs));
    }
    inline static bool exists(const std::string& type, const std::string& entry_name){
        const auto& list = get().reqs.at(type);
        for(auto& e:list){
            if(e.type_name == entry_name){
                return true;
            }
        }
        return false;
    }
    inline static bool get(const std::string& type, const std::string& entry_name, ObjectRequirements& entry){
        (void) entry;
        const auto& list = get().reqs.at(type);
        for(auto& e:list){
            if(e.type_name == entry_name){
                entry = e;
                return true;
            }
        }

        return false;
    }

    private:
    static Registry& get(){
        static Registry r;
        return r;
    }

    Registry() = default;
    std::unordered_map<std::string, std::list<ObjectRequirements>> reqs;
    DataEntry empty;
};

class ElementRegistry{
    public:
    ElementRegistry(const ElementRegistry&) = delete;
    ElementRegistry(ElementRegistry&&) = delete;
    ElementRegistry& operator=(const ElementRegistry&) = delete;
    ElementRegistry& operator=(ElementRegistry&&) = delete;

    inline static bool add(const std::string& type, std::unique_ptr<MeshElementFactory> elem){
        get().elems[type] = std::move(elem);
        return true;
    }
    inline static bool exists(const std::string& type){
        return get().elems.count(type) > 0;
    }
    inline static MeshElementFactory* get(const std::string& type){
        auto& instance = get();
        if(instance.exists(type)){
            return instance.elems.at(type).get();
        } else {
            return nullptr;
        }
    }

    private:
    static ElementRegistry& get(){
        static ElementRegistry r;
        return r;
    }

    ElementRegistry() = default;
    std::unordered_map<std::string, std::unique_ptr<MeshElementFactory>> elems;
    DataEntry empty;
};

template <class T>
class Factory{
    public:

    typedef std::function<std::unique_ptr<T>(const projspec::DataMap&)> ObjectMaker;

    Factory() = delete;

    static bool add(ObjectMaker funcCreate, ObjectRequirements reqs){
        auto& subclasses = get();
        auto it = subclasses.find(reqs.type_name);
        if(it == subclasses.end()){
            subclasses[reqs.type_name] = funcCreate;
            Registry::add(T::get_name(), std::move(reqs));
            return true;
        }
        return false;
    }

    static std::unique_ptr<T> construct(const std::string& name, const projspec::DataMap& data){
        auto& subclasses = get();
        auto it = subclasses.find(name);
        if(it != subclasses.end()){
            return it->second(data);
        }

        return nullptr;
    }

    private:
    static std::map<std::string, ObjectMaker>& get(){
        static std::map<std::string, ObjectMaker> subclasses;

        return subclasses;
    }
};

}

#endif
