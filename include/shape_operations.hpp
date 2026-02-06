/*
 *   Copyright (C) 2026 Tarcísio Ladeia de Oliveira.
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

#ifndef SHAPE_OPERATIONS_HPP
#define SHAPE_OPERATIONS_HPP

#include <cstddef>
#include <memory>
#include <limits>
#include <json/value.h>

namespace shape_op{

enum class Code{
    UNION,
    INTERSECTION,
    DIFFERENCE,
    GEOMETRY,
    SHELL
};

class ShapeOp;

class JsonRepresentation{
    public:
    JsonRepresentation() = default;
    JsonRepresentation(const Json::Value* doc):doc(doc){}

    inline std::unique_ptr<ShapeOp> construct() const{
        if(this->doc != nullptr){
            return this->get_shape_operations(*(this->doc));
        } else {
            return nullptr;
        }
    }

    private:
    std::unique_ptr<ShapeOp> get_shape_operations(const Json::Value& doc) const;
    size_t get_id(const Json::Value& doc, std::string name) const;
    std::string get_string(const Json::Value& doc, std::string name) const;
    const Json::Value& get_obj(const Json::Value& doc, std::string name) const;

    const Json::Value* doc = nullptr;
};

class ShapeOp{
    public:
    ShapeOp(std::unique_ptr<ShapeOp> s1, std::unique_ptr<ShapeOp> s2):
        s1(std::move(s1)), s2(std::move(s2)){}

    virtual ~ShapeOp() = default;

    virtual Code get_type() const = 0;
    inline ShapeOp* first() const{
        return s1.get();
    }
    inline ShapeOp* second() const{
        return s2.get();
    }
    virtual size_t get_id() const{
        return std::numeric_limits<size_t>::max();
    }

    private:
    std::unique_ptr<ShapeOp> s1, s2;
};

class Union : public ShapeOp{
    public:
    Union(std::unique_ptr<ShapeOp> s1, std::unique_ptr<ShapeOp> s2):
        ShapeOp(std::move(s1), std::move(s2)){}

    virtual ~Union() = default;
    virtual Code get_type() const override{
        return Code::UNION;
    }
};

class Intersection : public ShapeOp{
    public:
    Intersection(std::unique_ptr<ShapeOp> s1, std::unique_ptr<ShapeOp> s2):
        ShapeOp(std::move(s1), std::move(s2)){}

    virtual ~Intersection() = default;
    virtual Code get_type() const override{
        return Code::INTERSECTION;
    }
};

class Difference : public ShapeOp{
    public:
    Difference(std::unique_ptr<ShapeOp> s1, std::unique_ptr<ShapeOp> s2):
        ShapeOp(std::move(s1), std::move(s2)){}

    virtual ~Difference() = default;
    virtual Code get_type() const override{
        return Code::DIFFERENCE;
    }
};

class Geometry : public ShapeOp{
    public:
    Geometry(size_t id):
        ShapeOp(nullptr, nullptr), id(id){}

    virtual ~Geometry() = default;
    virtual Code get_type() const override{
        return Code::GEOMETRY;
    }
    virtual size_t get_id() const override{
        return this->id;
    }

    private:
    const size_t id;
};

class Shell : public ShapeOp{
    public:
    Shell(size_t id):
        ShapeOp(nullptr, nullptr), id(id){}

    virtual ~Shell() = default;
    virtual Code get_type() const override{
        return Code::SHELL;
    }
    virtual size_t get_id() const override{
        return this->id;
    }

    private:
    const size_t id;
};


};

#endif
