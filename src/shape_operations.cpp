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

#include "shape_operations.hpp"
#include "logger.hpp"

namespace shape_op{

std::unique_ptr<ShapeOp> JsonRepresentation::get_shape_operations(const Json::Value& doc) const{
    const std::string type = this->get_string(doc, "type");
    if(type == "geometry"){
        const size_t id = this->get_id(doc, "id");
        return std::make_unique<shape_op::Geometry>(id);
    } else if(type == "union"){
        return std::make_unique<shape_op::Union>(
                    this->get_shape_operations(this->get_obj(doc, "shape1")),
                    this->get_shape_operations(this->get_obj(doc, "shape2"))
                );
    } else if(type == "intersection"){
        return std::make_unique<shape_op::Intersection>(
                    this->get_shape_operations(this->get_obj(doc, "shape1")),
                    this->get_shape_operations(this->get_obj(doc, "shape2"))
                );
    } else if(type == "difference"){
        return std::make_unique<shape_op::Difference>(
                    this->get_shape_operations(this->get_obj(doc, "shape1")),
                    this->get_shape_operations(this->get_obj(doc, "shape2"))
                );
    }

    return nullptr;
}
size_t JsonRepresentation::get_id(const Json::Value& doc, std::string name) const{
    logger::log_assert(doc[name.c_str()].isInt(), logger::ERROR, "Value of key \"{}\" has wrong type, must be a non-negative integer.", name);
    int64_t test = doc[name.c_str()].asInt64();
    logger::log_assert(test >= 0, logger::ERROR, "Value of key \"{}\" has wrong type, must be a non-negative integer.", name);

    return test;
}
std::string JsonRepresentation::get_string(const Json::Value& doc, std::string name) const{
    logger::log_assert(doc[name.c_str()].isString(), logger::ERROR, "Value of key \"{}\" has wrong type, must be a string.", name);
    return doc[name.c_str()].asString();
}
const Json::Value& JsonRepresentation::get_obj(const Json::Value& doc, std::string name) const{
    logger::log_assert(doc[name.c_str()].isObject(), logger::ERROR, "Value of key \"{}\" has wrong type, must be an object.", name);
    return doc[name.c_str()];
}

}
