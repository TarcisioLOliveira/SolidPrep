/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   Foobar is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Foobar is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with SolidPrep.  If not, see <https://www.gnu.org/licenses/>.
 *
 */
#ifndef PROJECT_DATA_HPP
#define PROJECT_DATA_HPP

#include <string>
#include <vector>
#include "rapidjson/document.h"
#include "STEPCAFControl_Reader.hxx"
#include "logger.hpp"
#include "force.hpp"
#include "support.hpp"

/**
 * Reads and stores project data.
 */
class ProjectData {
    public:

    enum ProjectType{
        TYPE_2D,
        TYPE_3D
    };

    enum DataType{
        TYPE_NULL,
        TYPE_BOOL,
        TYPE_INT,
        TYPE_DOUBLE,
        TYPE_STRING,
        TYPE_ARRAY,
        TYPE_OBJECT
    };

    /**
     * Loads project data file.
     *
     * @param project_file Path to file.
     */
    ProjectData(std::string project_file);

    double step;
    std::vector<double> angles2D;
    std::vector<std::array<double,2>> angles3D;
    double restriction;
    double scale;
    double thickness;
    ProjectType type;
    TopoDS_Shape solid;
    std::vector<Force> forces;
    std::vector<Support> supports;
    
    private:
    /**
     * Checks for existence and type of data in JSON file.
     *
     * @param doc Document being used.
     * @param name JSON key.
     * @param type Data type to be checked for.
     * @param required Whether the parameter is required or optional;
     *
     * @return Whether the key exists and has the correct type.
     */
    template<typename A, typename B>
    bool log_data(const rapidjson::GenericValue<A, B>& doc, std::string name, DataType type, bool required) const;
};

template<typename A, typename B>
bool ProjectData::log_data(const rapidjson::GenericValue<A, B>& doc, std::string name, ProjectData::DataType type, bool required) const{
    logger::AssertType error = (required) ? logger::ERROR : logger::SILENT;
    bool exists = logger::log_assert(doc.HasMember(name.c_str()), error, "Missing member: {}", name);
    bool correct_type = false;
    switch (type){
        case TYPE_NULL:
            correct_type = logger::log_assert(doc[name.c_str()].IsNull(), error, "Value of key \"{}\" has wrong type, must be null.", name);
            break;
        case TYPE_BOOL:
            correct_type = logger::log_assert(doc[name.c_str()].IsBool(), error, "Value of key \"{}\" has wrong type, must be boolean.", name);
            break;
        case TYPE_INT:
            correct_type = logger::log_assert(doc[name.c_str()].IsInt(), error, "Value of key \"{}\" has wrong type, must be an integer.", name);
            break;
        case TYPE_DOUBLE:
            correct_type = logger::log_assert(doc[name.c_str()].IsNumber(), error, "Value of key \"{}\" has wrong type, must be a number.", name);
            break;
        case TYPE_STRING:
            correct_type = logger::log_assert(doc[name.c_str()].IsString(), error, "Value of key \"{}\" has wrong type, must be a string.", name);
            break;
        case TYPE_ARRAY:
            correct_type = logger::log_assert(doc[name.c_str()].IsArray(), error, "Value of key \"{}\" has wrong type, must be an array.", name);
            break;
        case TYPE_OBJECT:
            correct_type = logger::log_assert(doc[name.c_str()].IsObject(), error, "Value of key \"{}\" has wrong type, must be an object.", name);
            break;
    }
    return exists && correct_type;
}

#endif
