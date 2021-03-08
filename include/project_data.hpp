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
    std::vector<double> angles;
    double restriction;
    double scale;
    double thickness;
    ProjectType type;
    TopoDS_Shape solid;
    // std::vector<Force> forces;
    // std::vector<Support> supports;
    
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
    bool log_data(const rapidjson::Document& doc, std::string name, DataType type, bool required);
};


#endif
