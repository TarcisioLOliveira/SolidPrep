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
#include "project_data.hpp"
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/error/error.h"
#include "rapidjson/error/en.h"
#include "STEPCAFControl_Reader.hxx"
#include "BRepClass3d_SolidClassifier.hxx"
#include "logger.hpp"
#define _USE_MATH_DEFINES
#include <cmath>


ProjectData::ProjectData(std::string project_file){
#ifdef _WIN32
    FILE* fp = fopen(project_file.c_str(), "rb");
#else
    FILE* fp = fopen(project_file.c_str(), "r");
#endif

    char readBuffer[65536];
    rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));

    rapidjson::Document doc;
    rapidjson::ParseResult ok = doc.ParseStream<rapidjson::kParseCommentsFlag>(is);
    logger::log_assert(ok, logger::ERROR, "JSON parse error: {} ({}) \n", rapidjson::GetParseError_En(ok.Code()), ok.Offset());
    logger::log_assert(doc.IsObject(), logger::ERROR, "The root of the JSON file must be an object.");

    logger::log_assert(doc.HasMember("solid_type"), logger::ERROR, "Missing member: ");
    if(this->log_data(doc, "solid_type", TYPE_STRING, true)){
        std::string solid_type = doc["solid_type"].GetString();
        if(solid_type == "2D"){
            this->type = TYPE_2D;
        } else if(solid_type == "3D"){
            this->type = TYPE_3D;
        } else {
            logger::log_assert(false, logger::ERROR,  "Solid type incorrectly specified, must be \"2D\" or \"3D\".");
        }
    }
    if(this->log_data(doc, "geometry_path", TYPE_STRING, true)){
        STEPControl_Reader reader;
        IFSelect_ReturnStatus stat = reader.ReadFile("tmp/square.step");
        if(stat != IFSelect_RetDone){
            reader.PrintCheckLoad(false, IFSelect_ItemsByEntity);
            exit(EXIT_FAILURE);
        }

        Standard_Integer NbRoots = reader.NbRootsForTransfer();
        Standard_Integer num = reader.TransferRoots();
        this->solid = reader.OneShape();
        if(this->solid.IsNull()){
            reader.PrintCheckTransfer(true, IFSelect_ItemsByEntity);
            exit(EXIT_FAILURE);
        }
    }
    if(this->log_data(doc, "step", TYPE_DOUBLE, true)){
        this->step = doc["step"].GetDouble();
    }
    if(this->log_data(doc, "restriction_size", TYPE_DOUBLE, true)){
        this->restriction = doc["restriction_size"].GetDouble();
    }
    if(this->log_data(doc, "scale", TYPE_DOUBLE, true)){
        this->scale = doc["scale"].GetDouble();
    }
    if(this->type == TYPE_2D){
        if(this->log_data(doc, "thickness", TYPE_DOUBLE, true)){
            this->thickness = doc["thickness"].GetDouble();
        }
    }
    double angle = 0;
    int choices = 0;
    if(this->log_data(doc, "max_turn_angle", TYPE_DOUBLE, true)){
        angle = doc["max_turn_angle"].GetDouble();
    }
    if(this->log_data(doc, "turn_options", TYPE_INT, true)){
        choices = doc["turn_options"].GetInt();
    }
    if(this->type == TYPE_2D){
        double angle_step = (M_PI/360)*2*angle/(choices-1);
        for(int i = 0; i < choices; ++i){
            this->angles2D.push_back(-angle + i*angle_step);
        }
    } else if(this->type == TYPE_3D){
        double angle_step = (M_PI/360)*2*angle/(choices-1);
        double spin_step = 2*M_PI/(choices-1);
        bool odd = (choices % 2) == 1;
        if(odd){
            std::array<double, 2> arr = {0, 0};
            this->angles3D.push_back(arr);
            for(int j = 0; j < choices - 1; ++j){
                for(int i = 0; i < choices; ++i){
                    arr[0] = -angle + i*angle_step;
                    if(arr[0] != 0){
                        arr[1] = j*spin_step;
                        this->angles3D.push_back(arr);
                    }
                }
            }
        } else {
            std::array<double, 2> arr = {0, 0};
            for(int j = 0; j < choices+1; ++j){
                for(int i = 0; i < choices; ++i){
                    arr[0] = -angle + i*angle_step;
                    arr[1] = j*spin_step;
                    this->angles3D.push_back(arr);
                }
            }
        }
    }

    // Get forces
    // Get supports

    fclose(fp);
}

bool ProjectData::log_data(const rapidjson::Document& doc, std::string name, ProjectData::DataType type, bool required){
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
