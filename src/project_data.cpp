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
#include "logger.hpp"


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

    }


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
            correct_type = logger::log_assert(doc[name.c_str()].IsDouble(), error, "Value of key \"{}\" has wrong type, must be a number.", name);
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
