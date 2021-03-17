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
#include "pathfinding/meshless_astar.hpp"
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/error/error.h"
#include "rapidjson/error/en.h"
#include "STEPCAFControl_Reader.hxx"
#include "BRepClass3d_SolidClassifier.hxx"
#include "BRepBuilderAPI_Transform.hxx"
#include "logger.hpp"
#define _USE_MATH_DEFINES
#include <cmath>
#include "force.hpp"


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
    if(this->log_data(doc, "scale", TYPE_DOUBLE, false)){
        this->scale = doc["scale"].GetDouble();
        gp_Trsf t;
        t.SetScaleFactor(this->scale);
        BRepBuilderAPI_Transform transf(this->solid, t);
    } else {
        this->scale = 1;
    }
    if(this->type == TYPE_2D){
        if(this->log_data(doc, "thickness", TYPE_DOUBLE, true)){
            this->thickness = doc["thickness"].GetDouble();
        }
    }
    if(this->log_data(doc, "pathfinding", TYPE_OBJECT, true)){
        using namespace pathfinding;

        auto& pathf = doc["pathfinding"];
        this->log_data(pathf, "type", TYPE_STRING, true);
        if(pathf["type"] == "meshless_astar"){
            this->log_data(pathf, "step", TYPE_DOUBLE, true);
            this->log_data(pathf, "max_turn_angle", TYPE_DOUBLE, true);
            this->log_data(pathf, "turn_options", TYPE_INT, true);
            double step = pathf["step"].GetDouble();
            double angle = pathf["max_turn_angle"].GetDouble();
            int choices = pathf["turn_options"].GetInt();
            double restriction = 0;
            if(this->log_data(pathf, "restriction_size", TYPE_DOUBLE, false)){
                restriction = pathf["restriction_size"].GetDouble();
            }
            this->pathfinder = new MeshlessAStar(this->solid, step, angle, choices, restriction, MeshlessAStar::TYPE_2D);
        }
    }
    if(this->log_data(doc, "loads", TYPE_ARRAY, true)){
        if(this->type == TYPE_2D){
            for(auto& f : doc["loads"].GetArray()){
                logger::log_assert(f.IsObject(), logger::ERROR, "Each load must be stored as a JSON object");
                this->log_data(f, "vertices", TYPE_ARRAY, true);
                this->log_data(f, "load", TYPE_ARRAY, true);
                auto vertices = f["vertices"].GetArray();
                auto loads = f["load"].GetArray();
                logger::log_assert(loads.Size() == 2, logger::ERROR, "Load vector must have exactly two dimensions in 2D problems");

                std::array<double, 2> l{loads[0].GetDouble(), loads[1].GetDouble()};
                std::vector<std::array<double, 2>> vlist;
                for(auto& v : vertices){
                    logger::log_assert(v.Size() == 2, logger::ERROR, "Vertices must have exactly two dimensions in 2D problems");
                    vlist.push_back({v[0].GetDouble(), v[1].GetDouble()});
                }
                this->forces.emplace_back(vlist, this->thickness, l);
            }
        } else if(this->type == TYPE_3D) {
            for(auto& f : doc["loads"].GetArray()){
                logger::log_assert(f.IsObject(), logger::ERROR, "Each load must be stored as a JSON object");
                this->log_data(f, "vertices", TYPE_ARRAY, true);
                this->log_data(f, "load", TYPE_ARRAY, true);
                auto vertices = f["vertices"].GetArray();
                auto loads = f["load"].GetArray();
                logger::log_assert(loads.Size() == 2, logger::ERROR, "Load vector must have exactly three dimensions in 3D problems");

                std::array<double, 3> l{loads[0].GetDouble(), loads[1].GetDouble(), loads[2].GetDouble()};
                std::vector<std::array<double, 3>> vlist;
                for(auto& v : vertices){
                    logger::log_assert(v.Size() == 2, logger::ERROR, "Vertices must have exactly three dimensions in 3D problems");
                    vlist.push_back({v[0].GetDouble(), v[1].GetDouble(), v[2].GetDouble()});
                }
                this->forces.emplace_back(vlist, l);
            }
        }
    }
    if(this->log_data(doc, "supports", TYPE_ARRAY, true)){
        if(this->type == TYPE_2D){
            for(auto& f : doc["supports"].GetArray()){
                logger::log_assert(f.IsObject(), logger::ERROR, "Each support must be stored as a JSON object");
                this->log_data(f, "vertices", TYPE_ARRAY, true);
                this->log_data(f, "X", TYPE_BOOL, true);
                this->log_data(f, "Y", TYPE_BOOL, true);
                bool X = f["X"].GetBool();
                bool Y = f["Y"].GetBool();

                auto vertices = f["vertices"].GetArray();
                std::vector<std::array<double, 2>> vlist;
                for(auto& v : vertices){
                    logger::log_assert(v.Size() == 2, logger::ERROR, "Vertices must have exactly two dimensions in 2D problems");
                    vlist.push_back({v[0].GetDouble(), v[1].GetDouble()});
                }
                this->supports.emplace_back(X, Y, vlist);
            }
        } else if(this->type == TYPE_3D) {
            for(auto& f : doc["supports"].GetArray()){
                logger::log_assert(f.IsObject(), logger::ERROR, "Each support must be stored as a JSON object");
                this->log_data(f, "vertices", TYPE_ARRAY, true);
                this->log_data(f, "X", TYPE_BOOL, true);
                this->log_data(f, "Y", TYPE_BOOL, true);
                this->log_data(f, "Z", TYPE_BOOL, true);
                bool X = f["X"].GetBool();
                bool Y = f["Y"].GetBool();
                bool Z = f["Z"].GetBool();

                auto vertices = f["vertices"].GetArray();
                std::vector<std::array<double, 3>> vlist;
                for(auto& v : vertices){
                    logger::log_assert(v.Size() == 2, logger::ERROR, "Vertices must have exactly three dimensions in 3D problems");
                    vlist.push_back({v[0].GetDouble(), v[1].GetDouble(), v[2].GetDouble()});
                }
                this->supports.emplace_back(X, Y, Z, vlist);
            }
        }
    }


    fclose(fp);
}

ProjectData::~ProjectData(){
    delete this->pathfinder;
}
