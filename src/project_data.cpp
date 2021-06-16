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
#include "pathfinding/visibility_graph.hpp"
#include "material/linear_elastic_isotropic.hpp"
#include "material/linear_elastic_orthotropic.hpp"
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/error/error.h"
#include "rapidjson/error/en.h"
#include "logger.hpp"
#include "rapidjson/rapidjson.h"
#include "sizing/beam_sizing.hpp"
#include "utils.hpp"
#define _USE_MATH_DEFINES
#include <cmath>
#include "force.hpp"
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
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
            this->type = utils::PROBLEM_TYPE_2D;
        } else if(solid_type == "3D"){
            this->type = utils::PROBLEM_TYPE_3D;
        } else {
            logger::log_assert(false, logger::ERROR,  "Solid type incorrectly specified, must be \"2D\" or \"3D\".");
        }
    }
    if(this->log_data(doc, "geometry_path", TYPE_STRING, true)){
        std::string geom_path = doc["geometry_path"].GetString();
#ifdef _WIN32
        size_t last_slash = project_file.rfind("\\");
#else
        size_t last_slash = project_file.rfind("/");
#endif
        std::string absolute_path = project_file.substr(0, last_slash+1);
        absolute_path.append(geom_path);
        float scale;
        if(this->log_data(doc, "scale", TYPE_DOUBLE, false)){
            scale = doc["scale"].GetDouble();
        } else {
            scale = 1;
        }
        this->ground_structure.reset(new GroundStructure(absolute_path, scale, this->type));
    }
    if(this->type == utils::PROBLEM_TYPE_2D){
        if(this->log_data(doc, "thickness", TYPE_DOUBLE, true)){
            this->thickness = doc["thickness"].GetDouble()*1e-3;
        }
    }
    if(this->log_data(doc, "material", TYPE_OBJECT, true)){
        auto& mat = doc["material"];
        this->log_data(mat, "type", TYPE_STRING, true);
        if(mat["type"] == "linear_elastic_orthotropic"){
            std::vector<std::string> properties{"E", "nu", "G", "Smax", "Tmax"};
            for(auto& s:properties){
                logger::log_assert(mat.HasMember(s.c_str()), logger::ERROR, "missing material property: {}", s);
                logger::log_assert(mat[s.c_str()].IsArray() || mat[s.c_str()].IsDouble(), logger::ERROR, "material property {} must be either a number or an array of numbers", s);
            }
            std::vector<std::vector<float>> values(5);
            for(size_t i = 0; i < properties.size(); ++i){
                if(mat[properties[i].c_str()].IsArray()){
                    const auto& a = mat[properties[i].c_str()].GetArray();
                    if(a.Size() == 1){
                        values[i].resize(3, a[0].GetDouble());
                    } else {
                        values[i].resize(3, 0);
                        for(size_t j = 0; j < std::min(a.Size(), (rapidjson::SizeType) 3); ++j){
                            values[i][j] = a[j].GetDouble();
                        }
                    }
                } else {
                    values[i].resize(3, mat[properties[i].c_str()].GetDouble());
                }
            }
            for(auto& i:values[0]) i *= 1e9; // E
            for(auto& i:values[2]) i *= 1e9; // G
            for(auto& i:values[3]) i *= 1e6; // Smax
            for(auto& i:values[4]) i *= 1e6; // Tmax
            this->material.reset(new material::LinearElasticOrthotropic(values[0], values[1], values[2], values[3], values[4]));
        } else if(mat["type"] == "linear_elastic_isotropic"){
            std::vector<std::string> properties{"E", "nu", "Smax", "Tmax"};
            for(auto& s:properties){
                this->log_data(mat, s, TYPE_DOUBLE, true);
            }
            this->log_data(mat, "plane_stress", TYPE_BOOL, true);
            float E = mat["E"].GetDouble();
            float nu = mat["nu"].GetDouble();
            float Smax = mat["Smax"].GetDouble();
            float Tmax = mat["Tmax"].GetDouble();
            bool plane_stress = mat["plane_stress"].GetBool();
            this->material.reset(new material::LinearElasticIsotropic(E*1e9, nu, Smax*1e6, Tmax*1e6, plane_stress));
        }
    }
    if(this->log_data(doc, "sizing", TYPE_OBJECT, true)){
        auto& sizing = doc["sizing"];
        this->log_data(sizing, "type", TYPE_STRING, true);
        if(sizing["type"] == "beam_sizing"){
            if(this->log_data(sizing, "pathfinding", TYPE_OBJECT, true)){
                using namespace pathfinding;

                auto& pathf = sizing["pathfinding"];
                this->log_data(pathf, "type", TYPE_STRING, true);
                if(pathf["type"] == "meshless_astar"){
                    this->log_data(pathf, "step", TYPE_DOUBLE, true);
                    this->log_data(pathf, "max_turn_angle", TYPE_DOUBLE, true);
                    this->log_data(pathf, "turn_options", TYPE_INT, true);
                    float step = pathf["step"].GetDouble();
                    float angle = pathf["max_turn_angle"].GetDouble();
                    int choices = pathf["turn_options"].GetInt();
                    float restriction = 0;
                    if(this->log_data(pathf, "restriction_size", TYPE_DOUBLE, false)){
                        restriction = pathf["restriction_size"].GetDouble();
                    }
                    this->pathfinder.reset(new MeshlessAStar(this->ground_structure->shape, step, angle, choices, restriction, utils::PROBLEM_TYPE_2D));
                } else if(pathf["type"] == "visibility_graph"){
                    this->log_data(pathf, "step", TYPE_DOUBLE, true);
                    this->log_data(pathf, "max_turn_angle", TYPE_DOUBLE, true);
                    float step = pathf["step"].GetDouble();
                    float angle = pathf["max_turn_angle"].GetDouble();
                    float restriction = 0;
                    if(this->log_data(pathf, "restriction_size", TYPE_DOUBLE, false)){
                        restriction = pathf["restriction_size"].GetDouble();
                    }
                    this->pathfinder.reset(new VisibilityGraph(this->ground_structure.get(), step, angle, restriction, utils::PROBLEM_TYPE_2D));
                } else {
                    logger::log_assert(false, logger::ERROR, "unknown pathfinding algorithm inserted: {}.", pathf["type"].GetString());
                }
            }
            this->log_data(sizing, "element_type", TYPE_STRING, true);
            BeamElementFactory::BeamElementType t = BeamElementFactory::NONE;
            if(sizing["element_type"] == "beam_linear_2D"){
                t = BeamElementFactory::BEAM_LINEAR_2D;
            } else {
                logger::log_assert(false, logger::ERROR, "unknown element type for sizing algorithm: {}.", sizing["element_type"].GetString());
            }
            this->sizer.reset(new sizing::BeamSizing(this, t));
        }
    }
    if(this->log_data(doc, "loads", TYPE_ARRAY, true)){
        if(this->type == utils::PROBLEM_TYPE_2D){
            for(auto& f : doc["loads"].GetArray()){
                logger::log_assert(f.IsObject(), logger::ERROR, "Each load must be stored as a JSON object");
                this->log_data(f, "vertices", TYPE_ARRAY, true);
                this->log_data(f, "load", TYPE_ARRAY, true);
                auto vertices = f["vertices"].GetArray();
                auto loads = f["load"].GetArray();
                logger::log_assert(loads.Size() == 2, logger::ERROR, "Load vector must have exactly two dimensions in 2D problems");

                gp_Vec l(loads[0].GetDouble(), loads[1].GetDouble(), 0);
                std::vector<gp_Pnt> vlist;
                for(auto& v : vertices){
                    logger::log_assert(v.Size() == 2, logger::ERROR, "Vertices must have exactly two dimensions in 2D problems");
                    vlist.emplace_back(v[0].GetDouble(), v[1].GetDouble(), 0);
                }
                CrossSection S(vlist, this->thickness);
                this->forces.emplace_back(S, l);
            }
        } else if(this->type == utils::PROBLEM_TYPE_3D) {
            for(auto& f : doc["loads"].GetArray()){
                logger::log_assert(f.IsObject(), logger::ERROR, "Each load must be stored as a JSON object");
                this->log_data(f, "vertices", TYPE_ARRAY, true);
                this->log_data(f, "load", TYPE_ARRAY, true);
                auto vertices = f["vertices"].GetArray();
                auto loads = f["load"].GetArray();
                logger::log_assert(loads.Size() == 2, logger::ERROR, "Load vector must have exactly three dimensions in 3D problems");

                gp_Vec l(loads[0].GetDouble(), loads[1].GetDouble(), loads[2].GetDouble());
                std::vector<gp_Pnt> vlist;
                for(auto& v : vertices){
                    logger::log_assert(v.Size() == 2, logger::ERROR, "Vertices must have exactly three dimensions in 3D problems");
                    vlist.emplace_back(v[0].GetDouble(), v[1].GetDouble(), v[2].GetDouble());
                }
                CrossSection S(vlist);
                this->forces.emplace_back(S, l);
            }
        }
    }
    if(this->log_data(doc, "supports", TYPE_ARRAY, true)){
        if(this->type == utils::PROBLEM_TYPE_2D){
            for(auto& f : doc["supports"].GetArray()){
                logger::log_assert(f.IsObject(), logger::ERROR, "Each support must be stored as a JSON object");
                this->log_data(f, "vertices", TYPE_ARRAY, true);
                this->log_data(f, "X", TYPE_BOOL, true);
                this->log_data(f, "Y", TYPE_BOOL, true);
                this->log_data(f, "MZ", TYPE_BOOL, true);
                bool X = f["X"].GetBool();
                bool Y = f["Y"].GetBool();
                bool MZ = f["MZ"].GetBool();

                auto vertices = f["vertices"].GetArray();
                std::vector<gp_Pnt> vlist;
                for(auto& v : vertices){
                    logger::log_assert(v.Size() == 2, logger::ERROR, "Vertices must have exactly two dimensions in 2D problems");
                    vlist.emplace_back(v[0].GetDouble(), v[1].GetDouble(), 0);
                }
                CrossSection S(vlist, this->thickness);
                this->supports.emplace_back(X, Y, MZ, S);
            }
        } else if(this->type == utils::PROBLEM_TYPE_3D) {
            for(auto& f : doc["supports"].GetArray()){
                logger::log_assert(f.IsObject(), logger::ERROR, "Each support must be stored as a JSON object");
                this->log_data(f, "vertices", TYPE_ARRAY, true);
                this->log_data(f, "X", TYPE_BOOL, true);
                this->log_data(f, "Y", TYPE_BOOL, true);
                this->log_data(f, "Z", TYPE_BOOL, true);
                this->log_data(f, "MX", TYPE_BOOL, true);
                this->log_data(f, "MY", TYPE_BOOL, true);
                this->log_data(f, "MZ", TYPE_BOOL, true);
                bool X = f["X"].GetBool();
                bool Y = f["Y"].GetBool();
                bool Z = f["Z"].GetBool();
                bool MX = f["MX"].GetBool();
                bool MY = f["MY"].GetBool();
                bool MZ = f["MZ"].GetBool();

                auto vertices = f["vertices"].GetArray();
                std::vector<gp_Pnt> vlist;
                for(auto& v : vertices){
                    logger::log_assert(v.Size() == 2, logger::ERROR, "Vertices must have exactly three dimensions in 3D problems");
                    vlist.emplace_back(v[0].GetDouble(), v[1].GetDouble(), v[2].GetDouble());
                }
                CrossSection S(vlist);
                this->supports.emplace_back(X, Y, Z, MX, MY, MZ, S);
            }
        }
    }


    fclose(fp);
}
