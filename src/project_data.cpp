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
        double scale;
        if(this->log_data(doc, "scale", TYPE_DOUBLE, false)){
            scale = doc["scale"].GetDouble();
        } else {
            scale = 1;
        }
        this->ground_structure.reset(new GroundStructure(absolute_path, scale, this->type));
    }

    bool needs_sizing = true;
    bool needs_topopt = true;
    if(this->log_data(doc, "analysis", TYPE_STRING, true)){
        std::string a = doc["analysis"].GetString();
        if(a == "complete"){
            this->analysis = COMPLETE;
        } else if(a == "fea_only"){
            this->analysis = FEA_ONLY;
            needs_sizing = false;
            needs_topopt = false;
        } else if(a == "beams_only"){
            this->analysis = BEAMS_ONLY;
            needs_topopt = false;
        } else if(a == "topopt_only"){
            this->analysis = OPTIMIZE_ONLY;
            needs_sizing = false;
        }
    }
    if(this->type == utils::PROBLEM_TYPE_2D){
        if(this->log_data(doc, "thickness", TYPE_DOUBLE, true)){
            this->thickness = doc["thickness"].GetDouble();
        }
    }
    if(this->log_data(doc, "material", TYPE_OBJECT, true)){
        this->material = this->load_material(doc);
    }
    if(this->log_data(doc, "mesher", TYPE_OBJECT, true)){
        this->log_data(doc["mesher"], "element_type", TYPE_STRING, true);
        this->topopt_element = this->get_element_type(doc["mesher"]["element_type"]);
        this->topopt_mesher = this->load_mesher(doc);
    }
    if(this->log_data(doc, "finite_element", TYPE_OBJECT, false)){
        this->topopt_fea = this->load_fea(doc);
    }
    if(this->log_data(doc, "topopt", TYPE_OBJECT, needs_topopt)){
        this->topopt = this->load_topopt(doc);
    }
    if(this->log_data(doc, "sizing", TYPE_OBJECT, needs_sizing)){
        if(this->log_data(doc["sizing"], "pathfinding", TYPE_OBJECT, false)){
            this->pathfinder = this->load_pathfinder(doc["sizing"]);
        }
        if(this->log_data(doc["sizing"], "finite_element", TYPE_OBJECT, false)){
            this->sizer_fea = this->load_fea(doc["sizing"]);
        }
        this->sizer = this->load_sizer(doc);
    }
    if(this->log_data(doc, "loads", TYPE_ARRAY, true)){
        this->forces = this->get_loads(doc["loads"]);
    }
    if(this->log_data(doc, "supports", TYPE_ARRAY, true)){
        this->supports = this->get_support(doc["supports"]);
    }


    fclose(fp);
}

