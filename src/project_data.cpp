/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
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
#include "meshing/mesh_file.hpp"
#include "project_specification/data_map.hpp"
#include "project_specification/registry.hpp"
#include "shape_handler.hpp"
#include "utils/delayed_pointer.hpp"
#include <cstring>
#include <json/json.h>
#include <memory>
#include <vector>
#define _USE_MATH_DEFINES
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include "utils.hpp"
#include "force.hpp"
#include "logger.hpp"
#include "project_data.hpp"
#include "meshing/gmsh.hpp"
#include "element/GT9.hpp"
#include "element/TRI3.hpp"
#include "element/Q4.hpp"
#include "element/Q4S.hpp"
#include "element/TET4.hpp"
#include "element/H8.hpp"
#include "element/TET10.hpp"
#include "projection/none.hpp"

ProjectData::ProjectData(std::string project_file){
    std::ifstream f(project_file);
    Json::Value settings;
    f >> settings;

    this->folder_path = this->get_folder_path(project_file);
    if(this->log_data(settings, "contact_type", projspec::TYPE_OBJECT, false)){
        this->contact_data = this->get_contact_data(settings["contact_type"]);
    } else {
        this->contact_data = FiniteElement::ContactData{FiniteElement::ContactType::RIGID, 0.0};
    }

    logger::log_assert(settings.isMember("solid_type"), logger::ERROR, "Missing member: ");
    if(this->log_data(settings, "solid_type", projspec::TYPE_STRING, true)){
        std::string solid_type = settings["solid_type"].asString();
        if(solid_type == "2D"){
            this->type = utils::PROBLEM_TYPE_2D;
        } else if(solid_type == "3D"){
            this->type = utils::PROBLEM_TYPE_3D;
        } else {
            logger::log_assert(false, logger::ERROR,  "Solid type incorrectly specified, must be \"2D\" or \"3D\".");
        }
    }

    if(this->log_data(settings, "analysis", projspec::TYPE_OBJECT, true)){
        const auto& analysis = settings["analysis"];
        if(this->log_data(analysis, "simulation", projspec::TYPE_BOOL, false)){
            this->do_simulation = analysis["simulation"].asBool();
        }
        if(!this->do_simulation){
            if(this->log_data(analysis, "gen", projspec::TYPE_BOOL, false)){
                this->generate_beams = analysis["gen"].asBool();
            }
            if(this->log_data(analysis, "meshing", projspec::TYPE_BOOL, false)){
                this->do_meshing = analysis["meshing"].asBool();
            } else if(this->generate_beams){
                this->do_meshing = true;
            }
            if(this->generate_beams){
                logger::log_assert(this->do_meshing, logger::WARNING, "mesh loading not supported for geometry generation, ignoring setting");
                this->do_meshing = true;
            }
            if(this->log_data(analysis, "topopt", projspec::TYPE_BOOL, false)){
                this->do_topopt = analysis["topopt"].asBool();
            }
            if(this->log_data(analysis, "fea", projspec::TYPE_BOOL, false)){
                this->do_fea = analysis["fea"].asBool();
            }
            if(this->log_data(analysis, "shape_opt", projspec::TYPE_BOOL, false)){
                this->do_shape_opt = analysis["shape_opt"].asBool();
            }
            if(this->do_shape_opt){
                logger::log_assert(!this->do_topopt, logger::ERROR, "topology optimization followed by shape optimization is currently not supported");
                logger::log_assert(!this->generate_beams, logger::ERROR, "geometry generation followed by shape optimization is currently not supported");
            }
        }
        logger::log_assert(this->do_meshing || 
                           this->do_fea ||
                           this->generate_beams || 
                           this->do_topopt || 
                           this->do_shape_opt || 
                           this->do_simulation, 
                    logger::ERROR, "No analysis type was set");
    }
    if(this->generate_beams){
        logger::log_assert(this->type == utils::PROBLEM_TYPE_2D,
                           logger::ERROR, "sizing is currently only implemented for 2D elasticity problems.");
    }
    if(this->type == utils::PROBLEM_TYPE_2D){
        if(this->log_data(settings, "thickness", projspec::TYPE_DOUBLE, true)){
            this->thickness = settings["thickness"].asDouble();
        }
    }

    projspec::RequirementConditions conds{
        .do_meshing = this->do_meshing,
        .generate_beams = this->generate_beams,
        .do_topopt = this->do_topopt,
        .do_fea = this->do_fea,
        .do_shape_opt = this->do_shape_opt,
        .do_simulation = this->do_simulation,
        .problem_type = this->type
    };


    if(this->log_data(settings, "loads", projspec::TYPE_ARRAY, false)){
        this->forces = this->get_loads(settings["loads"]);
    }
    if(this->log_data(settings, "supports", projspec::TYPE_ARRAY, false)){
        this->supports = this->get_support(settings["supports"]);
    }
    if(this->log_data(settings, "mesher", projspec::TYPE_OBJECT, true)){
        this->log_data(settings["mesher"], "element_type", projspec::TYPE_STRING, true);
        if(this->do_meshing){
            this->element_name = settings["mesher"]["element_type"].asString();
        } else {
            std::string absolute_path = this->folder_path;
            std::string mesh_path = settings["mesher"]["file_path"].asString();
            absolute_path.append(mesh_path);

            std::ifstream file(absolute_path, std::ios::in);
            logger::log_assert(file.good(), logger::ERROR, "file not found: {}", absolute_path);
            std::getline(file, this->element_name);
            file.close();
        }
        if(this->element_name == "Q4S"){
            logger::log_assert(!this->do_shape_opt, logger::ERROR, "Q4S not supported for shape optimization as it assumes element shapes are rectangular");
        }
        this->topopt_element = projspec::ElementRegistry::get(this->element_name);
        logger::log_assert(this->topopt_element != nullptr, logger::ERROR, "element type not found: {}", this->element_name);
    }
    if(this->log_data(settings, Material::get_name(), projspec::TYPE_ARRAY, true)){
        this->materials = this->generate_material_stubs(settings, conds);
    }
    if(this->log_data(settings, Field::get_name(), projspec::TYPE_ARRAY, false)){
        this->fields = this->generate_field_stubs(settings, conds);
    }
    if(this->log_data(settings, "geometry", projspec::TYPE_ARRAY, true)){
        this->geometries = this->load_geometries(settings);
    }

    this->load_fields(settings, this->fields, conds);
    this->load_materials(settings, this->materials, conds);

    if(this->log_data(settings, "springs", projspec::TYPE_ARRAY, false)){
        this->springs = this->get_springs(settings["springs"]);
    }
    logger::log_assert(this->supports.size() > 0 || this->springs.size() > 0, logger::ERROR,
                       "a support or spring must be specified to prevent matrix singularity.");
    if(this->log_data(settings, "internal_loads", projspec::TYPE_ARRAY, false)){
        this->internal_loads = this->get_internal_loads(settings["internal_loads"]);
    }
    if(this->log_data(settings, "sub_problems", projspec::TYPE_ARRAY, false)){
        this->sub_problems = this->load_sub_problems(settings);
    } else {
        this->sub_problems.emplace_back();
        for(auto& i:this->supports){
            this->sub_problems[0].supports.push_back(&i);
        }
        for(auto& i:this->forces){
            this->sub_problems[0].forces.push_back(&i);
        }
        for(auto& i:this->springs){
            this->sub_problems[0].springs.push_back(&i);
        }
        for(auto& i:this->internal_loads){
            this->sub_problems[0].internal_loads.push_back(&i);
        }
    }
    if(this->log_data(settings, FiniteElement::get_name(), projspec::TYPE_OBJECT, true)){
        std::vector<std::unique_ptr<FiniteElement>> solvers(this->sub_problems.size());
        // Dirty, but quicker to implement and I don't have to mess with
        // `sizer_fea`down there.
        for(size_t i = 0; i < this->sub_problems.size(); ++i){
            solvers[i] = this->load_fea(settings, conds);
        }
        this->topopt_fea = std::make_unique<SolverManager>(std::move(solvers));
    }
    if(this->log_data(settings, "sizing", projspec::TYPE_OBJECT, this->generate_beams) && this->generate_beams){
        if(this->log_data(settings["sizing"], "pathfinding", projspec::TYPE_OBJECT, false)){
            this->pathfinder = this->load_pathfinder(settings["sizing"], conds);
        }
        if(this->log_data(settings["sizing"], "finite_element", projspec::TYPE_OBJECT, false)){
            this->sizer_fea = this->load_fea(settings["sizing"], conds);
        }
        this->sizer = this->load_sizer(settings, conds);
    }
    if(this->log_data(settings, "mesher", projspec::TYPE_OBJECT, true)){
        this->topopt_mesher = this->load_mesher(settings, conds);
    }
    if(this->log_data(settings, "topopt", projspec::TYPE_OBJECT, this->do_topopt) && this->do_topopt){
        this->topopt_optimizer = this->load_topopt_optimizer(settings, conds);
    }
    if(this->log_data(settings, "shape_opt", projspec::TYPE_OBJECT, this->do_shape_opt) && this->do_shape_opt){
        this->shopt_optimizer = this->load_shopt_optimizer(settings, conds);
    }
    if(this->log_data(settings, "simulation", projspec::TYPE_OBJECT, this->do_simulation) && this->do_simulation){
        this->simulator = this->load_simulation(settings, conds);
    }
}


std::string ProjectData::get_folder_path(const std::string& project_file_path) const{
#ifdef _WIN32
    size_t last_slash = project_file_path.rfind("\\");
#else
    size_t last_slash = project_file_path.rfind("/");
#endif

    return project_file_path.substr(0, last_slash+1);
}

utils::DelayedPointerView<Material> ProjectData::get_material(const std::string& name){
    if(name == ""){
        return nullptr;
    }
    for(auto& m:this->materials){
        if(m->material_name == name){
            return m.get_view();
        }
    }
    logger::log_assert(false, logger::ERROR, "material with name {} not found", name);

    return nullptr;
}

inline void add_stub(const projspec::DataEntry& r, projspec::DataMap* data);

inline void add_stub_array(ProjectData* proj, const std::shared_ptr<projspec::ArrayRequirements>& reqs, projspec::DataArray* array){
    const size_t size = std::max((size_t)1, reqs->size);
    switch(reqs->type){
        case projspec::TYPE_POINTER:
            break;
        case projspec::TYPE_NULL:
            break;
        case projspec::TYPE_BOOL:
            for(size_t i = 0; i < size; ++i){
                array->set_bool(i, false);
            }
            break;
        case projspec::TYPE_INT:
            for(size_t i = 0; i < size; ++i){
                array->set_int(i, 0);
            }
            break;
        case projspec::TYPE_DOUBLE:
            for(size_t i = 0; i < size; ++i){
                array->set_double(i, 0);
            }
            break;
        case projspec::TYPE_STRING:
        case projspec::TYPE_RELATIVE_PATH:
            for(size_t i = 0; i < size; ++i){
                array->set_string(i, "");
            }
            break;
        case projspec::TYPE_ARRAY:{
            for(size_t i = 0; i < size; ++i){
                std::unique_ptr<projspec::DataArray> subarray = std::make_unique<projspec::DataArray>();
                add_stub_array(proj, reqs->array_data, subarray.get());
                array->set_array(i, std::move(subarray));
            }
            break;
        }
        case projspec::TYPE_OBJECT:{
            for(size_t i = 0; i < size; ++i){
                std::unique_ptr<projspec::DataMap> subdata = std::make_unique<projspec::DataMap>(proj);
                for(auto& s:reqs->object_data){
                    add_stub(s, subdata.get());
                }
                array->set_object(i, std::move(subdata));
            }
            break;
        }
        case projspec::TYPE_CROSS_SECTION:
            for(size_t i = 0; i < size; ++i){
                array->set_cross_section(i, CrossSection());
            }
            break;
    }
}
inline void add_stub(const projspec::DataEntry& r, projspec::DataMap* data){
    switch(r.type){
        case projspec::TYPE_POINTER:
            break;
        case projspec::TYPE_NULL:
            break;
        case projspec::TYPE_BOOL:
            data->set_bool(r.name, false);
            break;
        case projspec::TYPE_INT:
            data->set_int(r.name, 0);
            break;
        case projspec::TYPE_DOUBLE:
            data->set_double(r.name, 0);
            break;
        case projspec::TYPE_STRING:
        case projspec::TYPE_RELATIVE_PATH:
            data->set_string(r.name, "");
            break;
        case projspec::TYPE_ARRAY:{
            std::unique_ptr<projspec::DataArray> array = std::make_unique<projspec::DataArray>();
            add_stub_array(data->proj, r.array_data, array.get());
            data->set_array(r.name, std::move(array));
            break;
        }
        case projspec::TYPE_OBJECT:{
            std::unique_ptr<projspec::DataMap> subdata = std::make_unique<projspec::DataMap>(data->proj);
            for(auto& s:r.object_data){
                add_stub(s, subdata.get());
            }
            data->set_object(r.name, std::move(subdata));
            break;
        }
        case projspec::TYPE_CROSS_SECTION:
            data->set_cross_section(r.name, CrossSection());
            break;
    }
}

projspec::DataMap ProjectData::generate_stub(const projspec::ObjectRequirements& reqs,
                                             const projspec::RequirementConditions& conds){
    projspec::DataMap data(this);

    for(auto& r:reqs.object_entries){
        if(r.required_if && !r.required_if(conds)){
            continue;
        } else if(!r.required) {
            continue;
        }
        add_stub(r, &data);
    }

    return data;
};

inline void add_data(const Json::Value& item, const projspec::DataEntry& r, projspec::DataMap* data);
inline void add_data_array(ProjectData* proj, const Json::Value& item, const std::shared_ptr<projspec::ArrayRequirements>& reqs, projspec::DataArray* array){
    size_t i = 0;
    switch(reqs->type){
        case projspec::TYPE_POINTER:
            break;
        case projspec::TYPE_NULL:
            break;
        case projspec::TYPE_BOOL:
            for(auto& it:item){
                array->set_bool(i, it.asBool());
                ++i;
            }
            break;
        case projspec::TYPE_INT:
            for(auto& it:item){
                array->set_int(i, it.asInt64());
                ++i;
            }
            break;
        case projspec::TYPE_DOUBLE:
            for(auto& it:item){
                array->set_double(i, it.asDouble());
                ++i;
            }
            break;
        case projspec::TYPE_STRING:
            for(auto& it:item){
                array->set_string(i, it.asString());
                ++i;
            }
            break;
        case projspec::TYPE_RELATIVE_PATH:
            for(auto& it:item){
                std::string absolute_path = proj->folder_path;
                std::string rel_path = it.asString();
                absolute_path.append(rel_path);
                array->set_string(i, absolute_path);
                ++i;
            }
            break;
        case projspec::TYPE_ARRAY:{
            for(auto& it:item){
                std::unique_ptr<projspec::DataArray> subarray = std::make_unique<projspec::DataArray>();
                add_data_array(proj, it, reqs->array_data, subarray.get());
                array->set_array(i, std::move(subarray));
                ++i;
            }
            break;
        }
        case projspec::TYPE_OBJECT:{
            for(auto& it:item){
                std::unique_ptr<projspec::DataMap> subdata = std::make_unique<projspec::DataMap>(proj);
                for(auto& s:reqs->object_data){
                    add_data(it[s.name], s, subdata.get());
                }
                array->set_object(i, std::move(subdata));
                ++i;
            }
            break;
        }
        case projspec::TYPE_CROSS_SECTION:
            for(auto& it:item){
                array->set_cross_section(i, proj->get_cross_section(it));
                ++i;
            }
            break;
    }
}
inline void add_data(const Json::Value& item, const projspec::DataEntry& r, projspec::DataMap* data){
    switch(r.type){
        case projspec::TYPE_POINTER:
            break;
        case projspec::TYPE_NULL:
            break;
        case projspec::TYPE_BOOL:
            data->set_bool(r.name, item.asBool());
            break;
        case projspec::TYPE_INT:
            data->set_int(r.name, item.asInt64());
            break;
        case projspec::TYPE_DOUBLE:
            data->set_double(r.name, item.asDouble());
            break;
        case projspec::TYPE_STRING:
            data->set_string(r.name, item.asString());
            break;
        case projspec::TYPE_RELATIVE_PATH:{
            std::string absolute_path = data->proj->folder_path;
            std::string rel_path = item.asString();
            absolute_path.append(rel_path);
            data->set_string(r.name, absolute_path);
            break;
        }
        case projspec::TYPE_ARRAY:{
            std::unique_ptr<projspec::DataArray> array = std::make_unique<projspec::DataArray>();
            if(r.array_data->size > 0){
                logger::log_assert(item.size() == r.array_data->size, logger::ERROR, "array {} has incorrect size: expected {}, got {}", r.name, r.array_data->size, item.size());
            }
            add_data_array(data->proj, item, r.array_data, array.get());
            data->set_array(r.name, std::move(array));
            break;
        }
        case projspec::TYPE_OBJECT:{
            std::unique_ptr<projspec::DataMap> subdata = std::make_unique<projspec::DataMap>(data->proj);
            for(auto& s:r.object_data){
                add_data(item[s.name], s, subdata.get());
            }
            data->set_object(r.name, std::move(subdata));
            break;
        }
        case projspec::TYPE_CROSS_SECTION:
            data->set_cross_section(r.name, data->proj->get_cross_section(item));
            break;
    }
}

projspec::DataMap ProjectData::generate_data_map(const Json::Value& doc,
                                                 const projspec::ObjectRequirements& reqs,
                                                 const projspec::RequirementConditions& conds){
    projspec::DataMap data(this);

    for(auto& r:reqs.object_entries){
        bool required = false;
        if(r.required_if){
            required = r.required_if(conds);
        } else {
            required = r.required;
        }
        if(this->log_data(doc, r.name, r.type, required)){
            add_data(doc[r.name], r, &data);
        }
    }

    return data;
};

bool ProjectData::log_data(const Json::Value& doc, std::string name, projspec::DataType type, bool required) const{
    logger::AssertType error = (required) ? logger::ERROR : logger::SILENT;
    bool exists = logger::log_assert(doc.isMember(name), error, "Missing member: {}", name);
    if(!exists){
        return false;
    }
    bool correct_type = false;
    switch (type){
        case projspec::TYPE_POINTER:
        case projspec::TYPE_NULL:
            correct_type = logger::log_assert(doc[name.c_str()].isNull(), error, "Value of key \"{}\" has wrong type, must be null.", name);
            break;
        case projspec::TYPE_BOOL:
            correct_type = logger::log_assert(doc[name.c_str()].isBool(), error, "Value of key \"{}\" has wrong type, must be boolean.", name);
            break;
        case projspec::TYPE_INT:
            correct_type = logger::log_assert(doc[name.c_str()].isInt(), error, "Value of key \"{}\" has wrong type, must be an integer.", name);
            break;
        case projspec::TYPE_DOUBLE:
            correct_type = logger::log_assert(doc[name.c_str()].isNumeric(), error, "Value of key \"{}\" has wrong type, must be a number.", name);
            break;
        case projspec::TYPE_STRING:
        case projspec::TYPE_RELATIVE_PATH:
            correct_type = logger::log_assert(doc[name.c_str()].isString(), error, "Value of key \"{}\" has wrong type, must be a string.", name);
            break;
        case projspec::TYPE_ARRAY:
            correct_type = logger::log_assert(doc[name.c_str()].isArray(), error, "Value of key \"{}\" has wrong type, must be an array.", name);
            break;
        case projspec::TYPE_CROSS_SECTION:
        case projspec::TYPE_OBJECT:
            correct_type = logger::log_assert(doc[name.c_str()].isObject(), error, "Value of key \"{}\" has wrong type, must be an object.", name);
            break;
    }
    return correct_type;
}

FiniteElement::ContactData ProjectData::get_contact_data(const Json::Value& doc) const{
    if(this->log_data(doc, "type", projspec::TYPE_STRING, true)){
        const std::string type = doc["type"].asString();
        FiniteElement::ContactType contact_type = FiniteElement::ContactType::RIGID;
        double rtol_abs = 0;
        double max_step = 1;
        double step_tol = 0;
        double EPS_DISPL = 0;
        if(type == "rigid"){
            contact_type = FiniteElement::ContactType::RIGID;
        } else {
            this->log_data(doc, "rtol_abs", projspec::TYPE_DOUBLE, true);
            rtol_abs = doc["rtol_abs"].asDouble();
            if(type == "frictionless_penalty"){
                contact_type = FiniteElement::ContactType::FRICTIONLESS_PENALTY;
            } else if(type == "frictionless_displ_log"){
                contact_type = FiniteElement::ContactType::FRICTIONLESS_DISPL_LOG;
            } else if(type == "frictionless_displ_simple"){
                contact_type = FiniteElement::ContactType::FRICTIONLESS_DISPL_SIMPLE;
            //} else if(type == "frictionless_displ_constr"){
            //    contact_type = FiniteElement::ContactType::FRICTIONLESS_DISPL_CONSTR;
            } else {
                logger::log_assert(false, logger::ERROR, "unknown contact type: {}", type);
            }
        }
        if(contact_type == FiniteElement::ContactType::FRICTIONLESS_PENALTY){
            this->log_data(doc, "opt_weight", projspec::TYPE_DOUBLE, true);
            EPS_DISPL = doc["opt_weight"].asDouble();
            if(this->log_data(doc, "max_step", projspec::TYPE_DOUBLE, false)){
                max_step = doc["max_step"].asDouble();
            }
        } else if(contact_type == FiniteElement::ContactType::FRICTIONLESS_DISPL_LOG){
            this->log_data(doc, "step_tol", projspec::TYPE_DOUBLE, true);
            step_tol = doc["step_tol"].asDouble();
            if(this->log_data(doc, "max_step", projspec::TYPE_DOUBLE, false)){
                max_step = doc["max_step"].asDouble();
            }
        } else if(contact_type >= FiniteElement::ContactType::FRICTIONLESS_DISPL_SIMPLE){
            this->log_data(doc, "step_tol", projspec::TYPE_DOUBLE, true);
            this->log_data(doc, "opt_weight", projspec::TYPE_DOUBLE, true);
            step_tol = doc["step_tol"].asDouble();
            EPS_DISPL = doc["opt_weight"].asDouble();
            if(this->log_data(doc, "max_step", projspec::TYPE_DOUBLE, false)){
                max_step = doc["max_step"].asDouble();
            }
        }
        return {contact_type, rtol_abs, max_step, step_tol, EPS_DISPL};
    }
    return FiniteElement::ContactData{FiniteElement::ContactType::RIGID};
}

std::vector<utils::DelayedPointer<Material>> ProjectData::generate_material_stubs(const Json::Value& doc,
                                const projspec::RequirementConditions& conds){
    std::vector<utils::DelayedPointer<Material>> vec;
    auto& mats = doc[Material::get_name()];
    vec.reserve(mats.size());
    for(const auto& it:mats){
        projspec::ObjectRequirements reqs;
        this->log_data(it, "type", projspec::TYPE_STRING, true);
        this->log_data(it, "name", projspec::TYPE_STRING, true);
        auto type_name = it["type"].asString();
        logger::log_assert(projspec::Registry::get(Material::get_name(), type_name, reqs),
                logger::ERROR, "unknown material type: {}", type_name);
        projspec::DataMap stub_data(this->generate_stub(reqs, conds));
        stub_data.set_string("name", it["name"].asString());
        stub_data.set_bool("STUB", true);
        vec.push_back(projspec::Factory<Material>::construct(type_name, stub_data));
    }

    return vec;
}

void ProjectData::load_materials(const Json::Value& doc, std::vector<utils::DelayedPointer<Material>>& vec,
                                const projspec::RequirementConditions& conds){

    auto& mats = doc[Material::get_name()];
    size_t i = 0;
    for(const auto& it:mats){
        projspec::ObjectRequirements reqs;
        this->log_data(it, "type", projspec::TYPE_STRING, true);
        auto type_name = it["type"].asString();
        logger::log_assert(projspec::Registry::get(Material::get_name(), type_name, reqs),
                logger::ERROR, "unknown material type: {}", type_name);
        projspec::DataMap data(this->generate_data_map(it, reqs, conds));
        vec[i] = projspec::Factory<Material>::construct(type_name, data);
        ++i;
    }
}

std::vector<utils::DelayedPointer<Field>> ProjectData::generate_field_stubs(const Json::Value& doc,
                                const projspec::RequirementConditions& conds){
    std::vector<utils::DelayedPointer<Field>> vec;
    auto& json_data = doc[Field::get_name()];
    vec.reserve(json_data.size());
    for(const auto& it:json_data){
        projspec::ObjectRequirements reqs;
        this->log_data(it, "type", projspec::TYPE_STRING, true);
        auto type_name = it["type"].asString();
        logger::log_assert(projspec::Registry::get(Field::get_name(), type_name, reqs),
                logger::ERROR, "unknown field type: {}", type_name);
        projspec::DataMap stub_data(this->generate_stub(reqs, conds));
        stub_data.set_bool("STUB", true);
        vec.push_back(projspec::Factory<Field>::construct(type_name, stub_data));
    }

    return vec;
}

void ProjectData::load_fields(const Json::Value& doc, std::vector<utils::DelayedPointer<Field>>& vec,
                                const projspec::RequirementConditions& conds){

    auto& json_data = doc[Field::get_name()];
    size_t i = 0;
    for(const auto& it:json_data){
        projspec::ObjectRequirements reqs;
        this->log_data(it, "type", projspec::TYPE_STRING, true);
        auto type_name = it["type"].asString();
        logger::log_assert(projspec::Registry::get(Field::get_name(), type_name, reqs),
                logger::ERROR, "unknown field type: {}", type_name);
        projspec::DataMap data(this->generate_data_map(it, reqs, conds));
        vec[i] = projspec::Factory<Field>::construct(type_name, data);
        ++i;
    }
}

std::vector<std::unique_ptr<Geometry>> ProjectData::load_geometries(const Json::Value& doc){
    const auto& geometries = doc["geometry"];
    std::vector<std::unique_ptr<Geometry>> geometry;
    size_t id = 0;
    for(const auto& geom:geometries){
        std::string absolute_path = folder_path;
        std::string geom_path = geom["file_path"].asString();
        absolute_path.append(geom_path);

        this->log_data(geom, "do_topopt", projspec::TYPE_BOOL, true);
        this->log_data(geom, "material", projspec::TYPE_STRING, true);

        double scale = 1;
        if(this->log_data(geom, "scale", projspec::TYPE_DOUBLE, false)){
            scale = geom["scale"].asDouble();
        }

        bool do_topopt = geom["do_topopt"].asBool();

        size_t alt_size = 0;
        if(this->log_data(geom, "alt_materials", projspec::TYPE_ARRAY, false)){
            const auto& alt = geom["alt_materials"];
            alt_size = alt.size();
        }

        bool with_void = ((alt_size > 0) ? false : true) && do_topopt;
        if(this->log_data(geom, "with_void", projspec::TYPE_BOOL, false) && alt_size > 0 && do_topopt){
            with_void = geom["with_void"].asBool();
        }
        if(this->do_simulation){
            do_topopt = false;
            with_void = false;
        }

        geometry.emplace_back(new Geometry(absolute_path, scale, this->type, this->topopt_element, do_topopt, with_void, id));

        std::string mat_name(geom["material"].asString());
        auto equal_name = [&mat_name](const utils::DelayedPointer<Material>& m)->bool{
            return mat_name == m->material_name;
        };
        auto it = std::find_if(this->materials.begin(), this->materials.end(), equal_name);
        logger::log_assert(it != this->materials.end(), logger::ERROR, "material with name '{}' not found", mat_name);

        std::vector<utils::DelayedPointerView<Material>> mats;
        mats.reserve(1 + alt_size);
        mats.push_back(it->get_view());

        if(this->log_data(geom, "alt_materials", projspec::TYPE_ARRAY, false)){
            const auto& alt = geom["alt_materials"];
            if(alt.size() > 0){
                for(const auto& mat:alt){
                    logger::log_assert(mat.isString(), logger::ERROR, "alt_materials must only contain the names of materials");
                    std::string mat_name(mat.asString());
                    auto equal_name = [&mat_name](const utils::DelayedPointer<Material>& m)->bool{
                        return mat_name == std::string(m->material_name);
                    };
                    auto it = std::find_if(this->materials.begin(), this->materials.end(), equal_name);
                    logger::log_assert(it != this->materials.end(), logger::ERROR, "material with name '{}' not found", mat_name);
                    mats.push_back(it->get_view());
                }
            }
        }
        geometry.back()->set_materials(std::move(mats));

        ++id;
    }
    return geometry;
}

std::unique_ptr<Pathfinding> ProjectData::load_pathfinder(const Json::Value& data,
                                            const projspec::RequirementConditions& conds){
    auto& json_data = data[Pathfinding::get_name()];
    projspec::ObjectRequirements reqs;
    this->log_data(json_data, "type", projspec::TYPE_STRING, true);
    auto type_name = json_data["type"].asString();
    logger::log_assert(projspec::Registry::get(Pathfinding::get_name(), type_name, reqs),
            logger::ERROR, "unknown pathfinder: {}", type_name);
    projspec::DataMap data_map(this->generate_data_map(json_data, reqs, conds));
    return projspec::Factory<Pathfinding>::construct(type_name, data_map);
}

std::unique_ptr<Sizing> ProjectData::load_sizer(const Json::Value& data,
                                            const projspec::RequirementConditions& conds){
    auto& json_data = data[Sizing::get_name()];
    projspec::ObjectRequirements reqs;
    this->log_data(json_data, "type", projspec::TYPE_STRING, true);
    auto type_name = json_data["type"].asString();
    logger::log_assert(projspec::Registry::get(Sizing::get_name(), type_name, reqs),
            logger::ERROR, "unknown pathfinder: {}", type_name);
    projspec::DataMap data_map(this->generate_data_map(json_data, reqs, conds));
    return projspec::Factory<Sizing>::construct(type_name, data_map);
}

std::unique_ptr<FiniteElement> ProjectData::load_fea(const Json::Value& data,
                                            const projspec::RequirementConditions& conds){
    auto& json_data = data[FiniteElement::get_name()];
    projspec::ObjectRequirements reqs;
    this->log_data(json_data, "type", projspec::TYPE_STRING, true);
    auto type_name = json_data["type"].asString();
    logger::log_assert(projspec::Registry::get(FiniteElement::get_name(), type_name, reqs),
            logger::ERROR, "unknown finite element solver: {}", type_name);
    projspec::DataMap data_map(this->generate_data_map(json_data, reqs, conds));
    return projspec::Factory<FiniteElement>::construct(type_name, data_map);
}

std::unique_ptr<Meshing> ProjectData::load_mesher(const Json::Value& data,
                                                  const projspec::RequirementConditions& conds){

    auto& json_data = data["mesher"];
    bool has_mesh_file_path = this->log_data(json_data , "file_path", projspec::TYPE_STRING, !this->do_meshing);
    if(has_mesh_file_path){
        std::string absolute_path = this->folder_path;
        std::string mesh_path = json_data["file_path"].asString();
        absolute_path.append(mesh_path);
        this->mesh_file_internal = std::make_unique<meshing::MeshFile>(this->geometries,
                this->topopt_element, this, this->thickness, 
                absolute_path, this->element_name);

        this->mesh_file = this->mesh_file_internal.get();
        // If meshing won't be done, use MeshFile instance as "mesher", so it
        // can load the raw mesh data and proceed from there.
        // As this is using unique_ptr, the mesh_file var will still
        // point to the correct place in memory.
        if(!this->do_meshing){
            return std::move(this->mesh_file_internal);
        } else {
            projspec::ObjectRequirements reqs;
            this->log_data(json_data, "type", projspec::TYPE_STRING, true);
            auto type_name = json_data["type"].asString();
            logger::log_assert(projspec::Registry::get(Meshing::get_name(), type_name, reqs),
                    logger::ERROR, "unknown mesher: {}", type_name);
            projspec::DataMap data_map(this->generate_data_map(json_data, reqs, conds));
            return projspec::Factory<Meshing>::construct(type_name, data_map);
        }
    }

    return nullptr;
}

std::unique_ptr<Projection> ProjectData::load_projection(const Json::Value& data,
                                                         const projspec::RequirementConditions& conds){
    if(!this->log_data(data, "projection", projspec::TYPE_OBJECT, false)){
        //return std::make_unique<projection::None>();
    }
    auto& json_data = data[Projection::get_name()];
    projspec::ObjectRequirements reqs;
    this->log_data(json_data, "type", projspec::TYPE_STRING, true);
    auto type_name = json_data["type"].asString();
    logger::log_assert(projspec::Registry::get(Projection::get_name(), type_name, reqs),
            logger::ERROR, "unknown projection method: {}", type_name);
    projspec::DataMap data_map(this->generate_data_map(json_data, reqs, conds));
    return projspec::Factory<Projection>::construct(type_name, data_map);
}

std::unique_ptr<DensityFilter> ProjectData::load_density_filter(const Json::Value& data,
                                                             const projspec::RequirementConditions& conds){

    auto& json_data = data[DensityFilter::get_name()];
    projspec::ObjectRequirements reqs;
    this->log_data(json_data, "type", projspec::TYPE_STRING, true);
    auto type_name = json_data["type"].asString();
    logger::log_assert(projspec::Registry::get(DensityFilter::get_name(), type_name, reqs),
            logger::ERROR, "unknown projection method: {}", type_name);
    projspec::DataMap data_map(this->generate_data_map(json_data, reqs, conds));
    return projspec::Factory<DensityFilter>::construct(type_name, data_map);
}

std::unique_ptr<DensityBasedOptimizer> ProjectData::load_topopt_optimizer(const Json::Value& data,
                                            const projspec::RequirementConditions& conds){
    auto& json_data = data[DensityBasedOptimizer::get_name()];
    projspec::ObjectRequirements reqs;
    if(this->log_data(json_data, "projection", projspec::TYPE_OBJECT, false)){
        auto& proj_data = json_data["projection"];
        this->log_data(proj_data, "beta", projspec::TYPE_OBJECT, true);
        this->projection_parameters = this->get_projection_parameter(proj_data["beta"]);

        this->projection = this->load_projection(json_data, conds);
    } else {
        this->projection = std::make_unique<projection::None>();
    }
    if(this->log_data(json_data, "density_filter", projspec::TYPE_OBJECT, true)){
        this->density_filter = this->load_density_filter(json_data, conds);
    }

    std::vector<std::unique_ptr<DensityBasedFunction>> objective;
    std::vector<double> weights;
    this->get_objective_functions(json_data, conds, objective, weights);

    std::vector<DensityBasedConstraint> constraints;
    this->get_constraints(json_data, conds, constraints);

    this->log_data(json_data, "type", projspec::TYPE_STRING, true);
    auto type_name = json_data["type"].asString();
    logger::log_assert(projspec::Registry::get(DensityBasedOptimizer::get_name(), type_name, reqs),
            logger::ERROR, "unknown topopt optimizer: {}", type_name);
    projspec::DataMap data_map(this->generate_data_map(json_data, reqs, conds));
    auto ptr(projspec::Factory<DensityBasedOptimizer>::construct(type_name, data_map));

    ptr->set_objective(std::move(objective), std::move(weights));
    ptr->set_constraints(std::move(constraints));

    return ptr;
}

std::unique_ptr<NodeShapeBasedOptimizer> ProjectData::load_shopt_optimizer(const Json::Value& data,
                                                                           const projspec::RequirementConditions& conds){
    auto& json_data = data[NodeShapeBasedOptimizer::get_name()];
    projspec::ObjectRequirements reqs;
   
    std::vector<Geometry*> geoms;
    geoms.reserve(this->geometries.size());
    for(auto& g:this->geometries){
        geoms.push_back(g.get());
    }

    this->log_data(json_data, "shape_set", projspec::TYPE_OBJECT, true); 
    auto shape_set = this->get_shape_operations(json_data["shape_set"]);
    this->shape_handler = ShapeHandler(this->topopt_mesher.get(), std::move(geoms), std::move(shape_set));

    std::vector<std::unique_ptr<NodeShapeBasedFunction>> objective;
    std::vector<double> weights;
    this->get_objective_functions(json_data, conds, objective, weights);

    std::vector<NodeShapeBasedConstraint> constraints;
    this->get_constraints(json_data, conds, constraints);

    this->log_data(json_data, "type", projspec::TYPE_STRING, true);
    auto type_name = json_data["type"].asString();
    logger::log_assert(projspec::Registry::get(NodeShapeBasedOptimizer::get_name(), type_name, reqs),
            logger::ERROR, "unknown topopt optimizer: {}", type_name);
    projspec::DataMap data_map(this->generate_data_map(json_data, reqs, conds));
    auto ptr(projspec::Factory<NodeShapeBasedOptimizer>::construct(type_name, data_map));

    ptr->set_objective(std::move(objective), std::move(weights));
    ptr->set_constraints(std::move(constraints));

    return ptr;
}

Projection::Parameter ProjectData::get_projection_parameter(const Json::Value& doc) const{
    this->log_data(doc, "initial", projspec::TYPE_DOUBLE, true);
    this->log_data(doc, "final", projspec::TYPE_DOUBLE, true);
    this->log_data(doc, "value_step", projspec::TYPE_DOUBLE, true);
    this->log_data(doc, "iteration_step", projspec::TYPE_INT, true);

    Projection::Parameter param;
    param.value = doc["initial"].asDouble();
    param.final_value = doc["final"].asDouble();
    param.value_step = doc["value_step"].asDouble();
    param.iteration_step = doc["iteration_step"].asInt();

    return param;
}

std::unique_ptr<MeshElementFactory> ProjectData::get_element_type(const std::string& name){
    this->element_name = name;
    if(name == "GT9"){
        return std::unique_ptr<MeshElementFactory>(static_cast<MeshElementFactory*>(
                    new MeshElementFactoryImpl<element::GT9>()
                ));
    } else if(name == "TRI3"){
        return std::unique_ptr<MeshElementFactory>(static_cast<MeshElementFactory*>(
                    new MeshElementFactoryImpl<element::TRI3>()
                ));
    } else if(name == "Q4"){
        return std::unique_ptr<MeshElementFactory>(static_cast<MeshElementFactory*>(
                    new MeshElementFactoryImpl<element::Q4>()
                ));
    } else if(name == "Q4S"){
        logger::log_assert(!this->do_shape_opt, logger::ERROR, "Q4S not supported for shape optimization, as it assumes element shapes are rectangular");
        return std::unique_ptr<MeshElementFactory>(static_cast<MeshElementFactory*>(
                    new MeshElementFactoryImpl<element::Q4S>()
                ));
    } else if(name == "TET4"){
        return std::unique_ptr<MeshElementFactory>(static_cast<MeshElementFactory*>(
                    new MeshElementFactoryImpl<element::TET4>()
                ));
    } else if(name == "H8"){
        return std::unique_ptr<MeshElementFactory>(static_cast<MeshElementFactory*>(
                    new MeshElementFactoryImpl<element::H8>()
                ));
    } else if(name == "TET10"){
        return std::unique_ptr<MeshElementFactory>(static_cast<MeshElementFactory*>(
                    new MeshElementFactoryImpl<element::TET10>()
                ));
    }
    return nullptr;
}

CrossSection ProjectData::get_cross_section(const Json::Value& doc) const{
    if(this->type == utils::PROBLEM_TYPE_2D){
        bool has_vertices = this->log_data(doc, "vertices", projspec::TYPE_ARRAY, false);
        bool has_point = this->log_data(doc, "point", projspec::TYPE_ARRAY, false);

        logger::log_assert(has_vertices || has_point, logger::ERROR, "invalid cross-section definition in configuration file.");
        if(has_vertices){
            auto vertices = doc["vertices"];

            std::vector<gp_Pnt> vlist;
            for(auto& v : vertices){
                logger::log_assert(v.size() == 2, logger::ERROR, "Vertices must have exactly two dimensions in 2D problems");
                vlist.emplace_back(v[0].asDouble(), v[1].asDouble(), 0);
            }
            return CrossSection(vlist, this->thickness);
        } else if(has_point){
            const auto& pa = doc["point"];
            logger::log_assert(pa.size() == 2, logger::ERROR, "Vertices must have exactly two dimensions in 2D problems");
            gp_Pnt p(pa[0].asDouble(), pa[1].asDouble(), 0);

            return CrossSection(p);
        }
    } else if(this->type == utils::PROBLEM_TYPE_3D) {
        bool has_rect = this->log_data(doc, "rectangle", projspec::TYPE_OBJECT, false);
        bool has_file = this->log_data(doc, "file", projspec::TYPE_OBJECT, false);
        bool has_point = this->log_data(doc, "point", projspec::TYPE_ARRAY, false);

        logger::log_assert(has_rect || has_file || has_point, logger::ERROR, "invalid cross-section definition in configuration file.");
        if(has_rect){
            const auto& rect = doc["rectangle"];
            this->log_data(rect, "center", projspec::TYPE_ARRAY, true);
            this->log_data(rect, "normal", projspec::TYPE_ARRAY, true);
            this->log_data(rect, "w", projspec::TYPE_DOUBLE, true);
            this->log_data(rect, "h", projspec::TYPE_DOUBLE, true);
            this->log_data(rect, "rotation", projspec::TYPE_DOUBLE, true);
            auto c = rect["center"];
            auto n = rect["normal"];
            CrossSection::Rectangle r{
                rect["w"].asDouble(),
                rect["h"].asDouble(),
                gp_Pnt(
                    c[0].asDouble(), 
                    c[1].asDouble(),
                    c[2].asDouble()
                ),
                gp_Dir(
                    n[0].asDouble(), 
                    n[1].asDouble(),
                    n[2].asDouble()
                ),
                rect["rotation"].asDouble()
            };
            return CrossSection(r);
        } else if(has_file) {
            logger::log_assert(!this->generate_beams,
                               logger::ERROR, "cross-section of type file is currently not fully compatible with beam generation");
            const auto& file = doc["file"];
            this->log_data(file, "path", projspec::TYPE_STRING, true);

            double scale = 1;
            if(this->log_data(file, "scale", projspec::TYPE_DOUBLE, false)){
                scale = file["scale"].asDouble();
            }

            std::string absolute_path = folder_path;
            std::string path = file["path"].asString();
            absolute_path.append(path);

            return CrossSection(absolute_path, scale);
        } else if(has_point){
            const auto& pa = doc["point"];
            logger::log_assert(pa.size() == 3, logger::ERROR, "Vertices must have exactly three dimensions in 3D problems");
            gp_Pnt p(pa[0].asDouble(), pa[1].asDouble(), pa[2].asDouble());

            return CrossSection(p);
        }
    }
    logger::log_assert(false, logger::ERROR, "unknown problem type detected in get_cross_section(), this shouldn't have happened.");
    return CrossSection(gp_Pnt(0,0,0));
}

std::vector<Force> ProjectData::get_loads(const Json::Value& doc) const{
    std::vector<Force> forces;
    if(this->type == utils::PROBLEM_TYPE_2D){
        for(auto& f : doc){
            logger::log_assert(f.isObject(), logger::ERROR, "Each load must be stored as a JSON object");
            this->log_data(f, "load", projspec::TYPE_ARRAY, true);
            this->log_data(f, "cut_into_shape", projspec::TYPE_BOOL, true);
            bool cut_into_shape = f["cut_into_shape"].asBool();

            auto loads = f["load"];
            this->log_data(loads, "region", projspec::TYPE_OBJECT, true);
            logger::log_assert(loads.size() == 2, logger::ERROR, "Load vector must have exactly two dimensions in 2D problems");

            gp_Vec l(loads[0].asDouble(), loads[1].asDouble(), 0);

            this->log_data(f, "region", projspec::TYPE_OBJECT, true);
            auto S = this->get_cross_section(f["region"]);
            forces.emplace_back(std::move(S), l, cut_into_shape);
        }
    } else if(this->type == utils::PROBLEM_TYPE_3D) {
        for(auto& f : doc){
            logger::log_assert(f.isObject(), logger::ERROR, "Each load must be stored as a JSON object");
            this->log_data(f, "load", projspec::TYPE_ARRAY, true);
            auto loads = f["load"];
            logger::log_assert(loads.size() == 3, logger::ERROR, "Load vector must have exactly three dimensions in 3D problems");
            this->log_data(f, "cut_into_shape", projspec::TYPE_BOOL, true);
            bool cut_into_shape = f["cut_into_shape"].asBool();

            gp_Vec l(loads[0].asDouble(), loads[1].asDouble(), loads[2].asDouble());

            this->log_data(f, "region", projspec::TYPE_OBJECT, true);
            auto S = this->get_cross_section(f["region"]);
            forces.emplace_back(std::move(S), l, cut_into_shape);
        }
    }
    return forces;
}

std::vector<Support> ProjectData::get_support(const Json::Value& doc) const{
    std::vector<Support> supports;
    if(this->type == utils::PROBLEM_TYPE_2D){
        for(auto& f : doc){
            logger::log_assert(f.isObject(), logger::ERROR, "Each support must be stored as a JSON object");
            this->log_data(f, "X", projspec::TYPE_BOOL, true);
            this->log_data(f, "Y", projspec::TYPE_BOOL, true);
            this->log_data(f, "MZ", projspec::TYPE_BOOL, true);
            bool X = f["X"].asBool();
            bool Y = f["Y"].asBool();
            bool MZ = f["MZ"].asBool();
            this->log_data(f, "cut_into_shape", projspec::TYPE_BOOL, true);
            bool cut_into_shape = f["cut_into_shape"].asBool();

            this->log_data(f, "region", projspec::TYPE_OBJECT, true);
            auto S = this->get_cross_section(f["region"]);
            supports.emplace_back(X, Y, MZ, S, cut_into_shape);
        }
    } else if(this->type == utils::PROBLEM_TYPE_3D) {
        for(auto& f : doc){
            logger::log_assert(f.isObject(), logger::ERROR, "Each support must be stored as a JSON object");
            this->log_data(f, "X", projspec::TYPE_BOOL, true);
            this->log_data(f, "Y", projspec::TYPE_BOOL, true);
            this->log_data(f, "Z", projspec::TYPE_BOOL, true);
            this->log_data(f, "MX", projspec::TYPE_BOOL, true);
            this->log_data(f, "MY", projspec::TYPE_BOOL, true);
            this->log_data(f, "MZ", projspec::TYPE_BOOL, true);
            bool X = f["X"].asBool();
            bool Y = f["Y"].asBool();
            bool Z = f["Z"].asBool();
            bool MX = f["MX"].asBool();
            bool MY = f["MY"].asBool();
            bool MZ = f["MZ"].asBool();
            this->log_data(f, "cut_into_shape", projspec::TYPE_BOOL, true);
            bool cut_into_shape = f["cut_into_shape"].asBool();

            this->log_data(f, "region", projspec::TYPE_OBJECT, true);
            auto S = this->get_cross_section(f["region"]);
            supports.emplace_back(X, Y, Z, MX, MY, MZ, S, cut_into_shape);
        }
    }
    return supports;
}

std::vector<Spring> ProjectData::get_springs(const Json::Value& doc){
    std::vector<Spring> springs;
    for(auto& f : doc){
        logger::log_assert(f.isObject(), logger::ERROR, "Each load must be stored as a JSON object");
        this->log_data(f, "L", projspec::TYPE_ARRAY, true);
        this->log_data(f, "normal", projspec::TYPE_ARRAY, true);
        this->log_data(f, "v", projspec::TYPE_ARRAY, true);
        if(this->type == utils::PROBLEM_TYPE_3D) {
            this->log_data(f, "w", projspec::TYPE_ARRAY, true);
        }
        this->log_data(f, "material", projspec::TYPE_STRING, true);
        auto L = f["L"];
        auto normal = f["normal"];
        auto v = f["v"];
        std::string mat_name(f["material"].asString());

        std::array<double, 3> l{L[0].asDouble(), L[1].asDouble(), 0};
        gp_Dir nv(normal[0].asDouble(), normal[1].asDouble(), 0);
        gp_Dir vv(v[0].asDouble(), v[1].asDouble(), 0);
        gp_Dir wv(0, 0, 1);

        if(this->type == utils::PROBLEM_TYPE_2D){
            logger::log_assert(L.size() == 2, logger::ERROR, "Length vector must have exactly two dimensions in 2D problems");
            logger::log_assert(normal.size() == 2, logger::ERROR, "Normal vector must have exactly two dimensions in 2D problems");
            logger::log_assert(v.size() == 2, logger::ERROR, "'v' vector must have exactly two dimensions in 2D problems");
        } else if(this->type == utils::PROBLEM_TYPE_3D) {
            logger::log_assert(L.size() == 3, logger::ERROR, "Length vector must have exactly three dimensions in 3D problems");
            logger::log_assert(normal.size() == 3, logger::ERROR, "Normal vector must have exactly three dimensions in 3D problems");
            logger::log_assert(v.size() == 3, logger::ERROR, "'v' vector must have exactly three dimensions in 3D problems");

            auto w = f["w"];
            logger::log_assert(w.size() == 3, logger::ERROR, "'w' vector must have exactly two dimensions in 3D problems");

            l[2] = L[2].asDouble();
            nv.SetZ(normal[2].asDouble());
            vv.SetZ(v[2].asDouble());

            wv = gp_Dir(w[0].asDouble(), w[1].asDouble(), w[2].asDouble());
        }

        auto mat(this->get_material(mat_name));

        this->log_data(f, "region", projspec::TYPE_OBJECT, true);
        auto S = this->get_cross_section(f["region"]);
        springs.emplace_back(S, this->thickness, nv, vv, wv, mat, l, this->topopt_element, this->type);
    }
    return springs;
}

std::vector<InternalLoads> ProjectData::get_internal_loads(const Json::Value& doc){
    std::vector<InternalLoads> internal_loads;
    for(auto& f : doc){
        logger::log_assert(f.isObject(), logger::ERROR, "Each load must be stored as a JSON object");
        this->log_data(f, "F", projspec::TYPE_ARRAY, true);
        this->log_data(f, "M", projspec::TYPE_ARRAY, true);
        this->log_data(f, "normal", projspec::TYPE_ARRAY, true);
        this->log_data(f, "v", projspec::TYPE_ARRAY, true);
        if(this->type == utils::PROBLEM_TYPE_3D) {
            this->log_data(f, "w", projspec::TYPE_ARRAY, true);
        }
        this->log_data(f, "material", projspec::TYPE_STRING, true);
        auto normal = f["normal"];
        auto v = f["v"];
        auto Fa = f["F"];
        auto Ma = f["M"];
        std::string mat_name(f["material"].asString());

        std::array<double, 3> F{Fa[0].asDouble(), Fa[1].asDouble(), 0};
        std::array<double, 3> M{Ma[0].asDouble(), Ma[1].asDouble(), 0};
        gp_Dir nv(normal[0].asDouble(), normal[1].asDouble(), 0);
        gp_Dir vv(v[0].asDouble(), v[1].asDouble(), 0);
        gp_Dir wv(0, 0, 1);

        if(this->type == utils::PROBLEM_TYPE_2D){
            logger::log_assert(normal.size() == 2, logger::ERROR, "Normal vector must have exactly two dimensions in 2D problems");
            logger::log_assert(v.size() == 2, logger::ERROR, "'v' vector must have exactly two dimensions in 2D problems");
            logger::log_assert(Fa.size() == 2, logger::ERROR, "Force vector must have exactly two dimensions in 2D problems");
            logger::log_assert(Ma.size() == 2, logger::ERROR, "Moment vector must have exactly two dimensions in 2D problems");
        } else if(this->type == utils::PROBLEM_TYPE_3D) {
            logger::log_assert(normal.size() == 3, logger::ERROR, "Normal vector must have exactly three dimensions in 3D problems");
            logger::log_assert(v.size() == 3, logger::ERROR, "'v' vector must have exactly three dimensions in 3D problems");
            logger::log_assert(Fa.size() == 3, logger::ERROR, "Force vector must have exactly three dimensions in 3D problems");
            logger::log_assert(Ma.size() == 3, logger::ERROR, "Moment vector must have exactly three dimensions in 3D problems");

            auto w = f["w"];
            logger::log_assert(w.size() == 3, logger::ERROR, "'w' vector must have exactly two dimensions in 3D problems");

            nv.SetZ(normal[2].asDouble());
            vv.SetZ(v[2].asDouble());
            F[2] = Fa[2].asDouble();
            M[2] = Ma[2].asDouble();

            wv = gp_Dir(w[0].asDouble(), w[1].asDouble(), w[2].asDouble());
        }

        auto mat(this->get_material(mat_name));

        this->log_data(f, "region", projspec::TYPE_OBJECT, true);
        auto S = this->get_cross_section(f["region"]);
        internal_loads.emplace_back(S, this->thickness, nv, vv, wv, mat, F, M, this->topopt_element, this->type);
    }
    return internal_loads;
}

std::vector<SubProblem> ProjectData::load_sub_problems(const Json::Value& doc){
    const auto& sub_problem_array = doc["sub_problems"];
    std::vector<SubProblem> sub_problems;
    for(const auto& sp:sub_problem_array){
        logger::log_assert(sp.isObject(), logger::ERROR, "\"sub_problems\" must be an array of objects");
        SubProblem curr_sp;
        if(this->log_data(sp, "loads", projspec::TYPE_ARRAY, false)){
            for(auto& f:sp["loads"]){
                curr_sp.forces.push_back(&this->forces[f.asInt()]);
            }
        }
        if(this->log_data(sp, "supports", projspec::TYPE_ARRAY, false)){
            for(auto& f:sp["supports"]){
                curr_sp.supports.push_back(&this->supports[f.asInt()]);
            }
        }
        if(this->log_data(sp, "springs", projspec::TYPE_ARRAY, false)){
            for(auto& f:sp["springs"]){
                curr_sp.springs.push_back(&this->springs[f.asInt()]);
            }
        }
        if(this->log_data(sp, "internal_loads", projspec::TYPE_ARRAY, false)){
            for(auto& f:sp["internal_loads"]){
                curr_sp.internal_loads.push_back(&this->internal_loads[f.asInt()]);
            }
        }
        sub_problems.push_back(std::move(curr_sp));
    }

    return sub_problems;
}

std::unique_ptr<shape_op::ShapeOp> ProjectData::get_shape_operations(const Json::Value& doc) const{
    this->log_data(doc, "type", projspec::TYPE_STRING, true);
    const std::string type = doc["type"].asString();
    if(type == "geometry"){
        this->log_data(doc, "id", projspec::TYPE_INT, true);
        const size_t id = doc["id"].asInt64();
        logger::log_assert(id < this->geometries.size(), logger::ERROR,
                "unknown geometry id: {}", id);
        return std::make_unique<shape_op::Geometry>(id);
    } else if(type == "union"){
        this->log_data(doc, "shape1", projspec::TYPE_OBJECT, true);
        this->log_data(doc, "shape2", projspec::TYPE_OBJECT, true);

        return std::make_unique<shape_op::Union>(
                    this->get_shape_operations(doc["shape1"]),
                    this->get_shape_operations(doc["shape2"])
                );
    } else if(type == "intersection"){
        this->log_data(doc, "shape1", projspec::TYPE_OBJECT, true);
        this->log_data(doc, "shape2", projspec::TYPE_OBJECT, true);

        return std::make_unique<shape_op::Intersection>(
                    this->get_shape_operations(doc["shape1"]),
                    this->get_shape_operations(doc["shape2"])
                );
    } else if(type == "difference"){
        this->log_data(doc, "shape1", projspec::TYPE_OBJECT, true);
        this->log_data(doc, "shape2", projspec::TYPE_OBJECT, true);

        return std::make_unique<shape_op::Difference>(
                    this->get_shape_operations(doc["shape1"]),
                    this->get_shape_operations(doc["shape2"])
                );
    }

    return nullptr;
}

std::unique_ptr<Simulation> ProjectData::load_simulation(const Json::Value& data,
                                                         const projspec::RequirementConditions& conds){
    auto& json_data = data[Simulation::get_name()];
    projspec::ObjectRequirements reqs;
    this->log_data(json_data, "type", projspec::TYPE_STRING, true);
    auto type_name = json_data["type"].asString();
    logger::log_assert(projspec::Registry::get(Simulation::get_name(), type_name, reqs),
            logger::ERROR, "unknown simulation: {}", type_name);
    projspec::DataMap data_map(this->generate_data_map(json_data, reqs, conds));
    return projspec::Factory<Simulation>::construct(type_name, data_map);
}
