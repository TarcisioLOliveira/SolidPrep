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
#include "material/mandible.hpp"
#include "meshing/mesh_file.hpp"
#include "shape_handler.hpp"
#include <cstring>
#include <memory>
#define _USE_MATH_DEFINES
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include "rapidjson/rapidjson.h"
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/error/error.h"
#include "rapidjson/error/en.h"

#include "utils.hpp"
#include "force.hpp"
#include "logger.hpp"
#include "project_data.hpp"
#include "pathfinding/meshless_astar.hpp"
#include "pathfinding/visibility_graph.hpp"
#include "material/linear_elastic_isotropic.hpp"
#include "material/linear_elastic_orthotropic.hpp"
#include "material/mandible.hpp"
#include "material/linear_elastic_orthotropic_field.hpp"
#include "sizing/standard_sizing.hpp"
//#include "finite_element/direct_solver.hpp"
//#include "finite_element/gradient_descent.hpp"
//#include "finite_element/PCG.hpp"
#include "finite_element/mumps_solver.hpp"
#include "finite_element/eigen_pcg.hpp"
#include "finite_element/petsc_pcg.hpp"
#include "meshing/gmsh.hpp"
#include "sizing/beam_sizing.hpp"
#include "element/GT9.hpp"
#include "element/TRI3.hpp"
#include "element/Q4.hpp"
#include "element/Q4S.hpp"
#include "element/TET4.hpp"
#include "element/H8.hpp"
#include "element/TET10.hpp"
#include "density_filter/convolution.hpp"
#include "density_filter/helmholtz.hpp"
#include "density_filter/averaging.hpp"
#include "projection/none.hpp"
#include "projection/threshold.hpp"
#include "projection/heaviside.hpp"
#include "optimizer/density_based/mma.hpp"
#include "optimizer/density_based/newton.hpp"
#include "optimizer/node_shape_based/mma.hpp"
#include "function/density_based/compliance.hpp"
#include "function/density_based/volume.hpp"
#include "function/density_based/global_stress_pnorm_normalized.hpp"
#include "function/density_based/global_stress_pnorm.hpp"
#include "function/density_based/global_stress_heaviside.hpp"
#include "function/density_based/omni_machining.hpp"
#include "function/density_based/am_support.hpp"
#include "function/density_based/mass.hpp"
#include "function/density_based/mass_first_material.hpp"
#include "function/density_based/mechanostat.hpp"
#include "function/node_shape_based/compliance.hpp"
#include "function/node_shape_based/volume.hpp"
#include "function/node_shape_based/global_stress_heaviside.hpp"
#include "function/node_shape_based/mechanostat.hpp"
#include "field/orthotropic_flow.hpp"
#include "field/principal_stress.hpp"
#include "simulation/marginal_bone_loss.hpp"

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

    this->folder_path = this->get_folder_path(project_file);
    if(this->log_data(doc, "contact_type", TYPE_OBJECT, false)){
        this->contact_data = this->get_contact_data(doc["contact_type"]);
    } else {
        this->contact_data = ProjectData::ContactData{FiniteElement::ContactType::RIGID, 0.0};
    }

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

    if(this->log_data(doc, "analysis", TYPE_OBJECT, true)){
        const auto& analysis = doc["analysis"];
        if(this->log_data(analysis, "simulation", TYPE_BOOL, false)){
            this->do_simulation = analysis["simulation"].GetBool();
        }
        if(!this->do_simulation){
            if(this->log_data(analysis, "gen", TYPE_BOOL, false)){
                this->generate_beams = analysis["gen"].GetBool();
            }
            if(this->log_data(analysis, "meshing", TYPE_BOOL, false)){
                this->do_meshing = analysis["meshing"].GetBool();
            } else if(this->generate_beams){
                this->do_meshing = true;
            }
            if(this->generate_beams){
                logger::log_assert(this->do_meshing, logger::WARNING, "mesh loading not supported for geometry generation, ignoring setting");
                this->do_meshing = true;
            }
            if(this->log_data(analysis, "topopt", TYPE_BOOL, false)){
                this->do_topopt = analysis["topopt"].GetBool();
            }
            if(this->log_data(analysis, "fea", TYPE_BOOL, false)){
                this->do_fea = analysis["fea"].GetBool();
            }
            if(this->log_data(analysis, "shape_opt", TYPE_BOOL, false)){
                this->do_shape_opt = analysis["shape_opt"].GetBool();
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
        if(this->log_data(doc, "thickness", TYPE_DOUBLE, true)){
            this->thickness = doc["thickness"].GetDouble();
        }
    }
    if(this->log_data(doc, "loads", TYPE_ARRAY, false)){
        this->forces = this->get_loads(doc["loads"]);
    }
    if(this->log_data(doc, "supports", TYPE_ARRAY, false)){
        this->supports = this->get_support(doc["supports"]);
    }
    if(this->log_data(doc, "mesher", TYPE_OBJECT, true)){
        this->log_data(doc["mesher"], "element_type", TYPE_STRING, true);
        if(this->do_meshing){
            this->topopt_element = this->get_element_type(doc["mesher"]["element_type"].GetString());
            this->topopt_boundary_element = this->topopt_element->get_boundary_element_info();
        } else {
            std::string absolute_path = this->folder_path;
            std::string mesh_path = doc["mesher"]["file_path"].GetString();
            absolute_path.append(mesh_path);

            std::ifstream file(absolute_path, std::ios::in);
            logger::log_assert(file.good(), logger::ERROR, "file not found: {}", absolute_path);
            std::getline(file, this->element_name);
            file.close();

            this->topopt_element = this->get_element_type(this->element_name);
            this->topopt_boundary_element = this->topopt_element->get_boundary_element_info();
        }
    }
    if(this->log_data(doc, "geometry", TYPE_ARRAY, true)){
        this->geometries = this->load_geometries(doc);
    }
    if(this->log_data(doc, "material", TYPE_ARRAY, true)){
        this->materials = this->load_materials_simple(doc);
    }
    if(this->log_data(doc, "fields", TYPE_ARRAY, false)){
        this->fields = this->load_fields(doc);
    }
    if(this->log_data(doc, "material", TYPE_ARRAY, true)){
        auto field_mat = this->load_materials_field(doc);
        const size_t old_size = this->materials.size();
        this->materials.resize(this->materials.size() + field_mat.size());
        std::move(field_mat.begin(), field_mat.end(), this->materials.begin() + old_size);
    }
    this->assign_materials(doc);
    if(this->log_data(doc, "springs", TYPE_ARRAY, false)){
        this->springs = this->get_springs(doc["springs"]);
    }
    logger::log_assert(this->supports.size() > 0 || this->springs.size() > 0, logger::ERROR,
                       "a support or spring must be specified to prevent matrix singularity.");
    if(this->log_data(doc, "internal_loads", TYPE_ARRAY, false)){
        this->internal_loads = this->get_internal_loads(doc["internal_loads"]);
    }
    if(this->log_data(doc, "sub_problems", TYPE_ARRAY, false)){
        this->sub_problems = this->load_sub_problems(doc);
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
    if(this->log_data(doc, "finite_element", TYPE_OBJECT, true)){
        std::vector<std::unique_ptr<FiniteElement>> solvers(this->sub_problems.size());
        // Dirty, but quicker to implement and I don't have to mess with
        // `sizer_fea`down there.
        for(size_t i = 0; i < this->sub_problems.size(); ++i){
            solvers[i] = this->load_fea(doc);
        }
        this->topopt_fea = std::make_unique<SolverManager>(std::move(solvers));
    }
    if(this->log_data(doc, "sizing", TYPE_OBJECT, this->generate_beams)){
        if(this->log_data(doc["sizing"], "pathfinding", TYPE_OBJECT, false)){
            this->pathfinder = this->load_pathfinder(doc["sizing"]);
        }
        if(this->log_data(doc["sizing"], "finite_element", TYPE_OBJECT, false)){
            this->sizer_fea = this->load_fea(doc["sizing"]);
        }
        this->sizer = this->load_sizer(doc);
    }
    if(this->log_data(doc, "mesher", TYPE_OBJECT, true)){
        this->topopt_mesher = this->load_mesher(doc);
    }
    if(this->log_data(doc, "topopt", TYPE_OBJECT, this->do_topopt)){
        this->topopt_optimizer = this->load_topopt_optimizer(doc);
    }
    if(this->log_data(doc, "shape_opt", TYPE_OBJECT, this->do_shape_opt)){
        this->shopt_optimizer = this->load_shopt_optimizer(doc);
    }
    if(this->log_data(doc, "simulation", TYPE_OBJECT, this->do_simulation)){
        this->simulator = this->load_simulation(doc["simulation"]);
    }


    fclose(fp);
}


std::string ProjectData::get_folder_path(const std::string& project_file_path) const{
#ifdef _WIN32
    size_t last_slash = project_file_path.rfind("\\");
#else
    size_t last_slash = project_file_path.rfind("/");
#endif

    return project_file_path.substr(0, last_slash+1);
}

bool ProjectData::log_data(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc, std::string name, ProjectData::DataType type, bool required) const{
    logger::AssertType error = (required) ? logger::ERROR : logger::SILENT;
    bool exists = logger::log_assert(doc.HasMember(name.c_str()), error, "Missing member: {}", name);
    if(!exists){
        return false;
    }
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
    return correct_type;
}

ProjectData::ContactData ProjectData::get_contact_data(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    if(this->log_data(doc, "type", TYPE_STRING, true)){
        const std::string type = doc["type"].GetString();
        FiniteElement::ContactType contact_type = FiniteElement::ContactType::RIGID;
        double rtol_abs = 0;
        double max_step = 0;
        double EPS_DISPL = 0;
        if(type == "rigid"){
            contact_type = FiniteElement::ContactType::RIGID;
        } else {
            this->log_data(doc, "rtol_abs", TYPE_DOUBLE, true);
            rtol_abs = doc["rtol_abs"].GetDouble();
            if(type == "frictionless_penalty"){
            contact_type = FiniteElement::ContactType::FRICTIONLESS_PENALTY;
            } else if(type == "frictionless_displ_simple"){
                contact_type = FiniteElement::ContactType::FRICTIONLESS_DISPL_SIMPLE;
            } else if(type == "frictionless_displ_constr"){
                contact_type = FiniteElement::ContactType::FRICTIONLESS_DISPL_CONSTR;
            } else {
                logger::log_assert(false, logger::ERROR, "unknown contact type: {}", type);
            }
        }
        if(contact_type == FiniteElement::ContactType::FRICTIONLESS_DISPL_SIMPLE){
            this->log_data(doc, "max_step", TYPE_DOUBLE, true);
            this->log_data(doc, "opt_weight", TYPE_DOUBLE, true);
            max_step = doc["max_step"].GetDouble();
            EPS_DISPL = doc["opt_weight"].GetDouble();
        }
        return {contact_type, rtol_abs, max_step, EPS_DISPL};
    }
    return ProjectData::ContactData{FiniteElement::ContactType::RIGID, 0.0};
}

std::vector<std::unique_ptr<Material>> ProjectData::load_materials_simple(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    const auto& materials = doc["material"].GetArray();
    std::vector<std::unique_ptr<Material>> material;

    std::vector<size_t> queue;
    for(const auto& mat:materials){
        logger::log_assert(mat.HasMember("name"), logger::ERROR, "missing material property: name");
        logger::log_assert(mat["name"].IsString(), logger::ERROR, "material property 'name' must be a string");
        std::string name(mat["name"].GetString());
        this->log_data(mat, "type", TYPE_STRING, true);
        if(mat["type"] == "linear_elastic_orthotropic"){
            std::vector<std::string> properties{"E", "nu", "G", "Smax", "Tmax"};
            for(auto& s:properties){
                logger::log_assert(mat.HasMember(s.c_str()), logger::ERROR, "missing material property: {}", s);
                logger::log_assert(mat[s.c_str()].IsArray() || mat[s.c_str()].IsDouble(), logger::ERROR, "material property {} must be either a number or an array of numbers", s);
            }
            std::vector<std::vector<double>> values(5);
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
            this->log_data(mat, "density", TYPE_DOUBLE, true);
            double density = mat["density"].GetDouble();
            this->log_data(mat, "nu_lower_half", TYPE_ARRAY, true);
            const auto& nuarr = mat["nu_lower_half"].GetArray();
            std::vector<bool> nu_lower_half
                {nuarr[0].GetBool(),
                 nuarr[1].GetBool(),
                 nuarr[2].GetBool()};

            for(auto& i:values[0]) i *= 1e3; // E
            for(auto& i:values[2]) i *= 1e3; // G
            // for(auto& i:values[0]) i *= 1e9; // E
            // for(auto& i:values[2]) i *= 1e9; // G
            // for(auto& i:values[3]) i *= 1e6; // Smax
            // for(auto& i:values[4]) i *= 1e6; // Tmax
            material.emplace_back(new material::LinearElasticOrthotropic(name, density, values[0], values[1], nu_lower_half, values[2], values[3], values[4]));
        } else if(mat["type"] == "linear_elastic_isotropic"){
            std::vector<std::string> properties{"density", "E", "nu", "Smax", "Tmax"};
            for(auto& s:properties){
                this->log_data(mat, s, TYPE_DOUBLE, true);
            }
            this->log_data(mat, "plane_stress", TYPE_BOOL, true);
            double E = mat["E"].GetDouble();
            double nu = mat["nu"].GetDouble();
            double Smax = mat["Smax"].GetDouble();
            double Tmax = mat["Tmax"].GetDouble();
            double density = mat["density"].GetDouble();
            bool plane_stress = mat["plane_stress"].GetBool();

            material.emplace_back(new material::LinearElasticIsotropic(name, density, E*1e3, nu, Smax, Tmax, plane_stress));
        }
    }

    return material;
}
std::vector<std::unique_ptr<Material>> ProjectData::load_materials_field(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    const auto& materials = doc["material"].GetArray();
    std::vector<std::unique_ptr<Material>> material;

    std::vector<size_t> queue;
    size_t mat_num = 0;
    for(const auto& mat:materials){
        logger::log_assert(mat.HasMember("name"), logger::ERROR, "missing material property: name");
        logger::log_assert(mat["name"].IsString(), logger::ERROR, "material property 'name' must be a string");
        std::string name(mat["name"].GetString());
        this->log_data(mat, "type", TYPE_STRING, true);
        if(mat["type"] == "linear_elastic_orthotropic_field"){
            std::vector<std::string> properties{"E", "nu", "G", "Smax", "Tmax"};
            for(auto& s:properties){
                logger::log_assert(mat.HasMember(s.c_str()), logger::ERROR, "missing material property: {}", s);
                logger::log_assert(mat[s.c_str()].IsArray() || mat[s.c_str()].IsDouble(), logger::ERROR, "material property {} must be either a number or an array of numbers", s);
            }
            std::vector<std::vector<double>> values(5);
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
            this->log_data(mat, "field", TYPE_INT, true);
            this->log_data(mat, "density", TYPE_DOUBLE, true);
            double density = mat["density"].GetDouble();
            int field_num = mat["field"].GetInt();
            this->log_data(mat, "nu_lower_half", TYPE_ARRAY, true);
            const auto& nuarr = mat["nu_lower_half"].GetArray();
            std::vector<bool> nu_lower_half
                {nuarr[0].GetBool(),
                 nuarr[1].GetBool(),
                 nuarr[2].GetBool()};

            for(auto& i:values[0]) i *= 1e3; // E
            for(auto& i:values[2]) i *= 1e3; // G
            // for(auto& i:values[0]) i *= 1e9; // E
            // for(auto& i:values[2]) i *= 1e9; // G
            // for(auto& i:values[3]) i *= 1e6; // Smax
            // for(auto& i:values[4]) i *= 1e6; // Tmax

            logger::log_assert(field_num >= 0, logger::ERROR, "field id must be non-negative");
            logger::log_assert((size_t)field_num < this->fields.size(), logger::ERROR, "field id is greater or equal to number of defined fields");
            Field* f = this->fields[field_num].get();
            logger::log_assert(f->get_type() == Field::Type::COORDINATE, logger::ERROR, "orthotropic_field material type requires a coordinate field");
            CoordinateField* cf = static_cast<CoordinateField*>(f);

            material.emplace_back(new material::LinearElasticOrthotropicField(name, density, values[0], values[1], nu_lower_half, values[2], values[3], values[4], cf));
        } else if(mat["type"] == "mandible"){
            queue.push_back(mat_num);
        }
        ++mat_num;
    }

    // Materials that depend on other materials.
    for(auto& q:queue){
        const auto& mat = materials[q];
        std::string name(mat["name"].GetString());
        if(mat["type"] == "mandible"){
            std::vector<std::string> prop_strings{"material_inner", "material_outer", "path_points1", "path_points2"};
            for(auto& s:prop_strings){
                this->log_data(mat, s, TYPE_STRING, true);
            }
            this->log_data(mat, "C", TYPE_DOUBLE, true);
            std::string inner = mat["material_inner"].GetString();
            std::string outer = mat["material_outer"].GetString();
            std::string path_points1 = this->folder_path;
            path_points1.append(mat["path_points1"].GetString());
            std::string path_points2 = this->folder_path;
            path_points2.append(mat["path_points2"].GetString());
            double C = mat["C"].GetDouble();
            auto equal_name_inner = [&inner](const std::unique_ptr<Material>& m)->bool{
                return inner == m->name;
            };
            auto equal_name_outer = [&outer](const std::unique_ptr<Material>& m)->bool{
                return outer == m->name;
            };
            Material* i = nullptr;
            Material* o = nullptr;

            auto it1 = std::find_if(material.begin(), material.end(), equal_name_inner);
            auto it2 = std::find_if(this->materials.begin(), this->materials.end(), equal_name_inner);
            logger::log_assert(it1 != material.end() || it2 != this->materials.end(), logger::ERROR, "material with name '{}' not found", inner);
            if(it1 != material.end()){
                i = it1->get();
            } else {
                i = it2->get();
            }

            it1 = std::find_if(material.begin(), material.end(), equal_name_outer);
            it2 = std::find_if(this->materials.begin(), this->materials.end(), equal_name_outer);
            logger::log_assert(it1 != material.end() || it2 != this->materials.end(), logger::ERROR, "material with name '{}' not found", inner);
            if(it1 != material.end()){
                o = it1->get();
            } else {
                o = it2->get();
            }

            bool has_implant = this->log_data(mat, "implant", TYPE_OBJECT, false);
            if(has_implant){
                const auto& impdata = mat["implant"];
                this->log_data(impdata, "center1", TYPE_ARRAY, true);
                this->log_data(impdata, "center2", TYPE_ARRAY, true);
                this->log_data(impdata, "r1", TYPE_DOUBLE, true);
                this->log_data(impdata, "r2", TYPE_DOUBLE, true);

                this->log_data(impdata, "decay_distance", TYPE_DOUBLE, true);
                this->log_data(impdata, "coefficients", TYPE_ARRAY, true);

                const auto a1 = impdata["center1"].GetArray();
                const auto a2 = impdata["center2"].GetArray();
                gp_Pnt center_1 = gp_Pnt(a1[0].GetDouble(), a1[1].GetDouble(), a1[2].GetDouble());
                gp_Pnt center_2 = gp_Pnt(a2[0].GetDouble(), a2[1].GetDouble(), a2[2].GetDouble());
                double r1 = impdata["r1"].GetDouble();
                double r2 = impdata["r2"].GetDouble();

                const double decay_distance = impdata["decay_distance"].GetDouble();
                const auto c = impdata["coefficients"].GetArray();
                std::vector<double> coeffs(c.Size());
                for(size_t i = 0; i < coeffs.size(); ++i){
                    coeffs[i] = c[i].GetDouble();
                }
                auto imp(material::Mandible::ImplantRegion(center_1, center_2, r1, r2, coeffs, decay_distance));
                material.emplace_back(new material::Mandible(name, o, i, path_points1, path_points2, C, has_implant, imp));
            } else {
                material.emplace_back(new material::Mandible(name, o, i, path_points1, path_points2, C, has_implant, material::Mandible::ImplantRegion()));
            }
        }
    }

    return material;
}

std::vector<std::unique_ptr<Geometry>> ProjectData::load_geometries(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    const auto& geometries = doc["geometry"].GetArray();
    std::vector<std::unique_ptr<Geometry>> geometry;
    size_t id = 0;
    for(const auto& geom:geometries){
        std::string absolute_path = folder_path;
        std::string geom_path = geom["file_path"].GetString();
        absolute_path.append(geom_path);

        this->log_data(geom, "do_topopt", TYPE_BOOL, true);
        this->log_data(geom, "material", TYPE_STRING, true);

        double scale = 1;
        if(this->log_data(geom, "scale", TYPE_DOUBLE, false)){
            scale = geom["scale"].GetDouble();
        }

        bool do_topopt = geom["do_topopt"].GetBool();

        size_t alt_size = 0;
        if(this->log_data(geom, "alt_materials", TYPE_ARRAY, false)){
            const auto& alt = geom["alt_materials"].GetArray();
            alt_size = alt.Size();
        }

        bool with_void = ((alt_size > 0) ? false : true) && do_topopt;
        if(this->log_data(geom, "with_void", TYPE_BOOL, false) && alt_size > 0 && do_topopt){
            with_void = geom["with_void"].GetBool();
        }
        if(this->do_simulation){
            do_topopt = false;
            with_void = false;
        }

        geometry.emplace_back(new Geometry(absolute_path, scale, this->type, this->topopt_element.get(), do_topopt, with_void, id));
        ++id;
    }
    return geometry;
}

void ProjectData::assign_materials(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    const auto& geometries = doc["geometry"].GetArray();
    std::vector<std::unique_ptr<Geometry>> geometry;
    size_t id = 0;
    for(const auto& geom:geometries){

        std::string mat_name(geom["material"].GetString());
        auto equal_name = [&mat_name](const std::unique_ptr<Material>& m)->bool{
            return mat_name == m->name;
        };
        auto it = std::find_if(this->materials.begin(), this->materials.end(), equal_name);
        logger::log_assert(it != this->materials.end(), logger::ERROR, "material with name '{}' not found", mat_name);
        std::vector<Material*> mats{it->get()};

        if(this->log_data(geom, "alt_materials", TYPE_ARRAY, false)){
            const auto& alt = geom["alt_materials"].GetArray();
            if(alt.Size() > 0){
                for(const auto& mat:alt){
                    logger::log_assert(mat.IsString(), logger::ERROR, "alt_materials must only contain the names of materials");
                    std::string mat_name(mat.GetString());
                    auto equal_name = [&mat_name](const std::unique_ptr<Material>& m)->bool{
                        return mat_name == std::string(m->name);
                    };
                    auto it = std::find_if(this->materials.begin(), this->materials.end(), equal_name);
                    logger::log_assert(it != this->materials.end(), logger::ERROR, "material with name '{}' not found", mat_name);
                    mats.push_back(it->get());
                }
            }
        }

        this->geometries[id]->set_materials(mats);

        ++id;
    }
}

std::unique_ptr<Pathfinding> ProjectData::load_pathfinder(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    using namespace pathfinding;

    std::unique_ptr<Pathfinding> pathfinder;

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
        pathfinder.reset(new MeshlessAStar(this->geometries[0]->shape, step, angle, choices, restriction, utils::PROBLEM_TYPE_2D));
    } else if(pathf["type"] == "visibility_graph"){
        this->log_data(pathf, "step", TYPE_DOUBLE, true);
        this->log_data(pathf, "max_turn_angle", TYPE_DOUBLE, true);
        double step = pathf["step"].GetDouble();
        double angle = pathf["max_turn_angle"].GetDouble();
        double restriction = 0;
        if(this->log_data(pathf, "restriction_size", TYPE_DOUBLE, false)){
            restriction = pathf["restriction_size"].GetDouble();
        }
        pathfinder.reset(new VisibilityGraph(this->geometries[0], step, angle, restriction, utils::PROBLEM_TYPE_2D));
    } else {
        logger::log_assert(false, logger::ERROR, "unknown pathfinding algorithm inserted: {}.", pathf["type"].GetString());
    }

    return pathfinder;
}

std::unique_ptr<Sizing> ProjectData::load_sizer(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    auto& sizing = doc["sizing"];
    std::unique_ptr<Sizing> sizer;
    this->log_data(sizing, "type", TYPE_STRING, true);
    if(sizing["type"] == "beam_sizing"){
        this->log_data(sizing, "element_type", TYPE_STRING, true);
        BeamElementFactory::BeamElementType t = BeamElementFactory::NONE;
        if(sizing["element_type"] == "beam_linear_2D"){
            t = BeamElementFactory::BEAM_LINEAR_2D;
        } else {
            logger::log_assert(false, logger::ERROR, "unknown element type for sizing algorithm: {}.", sizing["element_type"].GetString());
        }
        sizer.reset(new sizing::BeamSizing(this, t));
    } else if(sizing["type"] == "standard_sizing"){
        this->log_data(sizing, "element_size", TYPE_DOUBLE, true);
        double size = sizing["element_size"].GetDouble();
        double mult = 1.0;
        if(this->log_data(sizing, "oversizing", TYPE_DOUBLE, false)){
            mult = sizing["oversizing"].GetDouble();
        }
        sizer.reset(new sizing::StandardSizing(this, this->sizer_fea.get(), size, mult));
    }

    return sizer;
}

std::unique_ptr<FiniteElement> ProjectData::load_fea(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    auto& fea = doc["finite_element"];
    std::unique_ptr<FiniteElement> finite_element;
    //if(fea["type"] == "direct_solver"){
    //    finite_element.reset(new finite_element::DirectSolver());
    //} else if(fea["type"] == "gradient_descent"){
    //    this->log_data(fea, "eps", TYPE_DOUBLE, true);
    //    this->log_data(fea, "solver", TYPE_STRING, true);
    //    double eps = fea["eps"].GetDouble();
    //    std::string s = fea["solver"].GetString();
    //    finite_element::GradientDescent::Solver solver = finite_element::GradientDescent::Solver::STANDARD;
    //    if(s == "standard"){
    //        solver = finite_element::GradientDescent::Solver::STANDARD;
    //    } else if(s == "mma"){
    //        solver = finite_element::GradientDescent::Solver::MMA;
    //    } else if(s == "lagrange_mma"){
    //        solver = finite_element::GradientDescent::Solver::LAGRANGE_MMA;
    //    } else {
    //        logger::log_assert(false, logger::ERROR, "Unknown solver: {}", s);
    //    }
    //    finite_element.reset(new finite_element::GradientDescent(eps, solver));
    //} else if(fea["type"] == "PCG"){
    //    this->log_data(fea, "eps", TYPE_DOUBLE, true);
    //    this->log_data(fea, "preconditioner", TYPE_STRING, true);
    //    double eps = fea["eps"].GetDouble();
    //    std::string precond = fea["preconditioner"].GetString();
    //    finite_element::PCG::Preconditioner p = finite_element::PCG::Preconditioner::JACOBI;
    //    if(precond == "jacobi"){
    //        p = finite_element::PCG::Preconditioner::JACOBI;
    //    } else if(precond == "ssor"){
    //        p = finite_element::PCG::Preconditioner::SSOR;
    //    } else {
    //        logger::log_assert(false, logger::ERROR, "Unknown preconditioner: {}", precond);
    //    }
    //    finite_element.reset(new finite_element::PCG(eps, p));
    //} else if(fea["type"] == "mumps"){
    if(fea["type"] == "mumps"){
        finite_element.reset(new finite_element::MUMPSSolver(this->contact_data.contact_type, this->contact_data.rtol_abs, this->contact_data.max_step, this->contact_data.EPS_DISPL));
    } else if(fea["type"] == "eigen_pcg"){
        finite_element.reset(new finite_element::EigenPCG(this->contact_data.contact_type, this->contact_data.rtol_abs, this->contact_data.max_step, this->contact_data.EPS_DISPL));
    } else if(fea["type"] == "petsc_pcg"){
        this->log_data(fea, "backend", TYPE_STRING, true);
        std::string backend = fea["backend"].GetString();
        finite_element::PETScPCG::PETScBackend b = finite_element::PETScPCG::PETScBackend::CPU;
        if(backend == "cpu"){
            b = finite_element::PETScPCG::PETScBackend::CPU;
        } else if(backend == "cuda"){
            b = finite_element::PETScPCG::PETScBackend::CUDA;
        }

        finite_element.reset(new finite_element::PETScPCG(this->contact_data.contact_type, this->contact_data.rtol_abs, this->contact_data.max_step, this->contact_data.EPS_DISPL, b));
    }

    return finite_element;
}

std::unique_ptr<Meshing> ProjectData::load_mesher(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    auto& mesh = doc["mesher"];
    std::unique_ptr<Meshing> mesher;
    bool has_mesh_file_path = this->log_data(mesh, "file_path", TYPE_STRING, !this->do_meshing);
    if(has_mesh_file_path){
        std::string absolute_path = this->folder_path;
        std::string mesh_path = mesh["file_path"].GetString();
        absolute_path.append(mesh_path);
        this->mesh_file_internal = std::make_unique<meshing::MeshFile>(this->geometries,
                this->topopt_element.get(), this, this->thickness, 
                absolute_path, this->element_name);

        this->mesh_file = this->mesh_file_internal.get();
        // If meshing won't be done, use MeshFile instance as "mesher", so it
        // can load the raw mesh data and proceed from there.
        // As this is using unique_ptr, the mesh_file var will still
        // point to the correct place in memory.
        if(!this->do_meshing){
            return std::move(this->mesh_file_internal);
        }
    }
    if(mesh["type"] == "gmsh"){
        this->log_data(mesh, "element_size", TYPE_DOUBLE, true);
        size_t algorithm2D = 6;
        size_t algorithm3D = 4;
        double tmp_scale = 1;
        if(this->log_data(mesh, "algorithm2D", TYPE_INT, false)){
            algorithm2D = mesh["algorithm2D"].GetInt();
        }
        if(this->log_data(mesh, "algorithm3D", TYPE_INT, false)){
            algorithm3D = mesh["algorithm3D"].GetInt();
        }
        if(this->log_data(mesh, "tmp_scale", TYPE_DOUBLE, false)){
            tmp_scale = mesh["tmp_scale"].GetDouble();
        }
        double size = mesh["element_size"].GetDouble();
        mesher.reset(new meshing::Gmsh(this->geometries, this->topopt_element.get(), this, size, this->thickness, tmp_scale, algorithm2D, algorithm3D));
    }
    return mesher;
}

std::unique_ptr<DensityFilter> ProjectData::load_density_filter(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    this->log_data(doc, "density_filter", TYPE_OBJECT, true);

    auto& f = doc["density_filter"];
    std::unique_ptr<DensityFilter> filter;
    if(f["type"] == "convolution"){
        this->log_data(f, "radius", TYPE_DOUBLE, true);

        double radius = f["radius"].GetDouble();
        filter = std::make_unique<density_filter::Convolution>(radius);
    } else if(f["type"] == "helmholtz"){
        this->log_data(f, "radius", TYPE_DOUBLE, true);

        double radius = f["radius"].GetDouble();
        filter = std::make_unique<density_filter::Helmholtz>(radius);
    } else if(f["type"] == "averaging"){
        filter = std::make_unique<density_filter::Averaging>();
    }

    return filter;
}

std::unique_ptr<Projection> ProjectData::load_projection(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    if(!this->log_data(doc, "projection", TYPE_OBJECT, false)){
        return std::make_unique<projection::None>();
    }

    std::unique_ptr<Projection> filter;
    auto& f = doc["projection"];
    if(f["type"] == "none"){
        filter = std::make_unique<projection::None>();
    } else if(f["type"] == "threshold"){
        this->log_data(f, "beta", TYPE_OBJECT, false);
        this->log_data(f, "eta", TYPE_DOUBLE, false);

        auto beta = this->get_projection_parameter(f["beta"]);
        double eta = f["eta"].GetDouble();

        filter = std::make_unique<projection::Threshold>(beta, eta);
    } else if(f["type"] == "heaviside"){
        this->log_data(f, "beta", TYPE_OBJECT, false);
        this->log_data(f, "eta", TYPE_DOUBLE, false);

        auto beta = this->get_projection_parameter(f["beta"]);
        double eta = f["eta"].GetDouble();

        filter = std::make_unique<projection::Heaviside>(beta, eta);
    }

    return filter;
}

std::unique_ptr<DensityBasedOptimizer> ProjectData::load_topopt_optimizer(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    auto& to = doc["topopt"];
    this->log_data(to, "type", TYPE_STRING, true);
    if(to["type"] == "mma"){
        this->log_data(to, "rho_init", TYPE_DOUBLE, true);
        this->log_data(to, "xtol_abs", TYPE_DOUBLE, true);
        this->log_data(to, "ftol_rel", TYPE_DOUBLE, true);
        this->log_data(to, "result_threshold", TYPE_DOUBLE, true);
        this->log_data(to, "save_result", TYPE_BOOL, true);
        this->log_data(to, "pc", TYPE_INT, true);
        this->log_data(to, "psi", TYPE_DOUBLE, true);
        this->log_data(to, "objective", TYPE_ARRAY, true);
        this->log_data(to, "constraints", TYPE_ARRAY, true);
        this->log_data(to, "asyminit", TYPE_DOUBLE, true);
        this->log_data(to, "asymdec", TYPE_DOUBLE, true);
        this->log_data(to, "asyminc", TYPE_DOUBLE, true);
        this->log_data(to, "minfac", TYPE_DOUBLE, true);
        this->log_data(to, "maxfac", TYPE_DOUBLE, true);
        this->log_data(to, "c", TYPE_DOUBLE, true);

        this->density_filter = this->load_density_filter(to);
        this->projection = this->load_projection(to);

        double rho_init = to["rho_init"].GetDouble();
        double xtol_abs = to["xtol_abs"].GetDouble();
        double ftol_rel = to["ftol_rel"].GetDouble();
        double result_threshold = to["result_threshold"].GetDouble();
        bool save_result = to["save_result"].GetBool();
        int pc = to["pc"].GetInt();
        double psi = to["psi"].GetDouble();

        double asyminit = to["asyminit"].GetDouble();
        double asymdec = to["asymdec"].GetDouble();
        double asyminc = to["asyminc"].GetDouble();
        double minfac = to["minfac"].GetDouble();
        double maxfac = to["maxfac"].GetDouble();
        double c = to["c"].GetDouble();

        std::vector<std::unique_ptr<DensityBasedFunction>> objective;
        std::vector<double> weights;
        this->get_topopt_objective_functions(to["objective"], pc, psi, objective, weights);

        std::vector<DensityBasedConstraint> constraints;
        this->get_topopt_constraints(to["constraints"], pc, psi, constraints);

        return std::make_unique<optimizer::density_based::MMA>(this->density_filter.get(), this->projection.get(), this, std::move(objective), std::move(weights), std::move(constraints), asyminit, asymdec, asyminc, minfac, maxfac, c, pc, psi, rho_init, xtol_abs, ftol_rel, result_threshold, save_result);
    } else if(to["type"] == "newton"){
        this->log_data(to, "rho_init", TYPE_DOUBLE, true);
        this->log_data(to, "xtol_abs", TYPE_DOUBLE, true);
        this->log_data(to, "ftol_rel", TYPE_DOUBLE, true);
        this->log_data(to, "result_threshold", TYPE_DOUBLE, true);
        this->log_data(to, "save_result", TYPE_BOOL, true);
        this->log_data(to, "pc", TYPE_INT, true);
        this->log_data(to, "psi", TYPE_DOUBLE, true);
        this->log_data(to, "objective", TYPE_ARRAY, true);
        this->log_data(to, "constraints", TYPE_ARRAY, true);

        this->density_filter = this->load_density_filter(to);
        this->projection = this->load_projection(to);

        double rho_init = to["rho_init"].GetDouble();
        double xtol_abs = to["xtol_abs"].GetDouble();
        double ftol_rel = to["ftol_rel"].GetDouble();
        double result_threshold = to["result_threshold"].GetDouble();
        bool save_result = to["save_result"].GetBool();
        int pc = to["pc"].GetInt();
        double psi = to["psi"].GetDouble();

        std::vector<std::unique_ptr<DensityBasedFunction>> objective;
        std::vector<double> weights;
        this->get_topopt_objective_functions(to["objective"], pc, psi, objective, weights);

        std::vector<DensityBasedConstraint> constraints;
        this->get_topopt_constraints(to["constraints"], pc, psi, constraints);

        return std::make_unique<optimizer::density_based::Newton>(this->density_filter.get(), this->projection.get(), this, std::move(objective), std::move(weights), std::move(constraints), pc, psi, rho_init, xtol_abs, ftol_rel, result_threshold, save_result);
    }
    logger::log_assert(false, logger::ERROR, "optimization method \"{}\" not found.", type);

    return nullptr;
}

std::unique_ptr<NodeShapeBasedOptimizer> ProjectData::load_shopt_optimizer(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    auto& to = doc["shape_opt"];
    this->log_data(to, "type", TYPE_STRING, true);
    if(to["type"] == "mma"){
        this->log_data(to, "xtol_abs", TYPE_DOUBLE, true);
        this->log_data(to, "ftol_rel", TYPE_DOUBLE, true);
        this->log_data(to, "save_result", TYPE_BOOL, true);
        this->log_data(to, "objective", TYPE_ARRAY, true);
        this->log_data(to, "constraints", TYPE_ARRAY, true);
        this->log_data(to, "asyminit", TYPE_DOUBLE, true);
        this->log_data(to, "asymdec", TYPE_DOUBLE, true);
        this->log_data(to, "asyminc", TYPE_DOUBLE, true);
        this->log_data(to, "minfac", TYPE_DOUBLE, true);
        this->log_data(to, "maxfac", TYPE_DOUBLE, true);
        this->log_data(to, "c", TYPE_DOUBLE, true);
        this->log_data(to, "shape_set", TYPE_OBJECT, true);

        double xtol_abs = to["xtol_abs"].GetDouble();
        double ftol_rel = to["ftol_rel"].GetDouble();
        bool save_result = to["save_result"].GetBool();

        double asyminit = to["asyminit"].GetDouble();
        double asymdec = to["asymdec"].GetDouble();
        double asyminc = to["asyminc"].GetDouble();
        double minfac = to["minfac"].GetDouble();
        double maxfac = to["maxfac"].GetDouble();
        double c = to["c"].GetDouble();

        std::vector<std::unique_ptr<NodeShapeBasedFunction>> objective;
        std::vector<double> weights;
        this->get_shopt_objective_functions(to["objective"], objective, weights);

        std::vector<NodeShapeBasedConstraint> constraints;
        this->get_shopt_constraints(to["constraints"], constraints);

        std::vector<Geometry*> geoms;
        geoms.reserve(this->geometries.size());
        for(auto& g:this->geometries){
            geoms.push_back(g.get());
        }

        auto shape_set = this->get_shape_operations(to["shape_set"]);

        ShapeHandler sh(this->topopt_mesher.get(), geoms, std::move(shape_set));

        return std::make_unique<optimizer::node_shape_based::MMA>(std::move(sh), this, std::move(objective), std::move(weights), std::move(constraints), asyminit, asymdec, asyminc, minfac, maxfac, c, xtol_abs, ftol_rel, save_result);
    }
    logger::log_assert(false, logger::ERROR, "optimization method \"{}\" not found.", type);

    return nullptr;
}

std::unique_ptr<DensityBasedFunction> ProjectData::get_topopt_function(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc, double pc, double psiK){
    this->log_data(doc, "type", TYPE_STRING, true);
    std::string type = doc["type"].GetString();
    if(type == "compliance"){
        return std::make_unique<function::density_based::Compliance>(this->topopt_mesher.get(), pc, psiK);
    } else if(type == "volume"){
        return std::make_unique<function::density_based::Volume>(this->topopt_mesher.get());
    } else if(type == "mass"){
        return std::make_unique<function::density_based::Mass>(this->topopt_mesher.get());
    } else if(type == "mass_first_material"){
        return std::make_unique<function::density_based::MassFirstMaterial>(this->topopt_mesher.get());
    } else if(type == "global_stress_pnorm_normalized"){
        this->log_data(doc, "P", TYPE_DOUBLE, true);
        this->log_data(doc, "pt", TYPE_DOUBLE, true);
        this->log_data(doc, "psi", TYPE_DOUBLE, true);
        double P = doc["P"].GetDouble();
        double pt = doc["pt"].GetDouble();
        double psiS = doc["psi"].GetDouble();
        return std::make_unique<function::density_based::GlobalStressPnormNormalized>(this->topopt_mesher.get(), this->topopt_fea.get(), pc, P, pt, psiK, -psiS);
    } else if(type == "global_stress_pnorm"){
        this->log_data(doc, "P", TYPE_DOUBLE, true);
        this->log_data(doc, "pt", TYPE_DOUBLE, true);
        this->log_data(doc, "psi", TYPE_DOUBLE, true);
        double P = doc["P"].GetDouble();
        double pt = doc["pt"].GetDouble();
        double psiS = doc["psi"].GetDouble();
        return std::make_unique<function::density_based::GlobalStressPnorm>(this->topopt_mesher.get(), this->topopt_fea.get(), pc, P, pt, psiK, -psiS);
    } else if(type == "global_stress_heaviside"){
        this->log_data(doc, "max_stress", TYPE_DOUBLE, true);
        this->log_data(doc, "C", TYPE_DOUBLE, true);
        this->log_data(doc, "pt", TYPE_DOUBLE, true);
        this->log_data(doc, "psi", TYPE_DOUBLE, true);
        double C = doc["C"].GetDouble();
        double max_stress = doc["max_stress"].GetDouble();
        double pt = doc["pt"].GetDouble();
        double psiS = doc["psi"].GetDouble();
        return std::make_unique<function::density_based::GlobalStressHeaviside>(this->topopt_mesher.get(), this->topopt_fea.get(), max_stress, C, pc, pt, psiK, -psiS);
    } else if(type == "omni_machining"){
        this->log_data(doc, "beta1", TYPE_DOUBLE, true);
        this->log_data(doc, "beta2", TYPE_DOUBLE, true);
        this->log_data(doc, "v", TYPE_DOUBLE, true);
        this->log_data(doc, "center", TYPE_ARRAY, true);
        this->log_data(doc, "axis", TYPE_ARRAY, true);
        this->log_data(doc, "L", TYPE_DOUBLE, true);
        double L = doc["L"].GetDouble();
        double beta1 = doc["beta1"].GetDouble();
        double beta2 = doc["beta2"].GetDouble();
        double v = doc["v"].GetDouble();

        auto ca = doc["center"].GetArray();
        gp_Pnt c(ca[0].GetDouble(), ca[1].GetDouble(), ca[2].GetDouble());
        auto aa = doc["axis"].GetArray();
        gp_Dir a(aa[0].GetDouble(), aa[1].GetDouble(), aa[2].GetDouble());

        return std::make_unique<function::density_based::OmniMachining>(this->topopt_mesher.get(), this->density_filter.get(), c, a, v, beta1, beta2, L);
    } else if(type == "am_support"){
        this->log_data(doc, "beta", TYPE_DOUBLE, true);
        this->log_data(doc, "L", TYPE_DOUBLE, true);
        this->log_data(doc, "v", TYPE_DOUBLE, true);
        this->log_data(doc, "support_angle", TYPE_DOUBLE, true);
        this->log_data(doc, "axis", TYPE_ARRAY, true);
        double beta = doc["beta"].GetDouble();
        double L = doc["L"].GetDouble();
        double v = doc["v"].GetDouble();
        double angle = doc["support_angle"].GetDouble();

        auto aa = doc["axis"].GetArray();
        gp_Dir a(aa[0].GetDouble(), aa[1].GetDouble(), aa[2].GetDouble());

        return std::make_unique<function::density_based::AMSupport>(this->topopt_mesher.get(), this->density_filter.get(), this->projection.get(), a, v, L, beta, angle*M_PI/180);
    } else if(type == "mechanostat"){
        this->log_data(doc, "beta", TYPE_DOUBLE, true);
        this->log_data(doc, "traction", TYPE_ARRAY, true);
        this->log_data(doc, "compression", TYPE_ARRAY, true);
        this->log_data(doc, "shear", TYPE_ARRAY, true);
        double beta = doc["beta"].GetDouble();
        auto traction = doc["traction"].GetArray();
        auto compression = doc["compression"].GetArray();
        auto shear = doc["shear"].GetArray();
        logger::log_assert(traction.Size() == 2, logger::ERROR, "\"traction\" item must have size 2.");
        logger::log_assert(compression.Size() == 2, logger::ERROR, "\"compression\" item must have size 2.");
        logger::log_assert(shear.Size() == 2, logger::ERROR, "\"shear\" item must have size 2.");
        function::density_based::Mechanostat::Range t{traction[0].GetDouble(), traction[1].GetDouble()};
        function::density_based::Mechanostat::Range c{compression[0].GetDouble(), compression[1].GetDouble()};
        function::density_based::Mechanostat::Range s{shear[0].GetDouble(), shear[1].GetDouble()};
        return std::make_unique<function::density_based::Mechanostat>(this->topopt_mesher.get(), this->topopt_fea.get(), pc, psiK, beta, t, c, s, this->type);
    }
    logger::log_assert(false, logger::ERROR, "function \"{}\" not found.", type);

    return nullptr;
}

std::unique_ptr<NodeShapeBasedFunction> ProjectData::get_shopt_function(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    this->log_data(doc, "type", TYPE_STRING, true);
    std::string type = doc["type"].GetString();

    if(type == "compliance"){
        return std::make_unique<function::node_shape_based::Compliance>(this->topopt_mesher.get());
    } else if(type == "volume"){
        return std::make_unique<function::node_shape_based::Volume>(this->topopt_mesher.get());
    } else if(type == "global_stress_heaviside"){
        this->log_data(doc, "max_stress", TYPE_DOUBLE, true);
        this->log_data(doc, "C", TYPE_DOUBLE, true);
        double C = doc["C"].GetDouble();
        double max_stress = doc["max_stress"].GetDouble();
        return std::make_unique<function::node_shape_based::GlobalStressHeaviside>(this->topopt_mesher.get(), this->topopt_fea.get(), max_stress, C);
    } else if(type == "mechanostat"){
        this->log_data(doc, "beta", TYPE_DOUBLE, true);
        this->log_data(doc, "traction", TYPE_ARRAY, true);
        this->log_data(doc, "compression", TYPE_ARRAY, true);
        this->log_data(doc, "shear", TYPE_ARRAY, true);
        double beta = doc["beta"].GetDouble();
        auto traction = doc["traction"].GetArray();
        auto compression = doc["compression"].GetArray();
        auto shear = doc["shear"].GetArray();
        logger::log_assert(traction.Size() == 2, logger::ERROR, "\"traction\" item must have size 2.");
        logger::log_assert(compression.Size() == 2, logger::ERROR, "\"compression\" item must have size 2.");
        logger::log_assert(shear.Size() == 2, logger::ERROR, "\"shear\" item must have size 2.");
        function::density_based::Mechanostat::Range t{traction[0].GetDouble(), traction[1].GetDouble()};
        function::density_based::Mechanostat::Range c{compression[0].GetDouble(), compression[1].GetDouble()};
        function::density_based::Mechanostat::Range s{shear[0].GetDouble(), shear[1].GetDouble()};
        return std::make_unique<function::node_shape_based::Mechanostat>(this->topopt_mesher.get(), this->topopt_fea.get(), beta, t, c, s, this->type);
    }

    logger::log_assert(false, logger::ERROR, "function \"{}\" not found.", type);

    return nullptr;
}

void ProjectData::get_topopt_constraints(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc, double pc, double psi, std::vector<DensityBasedConstraint>& functions){
    auto array = doc.GetArray();
    for(auto& f:array){
        std::vector<Constraint::Type> types;
        std::vector<double> bounds;

        this->log_data(f, "type", TYPE_STRING, true);
        std::string type = f["type"].GetString();
        if(f.HasMember("less_than")){
            types.push_back(Constraint::Type::LESS_THAN);
            bounds.push_back(f["less_than"].GetDouble());
        }
        if(f.HasMember("greater_than")){
            types.push_back(Constraint::Type::GREATER_THAN);
            bounds.push_back(f["greater_than"].GetDouble());
        }
        if(f.HasMember("equals")){
            types.push_back(Constraint::Type::EQUAL);
            bounds.push_back(f["equals"].GetDouble());
        }
        logger::log_assert(types.size() > 0, logger::ERROR, "constraint {} is missing at least one of the following: \"equals\", \"greater_than\", \"less_than\"", f["type"].GetString());

        functions.emplace_back(this->get_topopt_function(f, pc, psi), types, bounds);
    }
}

void ProjectData::get_topopt_objective_functions(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc, double pc, double psi, std::vector<std::unique_ptr<DensityBasedFunction>>& functions, std::vector<double>& weights){
    auto array = doc.GetArray();
    for(auto& f:array){
        functions.push_back(this->get_topopt_function(f, pc, psi));
        weights.push_back(f["weight"].GetDouble());
    }
}

void ProjectData::get_shopt_constraints(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc, std::vector<NodeShapeBasedConstraint>& functions){
    auto array = doc.GetArray();
    for(auto& f:array){
        std::vector<Constraint::Type> types;
        std::vector<double> bounds;

        this->log_data(f, "type", TYPE_STRING, true);
        std::string type = f["type"].GetString();
        if(f.HasMember("less_than")){
            types.push_back(Constraint::Type::LESS_THAN);
            bounds.push_back(f["less_than"].GetDouble());
        }
        if(f.HasMember("greater_than")){
            types.push_back(Constraint::Type::GREATER_THAN);
            bounds.push_back(f["greater_than"].GetDouble());
        }
        if(f.HasMember("equals")){
            types.push_back(Constraint::Type::EQUAL);
            bounds.push_back(f["equals"].GetDouble());
        }
        logger::log_assert(types.size() > 0, logger::ERROR, "constraint {} is missing at least one of the following: \"equals\", \"greater_than\", \"less_than\"", f["type"].GetString());

        functions.emplace_back(this->get_shopt_function(f), types, bounds);
    }
}

void ProjectData::get_shopt_objective_functions(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc, std::vector<std::unique_ptr<NodeShapeBasedFunction>>& functions, std::vector<double>& weights){
    auto array = doc.GetArray();
    for(auto& f:array){
        functions.push_back(this->get_shopt_function(f));
        weights.push_back(f["weight"].GetDouble());
    }
}

Projection::Parameter ProjectData::get_projection_parameter(const rapidjson::GenericValue<rapidjson::UTF8<>>& p) const{
    this->log_data(p, "initial", TYPE_DOUBLE, true);
    this->log_data(p, "final", TYPE_DOUBLE, true);
    this->log_data(p, "value_step", TYPE_DOUBLE, true);
    this->log_data(p, "iteration_step", TYPE_INT, true);

    Projection::Parameter param;
    param.value = p["initial"].GetDouble();
    param.final_value = p["final"].GetDouble();
    param.value_step = p["value_step"].GetDouble();
    param.iteration_step = p["iteration_step"].GetInt();

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

std::vector<std::unique_ptr<Field>> ProjectData::load_fields(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    std::vector<std::unique_ptr<Field>> fields;
    const auto& fields_array = doc["fields"].GetArray();
    for(auto& f : fields_array){
        logger::log_assert(f.IsObject(), logger::ERROR, "Each field must be stored as a JSON object");
        this->log_data(f, "type", TYPE_STRING, true);
        std::string type = f["type"].GetString();
        if(type == "orthotropic_flow"){
            this->log_data(f, "geometries", TYPE_ARRAY, true);
            this->log_data(f, "boundary_conditions", TYPE_ARRAY, true);
            this->log_data(f, "display", TYPE_BOOL, true);
            this->log_data(f, "alpha", TYPE_DOUBLE, true);

            std::vector<Geometry*> geoms;
            std::vector<CrossSection> cs;
            std::vector<double> coeffs;
            bool show = f["display"].GetBool();
            const double alpha = f["alpha"].GetDouble();
            for(const auto& g:f["geometries"].GetArray()){
                logger::log_assert((size_t)g.GetInt() < this->geometries.size(), logger::ERROR, "there are only {} geometries being loaded, but geometry number {} was selected for field of type {}", this->geometries.size(), g.GetInt(), type);
                geoms.push_back(this->geometries[g.GetInt()].get());
            }
            for(const auto& c:f["boundary_conditions"].GetArray()){
                logger::log_assert(c.IsObject(), logger::ERROR, "Each boundary_condition element must be stored as a JSON object");
                cs.push_back(get_cross_section(c));
                this->log_data(c, "coeff", TYPE_DOUBLE, true);
                coeffs.push_back(c["coeff"].GetDouble());
            }

            fields.emplace_back(std::make_unique<field::OrthotropicFlow>(this->topopt_element.get(), geoms, cs, coeffs, alpha, this->thickness, show));
        } else if(type == "principal_stress"){
            this->log_data(f, "display", TYPE_BOOL, true);
            this->log_data(f, "initial_material", TYPE_STRING, true);
            this->log_data(f, "max_it", TYPE_INT, true);

            bool show = f["display"].GetBool();
            std::string mat_name = f["initial_material"].GetString();
            size_t max_it = f["max_it"].GetInt64();

            Material* init_mat = nullptr;
            for(auto& mat:this->materials){
                if(mat->name == mat_name){
                    init_mat = mat.get();
                    break;
                }
            }

            fields.emplace_back(std::make_unique<field::PrincipalStress>(this->topopt_element.get(), this, init_mat, max_it, this->thickness, show));
        }
    }

    return fields;
}

std::vector<Force> ProjectData::get_loads(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    std::vector<Force> forces;
    if(this->type == utils::PROBLEM_TYPE_2D){
        for(auto& f : doc.GetArray()){
            logger::log_assert(f.IsObject(), logger::ERROR, "Each load must be stored as a JSON object");
            this->log_data(f, "load", TYPE_ARRAY, true);
            auto loads = f["load"].GetArray();
            logger::log_assert(loads.Size() == 2, logger::ERROR, "Load vector must have exactly two dimensions in 2D problems");

            gp_Vec l(loads[0].GetDouble(), loads[1].GetDouble(), 0);

            auto S = this->get_cross_section(f);
            forces.emplace_back(std::move(S), l);
        }
    } else if(this->type == utils::PROBLEM_TYPE_3D) {
        for(auto& f : doc.GetArray()){
            logger::log_assert(f.IsObject(), logger::ERROR, "Each load must be stored as a JSON object");
            this->log_data(f, "load", TYPE_ARRAY, true);
            auto loads = f["load"].GetArray();
            logger::log_assert(loads.Size() == 3, logger::ERROR, "Load vector must have exactly three dimensions in 3D problems");

            gp_Vec l(loads[0].GetDouble(), loads[1].GetDouble(), loads[2].GetDouble());

            auto S = this->get_cross_section(f);
            forces.emplace_back(std::move(S), l);
        }
    }
    return forces;
}

std::vector<Support> ProjectData::get_support(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    std::vector<Support> supports;
    if(this->type == utils::PROBLEM_TYPE_2D){
        for(auto& f : doc.GetArray()){
            logger::log_assert(f.IsObject(), logger::ERROR, "Each support must be stored as a JSON object");
            this->log_data(f, "X", TYPE_BOOL, true);
            this->log_data(f, "Y", TYPE_BOOL, true);
            this->log_data(f, "MZ", TYPE_BOOL, true);
            bool X = f["X"].GetBool();
            bool Y = f["Y"].GetBool();
            bool MZ = f["MZ"].GetBool();

            auto S = this->get_cross_section(f);
            supports.emplace_back(X, Y, MZ, S);
        }
    } else if(this->type == utils::PROBLEM_TYPE_3D) {
        for(auto& f : doc.GetArray()){
            logger::log_assert(f.IsObject(), logger::ERROR, "Each support must be stored as a JSON object");
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

            auto S = this->get_cross_section(f);
            supports.emplace_back(X, Y, Z, MX, MY, MZ, S);
        }
    }
    return supports;
}

std::vector<Spring> ProjectData::get_springs(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    std::vector<Spring> springs;
    for(auto& f : doc.GetArray()){
        logger::log_assert(f.IsObject(), logger::ERROR, "Each load must be stored as a JSON object");
        this->log_data(f, "L", TYPE_ARRAY, true);
        this->log_data(f, "normal", TYPE_ARRAY, true);
        this->log_data(f, "v", TYPE_ARRAY, true);
        if(this->type == utils::PROBLEM_TYPE_3D) {
            this->log_data(f, "w", TYPE_ARRAY, true);
        }
        this->log_data(f, "material", TYPE_STRING, true);
        auto L = f["L"].GetArray();
        auto normal = f["normal"].GetArray();
        auto v = f["v"].GetArray();

        std::array<double, 3> l{L[0].GetDouble(), L[1].GetDouble(), 0};
        gp_Dir nv(normal[0].GetDouble(), normal[1].GetDouble(), 0);
        gp_Dir vv(v[0].GetDouble(), v[1].GetDouble(), 0);
        gp_Dir wv(0, 0, 1);

        if(this->type == utils::PROBLEM_TYPE_2D){
            logger::log_assert(L.Size() == 2, logger::ERROR, "Length vector must have exactly two dimensions in 2D problems");
            logger::log_assert(normal.Size() == 2, logger::ERROR, "Normal vector must have exactly two dimensions in 2D problems");
            logger::log_assert(v.Size() == 2, logger::ERROR, "'v' vector must have exactly two dimensions in 2D problems");
        } else if(this->type == utils::PROBLEM_TYPE_3D) {
            logger::log_assert(L.Size() == 3, logger::ERROR, "Length vector must have exactly three dimensions in 3D problems");
            logger::log_assert(normal.Size() == 3, logger::ERROR, "Normal vector must have exactly three dimensions in 3D problems");
            logger::log_assert(v.Size() == 3, logger::ERROR, "'v' vector must have exactly three dimensions in 3D problems");

            auto w = f["w"].GetArray();
            logger::log_assert(w.Size() == 3, logger::ERROR, "'w' vector must have exactly two dimensions in 3D problems");

            l[2] = L[2].GetDouble();
            nv.SetZ(normal[2].GetDouble());
            vv.SetZ(v[2].GetDouble());

            wv = gp_Dir(w[0].GetDouble(), w[1].GetDouble(), w[2].GetDouble());
        }

        std::string mat_name(f["material"].GetString());
        auto equal_name = [&mat_name](const std::unique_ptr<Material>& m)->bool{
            return mat_name == m->name;
        };
        auto it = std::find_if(this->materials.begin(), this->materials.end(), equal_name);
        logger::log_assert(it != this->materials.end(), logger::ERROR, "material with name '{}' not found", mat_name);

        Material* mat(it->get());

        auto S = this->get_cross_section(f);
        springs.emplace_back(S, this->thickness, nv, vv, wv, mat, l, this->topopt_element.get(), this->type);
    }
    return springs;
}

std::vector<InternalLoads> ProjectData::get_internal_loads(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    std::vector<InternalLoads> internal_loads;
    for(auto& f : doc.GetArray()){
        logger::log_assert(f.IsObject(), logger::ERROR, "Each load must be stored as a JSON object");
        this->log_data(f, "F", TYPE_ARRAY, true);
        this->log_data(f, "M", TYPE_ARRAY, true);
        this->log_data(f, "normal", TYPE_ARRAY, true);
        this->log_data(f, "v", TYPE_ARRAY, true);
        if(this->type == utils::PROBLEM_TYPE_3D) {
            this->log_data(f, "w", TYPE_ARRAY, true);
        }
        this->log_data(f, "material", TYPE_STRING, true);
        auto normal = f["normal"].GetArray();
        auto v = f["v"].GetArray();
        auto Fa = f["F"].GetArray();
        auto Ma = f["M"].GetArray();

        std::array<double, 3> F{Fa[0].GetDouble(), Fa[1].GetDouble(), 0};
        std::array<double, 3> M{Ma[0].GetDouble(), Ma[1].GetDouble(), 0};
        gp_Dir nv(normal[0].GetDouble(), normal[1].GetDouble(), 0);
        gp_Dir vv(v[0].GetDouble(), v[1].GetDouble(), 0);
        gp_Dir wv(0, 0, 1);

        if(this->type == utils::PROBLEM_TYPE_2D){
            logger::log_assert(normal.Size() == 2, logger::ERROR, "Normal vector must have exactly two dimensions in 2D problems");
            logger::log_assert(v.Size() == 2, logger::ERROR, "'v' vector must have exactly two dimensions in 2D problems");
            logger::log_assert(Fa.Size() == 2, logger::ERROR, "Force vector must have exactly two dimensions in 2D problems");
            logger::log_assert(Ma.Size() == 2, logger::ERROR, "Moment vector must have exactly two dimensions in 2D problems");
        } else if(this->type == utils::PROBLEM_TYPE_3D) {
            logger::log_assert(normal.Size() == 3, logger::ERROR, "Normal vector must have exactly three dimensions in 3D problems");
            logger::log_assert(v.Size() == 3, logger::ERROR, "'v' vector must have exactly three dimensions in 3D problems");
            logger::log_assert(Fa.Size() == 3, logger::ERROR, "Force vector must have exactly three dimensions in 3D problems");
            logger::log_assert(Ma.Size() == 3, logger::ERROR, "Moment vector must have exactly three dimensions in 3D problems");

            auto w = f["w"].GetArray();
            logger::log_assert(w.Size() == 3, logger::ERROR, "'w' vector must have exactly two dimensions in 3D problems");

            nv.SetZ(normal[2].GetDouble());
            vv.SetZ(v[2].GetDouble());
            F[2] = Fa[2].GetDouble();
            M[2] = Ma[2].GetDouble();

            wv = gp_Dir(w[0].GetDouble(), w[1].GetDouble(), w[2].GetDouble());
        }

        std::string mat_name(f["material"].GetString());
        auto equal_name = [&mat_name](const std::unique_ptr<Material>& m)->bool{
            return mat_name == m->name;
        };
        auto it = std::find_if(this->materials.begin(), this->materials.end(), equal_name);
        logger::log_assert(it != this->materials.end(), logger::ERROR, "material with name '{}' not found", mat_name);

        Material* mat(it->get());

        auto S = this->get_cross_section(f);
        internal_loads.emplace_back(S, this->thickness, nv, vv, wv, mat, F, M, this->topopt_element.get(), this->topopt_boundary_element.get(), this->type);
    }
    return internal_loads;
}

CrossSection ProjectData::get_cross_section(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc) const{
    if(this->type == utils::PROBLEM_TYPE_2D){
        bool has_vertices = this->log_data(doc, "vertices", TYPE_ARRAY, false);
        bool has_point = this->log_data(doc, "point", TYPE_ARRAY, false);

        logger::log_assert(has_vertices || has_point, logger::ERROR, "invalid cross-section definition in configuration file.");
        if(has_vertices){
            auto vertices = doc["vertices"].GetArray();

            std::vector<gp_Pnt> vlist;
            for(auto& v : vertices){
                logger::log_assert(v.Size() == 2, logger::ERROR, "Vertices must have exactly two dimensions in 2D problems");
                vlist.emplace_back(v[0].GetDouble(), v[1].GetDouble(), 0);
            }
            return CrossSection(vlist, this->thickness);
        } else if(has_point){
            const auto& pa = doc["point"].GetArray();
            logger::log_assert(pa.Size() == 2, logger::ERROR, "Vertices must have exactly two dimensions in 2D problems");
            gp_Pnt p(pa[0].GetDouble(), pa[1].GetDouble(), 0);

            return CrossSection(p);
        }
    } else if(this->type == utils::PROBLEM_TYPE_3D) {
        bool has_rect = this->log_data(doc, "rectangle", TYPE_OBJECT, false);
        bool has_file = this->log_data(doc, "file", TYPE_OBJECT, false);
        bool has_point = this->log_data(doc, "point", TYPE_ARRAY, false);

        logger::log_assert(has_rect || has_file || has_point, logger::ERROR, "invalid cross-section definition in configuration file.");
        if(has_rect){
            const auto& rect = doc["rectangle"];
            this->log_data(rect, "center", TYPE_ARRAY, true);
            this->log_data(rect, "normal", TYPE_ARRAY, true);
            this->log_data(rect, "w", TYPE_DOUBLE, true);
            this->log_data(rect, "h", TYPE_DOUBLE, true);
            this->log_data(rect, "rotation", TYPE_DOUBLE, true);
            auto c = rect["center"].GetArray();
            auto n = rect["normal"].GetArray();
            CrossSection::Rectangle r{
                rect["w"].GetDouble(),
                rect["h"].GetDouble(),
                gp_Pnt(
                    c[0].GetDouble(), 
                    c[1].GetDouble(),
                    c[2].GetDouble()
                ),
                gp_Dir(
                    n[0].GetDouble(), 
                    n[1].GetDouble(),
                    n[2].GetDouble()
                ),
                rect["rotation"].GetDouble()
            };
            return CrossSection(r);
        } else if(has_file) {
            logger::log_assert(!this->generate_beams,
                               logger::ERROR, "cross-section of type file is currently not fully compatible with beam generation");
            const auto& file = doc["file"];
            this->log_data(file, "path", TYPE_STRING, true);

            double scale = 1;
            if(this->log_data(file, "scale", TYPE_DOUBLE, false)){
                scale = file["scale"].GetDouble();
            }

            std::string absolute_path = folder_path;
            std::string path = file["path"].GetString();
            absolute_path.append(path);

            return CrossSection(absolute_path, scale);
        } else if(has_point){
            const auto& pa = doc["point"].GetArray();
            logger::log_assert(pa.Size() == 3, logger::ERROR, "Vertices must have exactly three dimensions in 3D problems");
            gp_Pnt p(pa[0].GetDouble(), pa[1].GetDouble(), pa[2].GetDouble());

            return CrossSection(p);
        }
    }
    logger::log_assert(false, logger::ERROR, "unknown problem type detected in get_cross_section(), this shouldn't have happened.");
    return CrossSection(gp_Pnt(0,0,0));
}


std::vector<SubProblem> ProjectData::load_sub_problems(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    const auto& sub_problem_array = doc["sub_problems"].GetArray();
    std::vector<SubProblem> sub_problems;
    for(const auto& sp:sub_problem_array){
        logger::log_assert(sp.IsObject(), logger::ERROR, "\"sub_problems\" must be an array of objects");
        SubProblem curr_sp;
        if(this->log_data(sp, "loads", TYPE_ARRAY, false)){
            for(auto& f:sp["loads"].GetArray()){
                curr_sp.forces.push_back(&this->forces[f.GetInt()]);
            }
        }
        if(this->log_data(sp, "supports", TYPE_ARRAY, false)){
            for(auto& f:sp["supports"].GetArray()){
                curr_sp.supports.push_back(&this->supports[f.GetInt()]);
            }
        }
        if(this->log_data(sp, "springs", TYPE_ARRAY, false)){
            for(auto& f:sp["springs"].GetArray()){
                curr_sp.springs.push_back(&this->springs[f.GetInt()]);
            }
        }
        if(this->log_data(sp, "internal_loads", TYPE_ARRAY, false)){
            for(auto& f:sp["internal_loads"].GetArray()){
                curr_sp.internal_loads.push_back(&this->internal_loads[f.GetInt()]);
            }
        }
        sub_problems.push_back(std::move(curr_sp));
    }

    return sub_problems;
}

std::unique_ptr<shape_op::ShapeOp> ProjectData::get_shape_operations(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc) const{
    this->log_data(doc, "type", TYPE_STRING, true);
    const std::string type = doc["type"].GetString();
    if(type == "geometry"){
        this->log_data(doc, "id", TYPE_INT, true);
        const size_t id = doc["id"].GetUint64();
        logger::log_assert(id < this->geometries.size(), logger::ERROR,
                "unknown geometry id: {}", id);
        return std::make_unique<shape_op::Geometry>(id);
    } else if(type == "union"){
        this->log_data(doc, "shape1", TYPE_OBJECT, true);
        this->log_data(doc, "shape2", TYPE_OBJECT, true);

        return std::make_unique<shape_op::Union>(
                    this->get_shape_operations(doc["shape1"]),
                    this->get_shape_operations(doc["shape2"])
                );
    } else if(type == "intersection"){
        this->log_data(doc, "shape1", TYPE_OBJECT, true);
        this->log_data(doc, "shape2", TYPE_OBJECT, true);

        return std::make_unique<shape_op::Intersection>(
                    this->get_shape_operations(doc["shape1"]),
                    this->get_shape_operations(doc["shape2"])
                );
    } else if(type == "difference"){
        this->log_data(doc, "shape1", TYPE_OBJECT, true);
        this->log_data(doc, "shape2", TYPE_OBJECT, true);

        return std::make_unique<shape_op::Difference>(
                    this->get_shape_operations(doc["shape1"]),
                    this->get_shape_operations(doc["shape2"])
                );
    }

    return nullptr;
}

std::unique_ptr<Simulation> ProjectData::load_simulation(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc){
    this->log_data(doc, "type", TYPE_STRING, true);
    std::string type = doc["type"].GetString();
    if(type == "marginal_bone_loss"){ 
        this->log_data(doc, "time_step", TYPE_DOUBLE, true);
        this->log_data(doc, "maximum_volume_variation", TYPE_DOUBLE, true);
        this->log_data(doc, "time_limit", TYPE_DOUBLE, true);
        this->log_data(doc, "maturation_rate", TYPE_DOUBLE, true);
        this->log_data(doc, "pc", TYPE_DOUBLE, true);

        this->log_data(doc, "mandible_geometry", TYPE_INT, true);

        this->log_data(doc, "traction", TYPE_ARRAY, true);
        this->log_data(doc, "compression", TYPE_ARRAY, true);
        this->log_data(doc, "shear", TYPE_ARRAY, true);

        this->log_data(doc, "a0", TYPE_ARRAY, true);
        this->log_data(doc, "a1", TYPE_ARRAY, true);
        this->log_data(doc, "a2", TYPE_ARRAY, true);
        this->log_data(doc, "a3", TYPE_ARRAY, true);

        double time_step = doc["time_step"].GetDouble();
        double maximum_volume_variation = doc["maximum_volume_variation"].GetDouble();
        double time_limit = doc["time_limit"].GetDouble();
        double maturation_rate = doc["maturation_rate"].GetDouble();
        double pc = doc["pc"].GetDouble();

        size_t mandible_geometry = doc["mandible_geometry"].GetInt();

        auto traction = doc["traction"].GetArray();
        auto compression = doc["compression"].GetArray();
        auto shear = doc["shear"].GetArray();
        auto a0 = doc["a0"].GetArray();
        auto a1 = doc["a1"].GetArray();
        auto a2 = doc["a2"].GetArray();
        auto a3 = doc["a3"].GetArray();

        const size_t RANGE_NUM = simulation::MarginalBoneLoss::RANGE_NUM;

        logger::log_assert(traction.Size() == RANGE_NUM, logger::ERROR, "\"traction\" item must have size {}.", RANGE_NUM);
        logger::log_assert(compression.Size() == RANGE_NUM, logger::ERROR, "\"compression\" item must have size {}.", RANGE_NUM);
        logger::log_assert(shear.Size() == RANGE_NUM, logger::ERROR, "\"shear\" item must have size {}.", RANGE_NUM);

        simulation::MarginalBoneLoss::Range t;
        simulation::MarginalBoneLoss::Range c;
        simulation::MarginalBoneLoss::Range s;
        for(size_t i = 0; i < RANGE_NUM; ++i){
            t[i] = traction[i].GetDouble();
            c[i] = compression[i].GetDouble();
            s[i] = shear[i].GetDouble();
        }

        std::vector<double> a0v(a0.Size());
        std::vector<double> a1v(a1.Size());
        std::vector<double> a2v(a2.Size());
        std::vector<double> a3v(a3.Size());

        const auto load_poly = [&](std::vector<double>& a, decltype(a0) array){
            for(size_t i = 0; i < array.Size(); ++i){
                a[i] = array[i].GetDouble();
            }
        };

        load_poly(a0v, a0);
        load_poly(a1v, a1);
        load_poly(a2v, a2);
        load_poly(a3v, a3);

        return std::make_unique<simulation::MarginalBoneLoss>(this->topopt_mesher.get(), this->topopt_fea.get(), this->geometries[mandible_geometry].get(), mandible_geometry, time_step, maximum_volume_variation, time_limit, maturation_rate, pc, t, c, s, a0v, a1v, a2v, a3v, this->type);
    }

    return nullptr;
}
