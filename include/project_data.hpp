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
#ifndef PROJECT_DATA_HPP
#define PROJECT_DATA_HPP

#include <json/value.h>
#include <memory>
#include <string>
#include <vector>
#include "meshing/mesh_file.hpp"
#include "project_specification/data_map.hpp"
#include "project_specification/registry.hpp"
#include "projection.hpp"
#include "utils/delayed_pointer.hpp"
#include "json/value.h"

#include "shape_handler.hpp"
#include "support.hpp"
#include "pathfinding.hpp"
#include "sizing.hpp"
#include "material.hpp"
#include "utils.hpp"
#include "finite_element.hpp"
#include "topology_optimization.hpp"
#include "geometry.hpp"
#include "density_filter.hpp"
#include "function.hpp"
#include "optimizer.hpp"
#include "spring.hpp"
#include "internal_loads.hpp"
#include "field.hpp"
#include "solver_manager.hpp"
#include "sub_problem.hpp"
#include "simulation.hpp"

/**
 * Reads and stores project data.
 */
class ProjectData {
    public:

    struct ContactData{
        FiniteElement::ContactType contact_type = FiniteElement::ContactType::RIGID;
        double rtol_abs = 0;
        double max_step = 0;
        double EPS_DISPL = 0;
    };

    //enum DataType{
    //    TYPE_NULL,
    //    TYPE_BOOL,
    //    TYPE_INT,
    //    TYPE_DOUBLE,
    //    TYPE_STRING,
    //    TYPE_ARRAY,
    //    TYPE_OBJECT
    //};
    /**
     * Loads project data file.
     *
     * @param project_file Path to file.
     */
    ProjectData(std::string project_file);

    // Analysis type
    bool do_meshing = false;
    bool generate_beams = false;
    bool do_topopt = false;
    bool do_fea = false;
    bool do_shape_opt = false;
    bool do_simulation = false;

    double thickness;
    std::unique_ptr<Pathfinding> pathfinder;
    std::unique_ptr<Sizing> sizer;
    std::vector<utils::DelayedPointer<Material>> materials;
    utils::ProblemType type;
    std::vector<std::unique_ptr<Geometry>> geometries;
    std::vector<Force> forces;
    std::vector<Support> supports;
    std::vector<Spring> springs;
    std::vector<InternalLoads> internal_loads;
    std::unique_ptr<FiniteElement> sizer_fea;
    std::unique_ptr<TopologyOptimization> topopt;
    std::unique_ptr<SolverManager> topopt_fea;
    std::unique_ptr<Meshing> topopt_mesher;
    MeshElementFactory* topopt_element;
    std::unique_ptr<DensityFilter> density_filter;
    std::unique_ptr<Projection> projection;
    std::unique_ptr<DensityBasedOptimizer> topopt_optimizer;
    std::unique_ptr<NodeShapeBasedOptimizer> shopt_optimizer;
    std::unique_ptr<Simulation> simulator;
    std::vector<utils::DelayedPointer<Field>> fields;
    std::vector<SubProblem> sub_problems;
    ContactData contact_data;
    meshing::MeshFile* mesh_file = nullptr;
    double topopt_penalization = 3; // Topopt with void
    double topopt_psi = 0.5; // Topopt with two materials
    Projection::Parameter projection_parameters;
    ShapeHandler shape_handler;
    
    std::string folder_path;
    std::string element_name;

    utils::DelayedPointerView<Material> get_material(const std::string& name);

    CrossSection get_cross_section(const Json::Value& doc) const;
    
    private:
    std::unique_ptr<meshing::MeshFile> mesh_file_internal = nullptr;

    std::string get_folder_path(const std::string& project_file_path) const;

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
    bool log_data(const Json::Value& doc, std::string name, projspec::DataType type, bool required) const;

    ContactData get_contact_data(const Json::Value& doc) const;

    std::vector<Force> get_loads(const Json::Value& doc) const;

    std::vector<Support> get_support(const Json::Value& doc) const;

    std::vector<std::unique_ptr<Geometry>> load_geometries(const Json::Value& doc);

    projspec::DataMap generate_stub(const projspec::ObjectRequirements& reqs,
                                    const projspec::RequirementConditions& conds);

    projspec::DataMap generate_data_map(const Json::Value& doc,
                                        const projspec::ObjectRequirements& reqs,
                                        const projspec::RequirementConditions& conds);

    std::vector<utils::DelayedPointer<Material>> generate_material_stubs(const Json::Value& doc,
                                    const projspec::RequirementConditions& conds);

    void load_materials(const Json::Value& doc, std::vector<utils::DelayedPointer<Material>>& vec,
                                    const projspec::RequirementConditions& conds);

    std::vector<utils::DelayedPointer<Field>> generate_field_stubs(const Json::Value& doc,
                                    const projspec::RequirementConditions& conds);

    void load_fields(const Json::Value& doc, std::vector<utils::DelayedPointer<Field>>& vec,
                                    const projspec::RequirementConditions& conds);

    std::vector<Spring> get_springs(const Json::Value& doc);

    std::vector<InternalLoads> get_internal_loads(const Json::Value& doc);

    std::vector<SubProblem> load_sub_problems(const Json::Value& doc);

    std::unique_ptr<FiniteElement> load_fea(const Json::Value& doc,
                                            const projspec::RequirementConditions& conds);

    std::unique_ptr<Pathfinding> load_pathfinder(const Json::Value& doc,
                                                 const projspec::RequirementConditions& conds);

    std::unique_ptr<Sizing> load_sizer(const Json::Value& doc,
                                       const projspec::RequirementConditions& conds);

    std::unique_ptr<Meshing> load_mesher(const Json::Value& doc,
                                         const projspec::RequirementConditions& conds);

    template<class T>
    std::unique_ptr<T> get_function(const Json::Value& doc,
                                    const projspec::RequirementConditions& conds);

    template<class T>
    void get_constraints(const Json::Value& doc,
                         const projspec::RequirementConditions& conds,
                         std::vector<T>& functions);

    template<class T>
    void get_objective_functions(const Json::Value& doc,
                                 const projspec::RequirementConditions& conds,
                                 std::vector<std::unique_ptr<T>>& functions,
                                 std::vector<double>& weights);

    std::unique_ptr<DensityBasedFunction> get_topopt_function(const Json::Value& doc,
                                                              const projspec::RequirementConditions& conds);

    void get_topopt_constraints(const Json::Value& doc,
                                const projspec::RequirementConditions& conds,
                                std::vector<DensityBasedConstraint>& functions);

    void get_topopt_objective_functions(const Json::Value& doc,
                                        const projspec::RequirementConditions& conds,
                                        std::vector<std::unique_ptr<DensityBasedFunction>>& functions,
                                        std::vector<double>& weights);

    std::unique_ptr<DensityBasedOptimizer> load_topopt_optimizer(const Json::Value& doc,
                                                                 const projspec::RequirementConditions& conds);

    std::unique_ptr<NodeShapeBasedFunction> get_shopt_function(const Json::Value& doc,
                                                               const projspec::RequirementConditions& conds);

    void get_shopt_constraints(const Json::Value& doc,
                               const projspec::RequirementConditions& conds,
                               std::vector<NodeShapeBasedConstraint>& functions);

    void get_shopt_objective_functions(const Json::Value& doc,
                                       const projspec::RequirementConditions& conds,
                                       std::vector<std::unique_ptr<NodeShapeBasedFunction>>& functions,
                                       std::vector<double>& weights);

    std::unique_ptr<NodeShapeBasedOptimizer> load_shopt_optimizer(const Json::Value& doc,
                                                                  const projspec::RequirementConditions& conds);

    std::unique_ptr<DensityFilter> load_density_filter(const Json::Value& doc,
                                                       const projspec::RequirementConditions& conds);

    std::unique_ptr<Projection> load_projection(const Json::Value& doc,
                                                const projspec::RequirementConditions& conds);

    Projection::Parameter get_projection_parameter(const Json::Value& doc) const;

    std::unique_ptr<shape_op::ShapeOp> get_shape_operations(const Json::Value& doc) const;

    std::unique_ptr<Simulation> load_simulation(const Json::Value& doc,
                                                const projspec::RequirementConditions& conds);

    std::unique_ptr<MeshElementFactory> get_element_type(const std::string& name);
};

template<class T>
std::unique_ptr<T> ProjectData::get_function(const Json::Value& data,
                                             const projspec::RequirementConditions& conds){
    auto& json_data = data;
    projspec::ObjectRequirements reqs;
    this->log_data(json_data, "type", projspec::TYPE_STRING, true);
    auto type_name = json_data["type"].asString();
    logger::log_assert(projspec::Registry::get(T::get_name(), type_name, reqs),
            logger::ERROR, "unknown function: {}", type_name);
    projspec::DataMap data_map(this->generate_data_map(json_data, reqs, conds));
    return projspec::Factory<T>::construct(type_name, data_map);
}



template<class T>
void ProjectData::get_constraints(const Json::Value& doc,
                                  const projspec::RequirementConditions& conds,
                                  std::vector<T>& functions){

    this->log_data(doc, "constraints", projspec::TYPE_ARRAY, true);
    auto array = doc["constraints"];
    typedef typename T::FunctionType F;
    functions.reserve(array.size());
    for(auto& f:array){
        std::vector<Constraint::Type> types;
        std::vector<double> bounds;

        this->log_data(f, "type", projspec::TYPE_STRING, true);
        std::string type = f["type"].asString();
        if(this->log_data(f, "less_than", projspec::TYPE_DOUBLE, false)){
            types.push_back(Constraint::Type::LESS_THAN);
            bounds.push_back(f["less_than"].asDouble());
        }
        if(this->log_data(f, "greater_than", projspec::TYPE_DOUBLE, false)){
            types.push_back(Constraint::Type::GREATER_THAN);
            bounds.push_back(f["greater_than"].asDouble());
        }
        if(this->log_data(f, "equals", projspec::TYPE_DOUBLE, false)){
            types.push_back(Constraint::Type::EQUAL);
            bounds.push_back(f["equals"].asDouble());
        }
        logger::log_assert(types.size() > 0, logger::ERROR, "constraint {} is missing at least one of the following: \"equals\", \"greater_than\", \"less_than\"", f["type"].asString());

        functions.emplace_back(this->get_function<F>(f, conds), types, bounds);
    }
}

template<class T>
void ProjectData::get_objective_functions(const Json::Value& doc,
                                          const projspec::RequirementConditions& conds,
                                          std::vector<std::unique_ptr<T>>& functions,
                                          std::vector<double>& weights){
    this->log_data(doc, "objective", projspec::TYPE_ARRAY, true);
    auto array = doc["objective"];
    functions.reserve(array.size());
    weights.reserve(array.size());
    for(auto& f:array){
        functions.push_back(this->get_function<T>(f, conds));
        if(this->log_data(f, "weight", projspec::TYPE_DOUBLE, false)){
            weights.push_back(f["weight"].asDouble());
        } else {
            weights.push_back(1.0);
        }
    }
}

#endif
