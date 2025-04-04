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

#include <memory>
#include <string>
#include <vector>
#include "meshing/mesh_file.hpp"
#include "projection.hpp"
#include "rapidjson/document.h"

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
    std::vector<std::unique_ptr<Material>> materials;
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
    std::unique_ptr<MeshElementFactory> topopt_element;
    std::unique_ptr<BoundaryMeshElementFactory> topopt_boundary_element;
    std::unique_ptr<DensityFilter> density_filter;
    std::unique_ptr<Projection> projection;
    std::unique_ptr<DensityBasedOptimizer> topopt_optimizer;
    std::unique_ptr<NodeShapeBasedOptimizer> shopt_optimizer;
    std::unique_ptr<Simulation> simulator;
    std::vector<std::unique_ptr<Field>> fields;
    std::vector<SubProblem> sub_problems;
    ContactData contact_data;
    meshing::MeshFile* mesh_file = nullptr;
    
    std::string folder_path;
    std::string element_name;
    
    private:
    std::unique_ptr<meshing::MeshFile> mesh_file_internal = nullptr;
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
    bool log_data(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc, std::string name, DataType type, bool required) const;

    std::string get_folder_path(const std::string& project_file_path) const;

    ContactData get_contact_data(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::vector<std::unique_ptr<Material>> load_materials_simple(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::vector<std::unique_ptr<Material>> load_materials_field(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::vector<std::unique_ptr<Geometry>> load_geometries(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    void assign_materials(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<Pathfinding> load_pathfinder(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<Sizing> load_sizer(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<FiniteElement> load_fea(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::vector<SubProblem> load_sub_problems(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<Meshing> load_mesher(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<DensityFilter> load_density_filter(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<Projection> load_projection(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<DensityBasedOptimizer> load_topopt_optimizer(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<NodeShapeBasedOptimizer> load_shopt_optimizer(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::vector<std::unique_ptr<Field>> load_fields(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<DensityBasedFunction> get_topopt_function(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc, double pc, double psiK);

    std::unique_ptr<NodeShapeBasedFunction> get_shopt_function(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    void get_topopt_constraints(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc, double pc, double psi, std::vector<DensityBasedConstraint>& functions);

    void get_topopt_objective_functions(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc, double pc, double psi, std::vector<std::unique_ptr<DensityBasedFunction>>& functions, std::vector<double>& weights);

    void get_shopt_constraints(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc, std::vector<NodeShapeBasedConstraint>& functions);

    void get_shopt_objective_functions(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc, std::vector<std::unique_ptr<NodeShapeBasedFunction>>& functions, std::vector<double>& weights);

    std::unique_ptr<MeshElementFactory> get_element_type(const std::string& name);

    std::vector<Force> get_loads(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::vector<Support> get_support(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::vector<Spring> get_springs(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::vector<InternalLoads> get_internal_loads(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    Projection::Parameter get_projection_parameter(const rapidjson::GenericValue<rapidjson::UTF8<>>& p) const;

    CrossSection get_cross_section(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc) const;

    std::unique_ptr<shape_op::ShapeOp> get_shape_operations(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc) const;

    std::unique_ptr<Simulation> load_simulation(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);
};

#endif
