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

#include <algorithm>
#include <memory>
#include <string>
#include <vector>
#include "finite_element/direct_solver.hpp"
#include "meshing/gmsh.hpp"
#include "rapidjson/allocators.h"
#include "rapidjson/document.h"
#include "STEPCAFControl_Reader.hxx"
#include "logger.hpp"
#include "force.hpp"
#include "rapidjson/encodings.h"
#include "sizing/standard_sizing.hpp"
#include "support.hpp"
#include "pathfinding.hpp"
#include "sizing.hpp"
#include "ground_structure.hpp"
#include "material.hpp"
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
#include "finite_element.hpp"
#include "topology_optimization.hpp"
#include "topology_optimization/minimal_volume.hpp"
#include "topology_optimization/minimal_compliance.hpp"

/**
 * Reads and stores project data.
 */
class ProjectData {
    public:

    enum DataType{
        TYPE_NULL,
        TYPE_BOOL,
        TYPE_INT,
        TYPE_DOUBLE,
        TYPE_STRING,
        TYPE_ARRAY,
        TYPE_OBJECT
    };
    enum AnalysisType{
        COMPLETE,
        FEA_ONLY,
        OPTIMIZE_ONLY,
        BEAMS_ONLY
    };
    /**
     * Loads project data file.
     *
     * @param project_file Path to file.
     */
    ProjectData(std::string project_file);

    double thickness;
    std::unique_ptr<Pathfinding> pathfinder;
    std::unique_ptr<Sizing> sizer;
    std::unique_ptr<Material> material;
    utils::ProblemType type;
    std::unique_ptr<GroundStructure> ground_structure;
    std::vector<Force> forces;
    std::vector<Support> supports;
    std::unique_ptr<FiniteElement> sizer_fea;
    std::unique_ptr<TopologyOptimization> topopt;
    std::unique_ptr<FiniteElement> topopt_fea;
    std::unique_ptr<Meshing> topopt_mesher;
    std::unique_ptr<MeshElementFactory> topopt_element;
    AnalysisType analysis;
    
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
    bool log_data(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc, std::string name, DataType type, bool required) const;

    std::unique_ptr<Material> load_material(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<Pathfinding> load_pathfinder(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<Sizing> load_sizer(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<FiniteElement> load_fea(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<Meshing> load_mesher(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<TopologyOptimization> load_topopt(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<MeshElementFactory> get_element_type(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::vector<Force> get_loads(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::vector<Support> get_support(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);
};




#endif
