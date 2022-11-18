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
#include "projection.hpp"
#include "rapidjson/document.h"
#include "rapidjson/rapidjson.h"

#include "support.hpp"
#include "pathfinding.hpp"
#include "sizing.hpp"
#include "material.hpp"
#include "utils.hpp"
#include "finite_element.hpp"
#include "topology_optimization.hpp"
#include "geometry.hpp"
#include "density_filter.hpp"

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
    std::vector<std::unique_ptr<Material>> materials;
    utils::ProblemType type;
    std::vector<std::unique_ptr<Geometry>> geometries;
    std::vector<Force> forces;
    std::vector<Support> supports;
    std::unique_ptr<FiniteElement> sizer_fea;
    std::unique_ptr<TopologyOptimization> topopt;
    std::unique_ptr<FiniteElement> topopt_fea;
    std::unique_ptr<Meshing> topopt_mesher;
    std::unique_ptr<MeshElementFactory> topopt_element;
    std::unique_ptr<DensityFilter> density_filter;
    std::unique_ptr<Projection> projection;
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

    std::vector<std::unique_ptr<Material>> load_materials(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::vector<std::unique_ptr<Geometry>> load_geometries(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc, const std::string& folder_path);

    std::unique_ptr<Pathfinding> load_pathfinder(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<Sizing> load_sizer(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<FiniteElement> load_fea(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<Meshing> load_mesher(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<TopologyOptimization> load_topopt(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<DensityFilter> load_density_filter(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<Projection> load_projection(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::unique_ptr<MeshElementFactory> get_element_type(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::vector<Force> get_loads(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    std::vector<Support> get_support(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc);

    Projection::Parameter get_projection_parameter(const rapidjson::GenericValue<rapidjson::UTF8<>>& p) const;

    CrossSection get_cross_section(const rapidjson::GenericValue<rapidjson::UTF8<>>& doc) const;
};




#endif
