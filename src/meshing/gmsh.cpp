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

#include "meshing/gmsh.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include <BRepBuilderAPI_Transform.hxx>
#include <cmath>
#include <gmsh.h>
#include <algorithm>
#include <iomanip>
#include <limits>
#include "project_data.hpp"
#include <unordered_map>
#include <BOPAlgo_Splitter.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <set>

namespace meshing{

Gmsh::Gmsh(const std::vector<std::unique_ptr<Geometry>>& geometries,
           const MeshElementFactory* const elem_type,
           double size, double thickness, double tmp_scale, int algorithm2D, int algorithm3D):
    Meshing(geometries, elem_type, thickness),
    size(tmp_scale*size), algorithm2D(algorithm2D), algorithm3D(algorithm3D), tmp_scale(tmp_scale){
}

void Gmsh::mesh(const std::vector<Force>& forces, 
                const std::vector<Support>& supports){
    TopoDS_Shape shape = this->make_compound(this->geometries);

    bool has_condition_inside = false;
    auto problem_type = this->elem_info->get_problem_type();
    TopoDS_Shape sh;
    if(this->tmp_scale == 1.0){
        sh = BRepBuilderAPI_Copy(shape);
    } else {
        sh = this->make_compound(this->geometries, this->tmp_scale);
    }
    if(problem_type == utils::PROBLEM_TYPE_2D){
        has_condition_inside = this->adapt_for_boundary_condition_inside(shape, forces, supports);
    } else {
        // Maybe disable for cubic elements?
        for(auto& f:forces){
            has_condition_inside = this->is_strictly_inside3D(f.S.get_centroid(), shape);
            BOPAlgo_Splitter splitter;
            splitter.SetNonDestructive(true);
            splitter.AddArgument(sh);
            splitter.AddTool(f.S.get_shape());
            splitter.Perform();
            sh = splitter.Shape();
        }
        for(auto& s:supports){
            has_condition_inside = this->is_strictly_inside3D(s.S.get_centroid(), shape);
            BOPAlgo_Splitter splitter;
            splitter.SetNonDestructive(true);
            splitter.AddArgument(sh);
            splitter.AddTool(s.S.get_shape());
            splitter.Perform();
            sh = splitter.Shape();
        }
    }
    // if(has_condition_inside){
    //     // Workaround so that this does not break current (faster) method of
    //     // distributing elements to the different geometry instances, at least
    //     // considering global meshing for linear analysis.
    //     //
    //     // Its current use is mostly for simple topopt/beam sizing, so it's a good
    //     // workaround for now.
    //     //
    //     // Otherwise, the simpler idea would be to test for the center of mass
    //     // of each geometry, which would be a problem if the geometry has a 
    //     // hole in its center.
    //     //
    //     // Seemed better to sacrifice this gimmick than something more useful
    //     // such as support for complex geometries.
    //     logger::log_assert(this->geometries.size() == 1, logger::ERROR, "applying boundary conditions inside a geometry is currently not supported when the number of geometries is greater than 1.");
    // }

    std::vector<size_t> geom_elem_mapping, elem_node_tags, bound_elem_node_tags;
    auto id_map = this->gmsh_meshing(has_condition_inside, sh, geom_elem_mapping, elem_node_tags, bound_elem_node_tags, this->elem_info);

    bool deduplicate = false;
    if(geometries.size() > 1){
        deduplicate = true;
    }
    this->generate_elements(shape,
                            geom_elem_mapping, 
                            elem_node_tags, 
                            bound_elem_node_tags,
                            id_map,
                            forces, 
                            supports,
                            deduplicate,
                            has_condition_inside);
}

std::unordered_map<size_t, MeshNode*> Gmsh::gmsh_meshing(bool has_condition_inside, TopoDS_Shape sh, std::vector<size_t>& geom_elem_mapping, std::vector<size_t>& elem_node_tags, std::vector<size_t>& bound_elem_node_tags, const MeshElementFactory* const elem_type){
    gmsh::initialize();

    gmsh::model::add("base");

    gmsh::vectorpair vec;
    gmsh::model::occ::importShapesNativePointer(static_cast<const void*>(&sh), vec);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMin", 0);
    gmsh::option::setNumber("Mesh.MeshSizeMax", this->size);

    gmsh::option::setNumber("Mesh.Algorithm", this->algorithm2D);
    gmsh::option::setNumber("Mesh.Algorithm3D", this->algorithm3D);

    gmsh::option::setNumber("Mesh.ElementOrder", this->elem_info->get_element_order());
    gmsh::option::setNumber("Mesh.Optimize", 1);

    gmsh::option::setNumber("Mesh.AngleToleranceFacetOverlap", 0.00001);

    // Quad/hex recombination
    auto problem_type = this->elem_info->get_problem_type();
    if(elem_type->get_shape_type() == Element::Shape::QUAD){
        gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 2);
        gmsh::option::setNumber("Mesh.RecombineAll", 1);
        gmsh::option::setNumber("Mesh.Recombine3DAll", 1);
    }

    size_t dim = 0;
    if(problem_type == utils::PROBLEM_TYPE_2D){
        dim = 2;
    } else if(problem_type == utils::PROBLEM_TYPE_3D){
        dim = 3;
    }

    gmsh::model::mesh::generate(dim);

    gmsh::write("mesh.msh");

    size_t type = elem_type->get_gmsh_element_type();
    size_t bound_type = elem_type->get_boundary_gmsh_element_type();
    std::vector<int> elem_types;
    gmsh::model::mesh::getElementTypes(elem_types);
    // Check if meshing went well
    logger::log_assert(std::find(elem_types.begin(), elem_types.end(), type) != elem_types.end(), logger::ERROR,
                        "element type not found in mesh's list of element types (this shouldn't happen).");
    logger::log_assert(std::find(elem_types.begin(), elem_types.end(), bound_type) != elem_types.end(), logger::ERROR,
                        "element type of boundary elements not found in mesh's list of element types (this shouldn't happen).");

    // Get nodes and their tags
    std::vector<std::size_t> node_tags, boundary_node_tags;
    std::vector<double> node_coords, node_params;
    if(has_condition_inside){
        gmsh::model::mesh::getNodes(node_tags, node_coords, node_params, -1, -1, true);
    } else {
        gmsh::model::mesh::getNodes(node_tags, node_coords, node_params, dim, -1, true);
        std::vector<double> bnode_coords, bnode_params;
        gmsh::model::mesh::getNodes(boundary_node_tags, bnode_coords, bnode_params, dim-1, -1, true);
    }

    // Would need to be changed to support multiple elements
    geom_elem_mapping.resize(geometries.size(),0);
    if(geom_elem_mapping.size() == 1){
        std::vector<size_t> elem_tags;
        gmsh::model::mesh::getElementsByType(type, elem_tags, elem_node_tags, -1);
        geom_elem_mapping[0] = elem_tags.size();
        elem_tags.clear();
        gmsh::model::mesh::getElementsByType(bound_type, elem_tags, bound_elem_node_tags);
    } else {
        std::vector<size_t> elem_tags;
        // Boundary elements
        gmsh::model::mesh::getElementsByType(bound_type, elem_tags, bound_elem_node_tags);
        elem_tags.clear();
        // I'm assuming here that the geometry's id and the internal model id
        // in Gmsh/OCCT are the same.
        size_t end = gmsh::model::occ::getMaxTag(dim);
        if(end == geometries.size()){
            for(size_t i = 1; i <= end; ++i){
                std::vector<size_t> elem_tags_tmp;
                std::vector<size_t> elem_node_tags_tmp;
                gmsh::model::mesh::getElementsByType(type, elem_tags_tmp, elem_node_tags_tmp, i);
                elem_node_tags.insert(elem_node_tags.end(), elem_node_tags_tmp.begin(), elem_node_tags_tmp.end());
                geom_elem_mapping[i-1] = elem_tags_tmp.size();
            }
            for(size_t j = 1; j < geom_elem_mapping.size(); ++j){
                geom_elem_mapping[j] += geom_elem_mapping[j-1];
            }
        } else {
            const size_t nodes_per_elem = elem_type->get_nodes_per_element();
            // Sometimes Gmsh finds more volumes than there are input geometries,
            // so this is used in that case
            for(size_t i = 1; i <= end; ++i){
                std::vector<size_t> elem_tags_tmp;
                std::vector<size_t> elem_node_tags_tmp;
                double x, y, z;
                gmsh::model::occ::getCenterOfMass(dim, i, x, y, z);
                gp_Pnt p(x, y, z);
                for(size_t j = 0; j < geometries.size(); ++j){
                    const auto& g = geometries[j];
                    if(g->is_inside(p)){
                        gmsh::model::mesh::getElementsByType(type, elem_tags_tmp, elem_node_tags_tmp, i);
                        elem_node_tags.insert(elem_node_tags.begin()+geom_elem_mapping[j]*nodes_per_elem, elem_node_tags_tmp.begin(), elem_node_tags_tmp.end());
                        geom_elem_mapping[j] += elem_tags_tmp.size();
                        break;
                    }
                }
            }
            for(size_t j = 1; j < geom_elem_mapping.size(); ++j){
                geom_elem_mapping[j] += geom_elem_mapping[j-1];
            }
        }
    }

    gmsh::clear();
    gmsh::finalize();

    const size_t dof = elem_type->get_dof_per_node();

    const auto node_comp = [](const MeshNode& n1, const MeshNode& n2) -> bool{
        return n1.id < n2.id;
    };
    std::set<MeshNode, decltype(node_comp)> node_set(node_comp);
    for(size_t i = 0; i < node_tags.size(); ++i){
        gp_Pnt p(node_coords[i*3]/this->tmp_scale, node_coords[i*3+1]/this->tmp_scale, node_coords[i*3+2]/this->tmp_scale);

        node_set.emplace(p, node_tags[i], dof);
    }

    std::unordered_map<size_t, MeshNode*> id_map;
    id_map.reserve(node_set.size());

    this->node_list.clear();
    this->node_list.reserve(node_set.size());
    for(auto it = node_set.cbegin(); it != node_set.cend(); ++it){

        this->node_list.emplace_back(std::make_unique<MeshNode>(it->point, it->id, dof));
        id_map.emplace(it->id, this->node_list.back().get());
    }
    node_set.clear();

    std::set<size_t> boundary_node_set(boundary_node_tags.begin(), boundary_node_tags.end());

    this->boundary_node_list.clear();
    this->boundary_node_list.reserve(boundary_node_set.size());
    for(auto it = boundary_node_set.cbegin(); it != boundary_node_set.cend(); ++it){
        this->boundary_node_list.push_back(id_map.at(*it));
    }

    return id_map;
}

}
