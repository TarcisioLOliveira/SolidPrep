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
#include <gmsh.h>
#include <algorithm>
#include <limits>
#include "project_data.hpp"
#include <unordered_map>
#include <BOPAlgo_Splitter.hxx>
#include <BRepBuilderAPI_Copy.hxx>

namespace meshing{

Gmsh::Gmsh(const std::vector<std::unique_ptr<Geometry>>& geometries,
           const MeshElementFactory* const elem_type,
           double size, double thickness, int algorithm2D, int algorithm3D):
    Meshing(geometries, elem_type, thickness),
    size(size), algorithm2D(algorithm2D), algorithm3D(algorithm3D){
}

void Gmsh::mesh(const std::vector<Force>& forces, 
                const std::vector<Support>& supports){
    TopoDS_Shape shape = this->make_compound(this->geometries);

    bool has_condition_inside = false;
    // Workaround so that this does not break current (faster) method of
    // distributing elements to the different geometry instances, at least
    // considering global meshing for linear analysis.
    //
    // Its current use is mostly for simple topopt/beam sizing, so it's a good
    // workaround for now.
    //
    // Otherwise, the simpler idea would be to test for the center of mass
    // of each geometry, which would be a problem if the geometry has a 
    // hole in its center.
    //
    // Seemed better to sacrifice this gimmick than something more useful
    // such as support for complex geometries.
    auto problem_type = this->elem_info->get_problem_type();
    if(this->geometries.size() == 1 && problem_type == utils::PROBLEM_TYPE_2D){
        has_condition_inside = this->adapt_for_boundary_condition_inside(shape, forces, supports);
    } else {
        TopoDS_Shape sh = BRepBuilderAPI_Copy(shape);
        for(auto& f:forces){
            BOPAlgo_Splitter splitter;
            splitter.SetNonDestructive(true);
            splitter.AddArgument(sh);
            splitter.AddTool(f.S.get_shape());
            splitter.Perform();
            sh = splitter.Shape();
        }
        for(auto& s:supports){
            BOPAlgo_Splitter splitter;
            splitter.SetNonDestructive(true);
            splitter.AddArgument(sh);
            splitter.AddTool(s.S.get_shape());
            splitter.Perform();
            sh = splitter.Shape();
        }
        shape = sh;
    }

    std::vector<size_t> geom_elem_mapping, elem_tags, elem_node_tags;
    auto id_map = this->gmsh_meshing(has_condition_inside, shape, geom_elem_mapping, elem_tags, elem_node_tags, this->elem_info);

    if(geometries.size() > 1){
       this->find_duplicates(id_map);
    }

    size_t nodes_per_elem = this->elem_info->get_nodes_per_element();

    auto list = this->generate_element_shapes(elem_tags, elem_node_tags, nodes_per_elem, id_map);

    this->optimize(list, has_condition_inside);

    this->prepare_for_FEM(shape, geom_elem_mapping, list, forces, supports);
}

std::unordered_map<size_t, MeshNode*> Gmsh::gmsh_meshing(bool has_condition_inside, TopoDS_Shape sh, std::vector<size_t>& geom_elem_mapping, std::vector<size_t>& elem_tags, std::vector<size_t>& elem_node_tags, const MeshElementFactory* const elem_type){
    gmsh::initialize();

    gmsh::model::add("base");

    gmsh::vectorpair vec;
    gmsh::model::occ::importShapesNativePointer(static_cast<const void*>(&sh), vec);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMin", this->size);
    gmsh::option::setNumber("Mesh.MeshSizeMax", this->size);

    gmsh::option::setNumber("Mesh.Algorithm", this->algorithm2D);
    gmsh::option::setNumber("Mesh.Algorithm3D", this->algorithm3D);

    gmsh::option::setNumber("Mesh.ElementOrder", this->elem_info->get_element_order());
    gmsh::option::setNumber("Mesh.Optimize", 1);

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

    size_t type = elem_type->get_gmsh_element_type();
    std::vector<int> elem_types;
    gmsh::model::mesh::getElementTypes(elem_types);
    // Check if meshing went well
    logger::log_assert(std::find(elem_types.begin(), elem_types.end(), type) != elem_types.end(), logger::ERROR,
                        "element type not found in mesh's list of element types (this shouldn't happen).");

    // Get nodes and their tags
    std::vector<std::size_t> node_tags;
    std::vector<double> node_coords, node_params;
    if(has_condition_inside){
        gmsh::model::mesh::getNodes(node_tags, node_coords, node_params, -1, -1, true);
    } else {
        gmsh::model::mesh::getNodes(node_tags, node_coords, node_params, dim, -1, true);
    }

    // Would need to be changed to support multiple elements
    geom_elem_mapping.resize(geometries.size());
    if(geometries.size() == 1){
        gmsh::model::mesh::getElementsByType(type, elem_tags, elem_node_tags, -1);
        geom_elem_mapping[0] = elem_tags.size();
    } else {
        // I'm assuming here that the geometry's id and the internal model id
        // in Gmsh/OCCT are the same.
        size_t end = gmsh::model::occ::getMaxTag(dim);
        for(size_t i = 1; i <= end; ++i){
            std::vector<size_t> elem_tags_tmp;
            std::vector<size_t> elem_node_tags_tmp;
            gmsh::model::mesh::getElementsByType(type, elem_tags_tmp, elem_node_tags_tmp, i);
            elem_tags.insert(elem_tags.end(), elem_tags_tmp.begin(), elem_tags_tmp.end());
            elem_node_tags.insert(elem_node_tags.end(), elem_node_tags_tmp.begin(), elem_node_tags_tmp.end());
            geom_elem_mapping[i-1] = elem_tags.size();
        }
    }

    gmsh::clear();
    gmsh::finalize();

    size_t dof = elem_type->get_dof_per_node();

    std::unordered_map<size_t, MeshNode*> id_map;
    id_map.reserve(node_tags.size());

    this->node_list.clear();
    this->node_list.reserve(node_tags.size());
    for(size_t i = 0; i < node_tags.size(); ++i){
        gp_Pnt p(node_coords[i*3], node_coords[i*3+1], node_coords[i*3+2]);
        this->node_list.emplace_back(std::make_unique<MeshNode>(p, node_tags[i], dof));
        id_map.emplace(node_tags[i], this->node_list[i].get());
    }

    return id_map;
}

}
