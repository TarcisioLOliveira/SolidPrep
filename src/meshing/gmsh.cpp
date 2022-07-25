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

namespace meshing{

Gmsh::Gmsh(const std::vector<std::unique_ptr<Geometry>>& geometries,
           const MeshElementFactory* const elem_type,
           double size, double thickness, int algorithm):
    Meshing(geometries, elem_type, thickness),
    size(size), algorithm(algorithm){
}

void Gmsh::mesh(const std::vector<Force>& forces, 
                const std::vector<Support>& supports){
    TopoDS_Shape shape = this->make_compound(this->geometries);
    bool has_condition_inside = this->adapt_for_boundary_condition_inside(shape, forces, supports);

    std::vector<size_t> elem_tags, elem_node_tags;
    this->gmsh_meshing(has_condition_inside, shape, elem_tags, elem_node_tags, this->elem_info);

    std::unordered_map<size_t, size_t> duplicate_map;
    if(geometries.size() > 1){
        duplicate_map = this->find_duplicates();
    }

    size_t nodes_per_elem = this->elem_info->get_nodes_per_element();

    auto list = this->generate_element_shapes(elem_tags, elem_node_tags, nodes_per_elem, duplicate_map);

    this->prune(list);

    this->prepare_for_FEM(shape, list, forces, supports);
}

void Gmsh::gmsh_meshing(bool has_condition_inside, TopoDS_Shape sh, std::vector<size_t>& elem_tags, std::vector<size_t>& elem_node_tags, const MeshElementFactory* const elem_type){
    gmsh::initialize();

    gmsh::model::add("base");

    gmsh::vectorpair vec;
    gmsh::model::occ::importShapesNativePointer(static_cast<const void*>(&sh), vec);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMin", this->size);
    gmsh::option::setNumber("Mesh.MeshSizeMax", this->size);

    gmsh::option::setNumber("Mesh.Algorithm", this->algorithm);

    gmsh::option::setNumber("Mesh.ElementOrder", this->elem_info->get_element_order());
    gmsh::option::setNumber("Mesh.HighOrderOptimize", 2);
    gmsh::option::setNumber("Mesh.Optimize", 1);
    gmsh::option::setNumber("Mesh.OptimizeNetgen", 1);

    size_t dim = 0;
    auto problem_type = this->elem_info->get_problem_type();
    if(problem_type == utils::PROBLEM_TYPE_2D){
        dim = 2;
    } else if(problem_type == utils::PROBLEM_TYPE_3D){
        dim = 3;
    }

    gmsh::model::mesh::generate(dim);

    std::vector<std::size_t> node_tags;
    std::vector<double> node_coords, node_params;
    if(has_condition_inside){
        gmsh::model::mesh::getNodes(node_tags, node_coords, node_params, -1, -1, true);
    } else {
        gmsh::model::mesh::getNodes(node_tags, node_coords, node_params, dim, -1, true);
    }

    // Would need to be changed to support multiple elements
    size_t type = elem_type->get_gmsh_element_type();
    gmsh::model::mesh::getElementsByType(type, elem_tags, elem_node_tags, -1);

    gmsh::clear();
    gmsh::finalize();

    size_t dof = elem_type->get_dof_per_node();

    this->node_list.clear();
    this->node_list.reserve(node_tags.size());
    for(size_t i = 0; i < node_tags.size(); ++i){
        gp_Pnt p(node_coords[i*3], node_coords[i*3+1], node_coords[i*3+2]);
        this->node_list.emplace_back(std::make_unique<MeshNode>(p, node_tags[i], dof));
    }

}

}
