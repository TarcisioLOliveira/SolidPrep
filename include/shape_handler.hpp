/*
 *   Copyright (C) 2024 Tarcísio Ladeia de Oliveira.
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

#ifndef SHAPE_HANDLER_HPP
#define SHAPE_HANDLER_HPP

#include "meshing.hpp"
#include "general_solver/mumps_general.hpp"

class ShapeHandler{
    public:
    struct AffectedElement{
        public:
        const MeshElement* e;
        const size_t node_num;
    };
    struct AffectedNode{
        public:
        const std::vector<size_t> node_ids;
        const std::vector<AffectedElement> elements;
    };

    ShapeHandler(Meshing* mesh, std::vector<Geometry*> geometries);

    void obtain_affected_nodes();
    void update_nodes(const std::vector<double>& dx);

    inline const std::vector<double>& get_shape_displacement() const{
        return this->shape_displacement;
    }
    inline const std::vector<AffectedNode>& get_nodes() const{
        return this->optimized_nodes;
    }
    inline size_t get_number_of_nodes() const{
        return this->optimized_nodes.size();
    }
    size_t get_number_of_variables() const{
        return this->optimized_nodes.size()*this->mesh->elem_info->get_dof_per_node();
    }

    private:
    Meshing* mesh;
    std::vector<Geometry*> geometries;

    std::vector<AffectedNode> optimized_nodes;
    std::vector<std::unique_ptr<ShapeMeshElement>> shape_elements;
    std::vector<BoundaryElement*> boundary_elements;
    std::map<size_t, long> id_mapping;
    std::map<size_t, size_t> bound_to_shape_mapping;
    std::map<size_t, size_t> optimized_nodes_mapping;
    std::map<size_t, MeshElement*> node_to_elem_unique_mapping;
    std::map<MeshElement*, std::vector<size_t>> elem_to_affected_node_mapping;
    std::vector<MeshNode*> domain_nodes;
    size_t matrix_width;
    // Full boundary except for boundary conditions
    bool full_boundary_optimization = true;

    std::vector<double> original_points;
    std::vector<double> shape_displacement;

    std::unique_ptr<general_solver::MUMPSGeneral> solver;
    std::vector<double> b;
};

#endif
