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

#include <set>
#include "general_solver/petsc_general_pcg.hpp"
#include "meshing.hpp"
#include "general_solver/mumps_general.hpp"
#include "shape_operations.hpp"

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
    struct AffectedPairedBoundary{
        public:
        PairedBoundaryElements* e;
        size_t en1;
        size_t en2;
        size_t bn1;
        size_t bn2;
    };
    struct AffectedContactNode{
        public:
        const size_t id;
        const std::vector<AffectedPairedBoundary> elements;
    };

    struct SuperimposedNodes{
        public:
        size_t id;
        std::vector<MeshNode*> nodes;
    };

    ShapeHandler() = default;
    ShapeHandler(Meshing* mesh, std::vector<Geometry*> geometries, std::unique_ptr<shape_op::ShapeOp> root_op, double smoothing_radius);

    void obtain_affected_nodes();
    void update_nodes(const std::vector<double>& dx);
    void filter_displacement(std::vector<double>& dx);
    void filter_gradient(std::vector<double>& df);

    inline const std::vector<double>& get_shape_displacement() const{
        return this->shape_displacement;
    }
    inline const std::vector<AffectedNode>& get_nodes() const{
        return this->optimized_nodes;
    }
    inline const std::vector<AffectedContactNode>& get_contact_nodes() const{
        return this->optimized_contact_nodes;
    }
    inline size_t get_number_of_nodes() const{
        return this->optimized_nodes.size();
    }
    size_t get_number_of_variables() const{
        return this->optimized_nodes.size()*this->mesh->elem_info->get_dof_per_node();
    }

    private:
    struct GeometryCluster{
        std::set<Geometry*> geometries;
        std::map<size_t, long> id_mapping;
        size_t matrix_width;
        std::unique_ptr<general_solver::PETScGeneralPCG> solver;
        std::vector<double> b;
    };
    std::set<SuperimposedNodes*> apply_op(shape_op::ShapeOp* op) const;

    std::unique_ptr<general_solver::MUMPSGeneral> helmholtz_solver;

    void generate_helmholtz();

    Meshing* mesh;
    std::vector<Geometry*> geometries;

    std::map<size_t, size_t> optimized_nodes_mapping;
    std::map<const MeshElement*, std::vector<size_t>> elem_to_affected_node_mapping;
    std::map<size_t, size_t> bound_to_shape_mapping;

    std::vector<AffectedNode> optimized_nodes;
    std::vector<AffectedContactNode> optimized_contact_nodes;
    std::vector<BoundaryElement*> boundary_elements;
    std::vector<GeometryCluster> clusters;
    std::vector<ShapeMeshElement*> shape_mesh;

    std::map<size_t, SuperimposedNodes*> merged_nodes_mapping;
    std::vector<SuperimposedNodes> merged_nodes;

    std::vector<double> original_points;
    std::vector<double> shape_displacement;

    std::unique_ptr<shape_op::ShapeOp> root_op;

    double r = 1;

    constexpr static double MU = 1e6;
};

#endif
