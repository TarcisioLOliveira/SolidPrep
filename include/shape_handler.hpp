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
#include <limits>
#include <set>

namespace shape_op{

enum class Code{
    UNION,
    INTERSECTION,
    DIFFERENCE,
    GEOMETRY,
    SHELL
};

class ShapeOp{
    public:
    ShapeOp(std::unique_ptr<ShapeOp> s1, std::unique_ptr<ShapeOp> s2):
        s1(std::move(s1)), s2(std::move(s2)){}

    virtual ~ShapeOp() = default;

    virtual Code get_type() const = 0;
    inline ShapeOp* first() const{
        return s1.get();
    }
    inline ShapeOp* second() const{
        return s2.get();
    }
    virtual size_t get_id() const{
        return std::numeric_limits<size_t>::max();
    }

    private:
    std::unique_ptr<ShapeOp> s1, s2;
};

class Union : public ShapeOp{
    public:
    Union(std::unique_ptr<ShapeOp> s1, std::unique_ptr<ShapeOp> s2):
        ShapeOp(std::move(s1), std::move(s2)){}

    virtual ~Union() = default;
    virtual Code get_type() const override{
        return Code::UNION;
    }
};

class Intersection : public ShapeOp{
    public:
    Intersection(std::unique_ptr<ShapeOp> s1, std::unique_ptr<ShapeOp> s2):
        ShapeOp(std::move(s1), std::move(s2)){}

    virtual ~Intersection() = default;
    virtual Code get_type() const override{
        return Code::INTERSECTION;
    }
};

class Difference : public ShapeOp{
    public:
    Difference(std::unique_ptr<ShapeOp> s1, std::unique_ptr<ShapeOp> s2):
        ShapeOp(std::move(s1), std::move(s2)){}

    virtual ~Difference() = default;
    virtual Code get_type() const override{
        return Code::DIFFERENCE;
    }
};

class Geometry : public ShapeOp{
    public:
    Geometry(size_t id):
        ShapeOp(nullptr, nullptr), id(id){}

    virtual ~Geometry() = default;
    virtual Code get_type() const override{
        return Code::GEOMETRY;
    }
    virtual size_t get_id() const override{
        return this->id;
    }

    private:
    const size_t id;
};

class Shell : public ShapeOp{
    public:
    Shell(size_t id):
        ShapeOp(nullptr, nullptr), id(id){}

    virtual ~Shell() = default;
    virtual Code get_type() const override{
        return Code::SHELL;
    }
    virtual size_t get_id() const override{
        return this->id;
    }

    private:
    const size_t id;
};


};

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

    struct SuperimposedNodes{
        public:
        size_t id;
        std::vector<MeshNode*> nodes;
    };

    ShapeHandler(Meshing* mesh, std::vector<Geometry*> geometries, std::unique_ptr<shape_op::ShapeOp> root_op);

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
    std::set<SuperimposedNodes*> apply_op(shape_op::ShapeOp* op) const;

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

    std::map<size_t, SuperimposedNodes*> merged_nodes_mapping;
    std::vector<SuperimposedNodes> merged_nodes;

    std::vector<double> original_points;
    std::vector<double> shape_displacement;

    std::unique_ptr<general_solver::MUMPSGeneral> solver;
    std::vector<double> b;

    std::unique_ptr<shape_op::ShapeOp> root_op;
};

#endif
