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
    const std::vector<AffectedNode>& nodes = this->optimized_nodes;

    ShapeHandler(Meshing* mesh, std::vector<Geometry*> geometries);

    void obtain_affected_nodes();
    void update_nodes(const std::vector<double>& dx);

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
};

#endif
