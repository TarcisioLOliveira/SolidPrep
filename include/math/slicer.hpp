/*
 *   Copyright (C) 2025 Tarcísio Ladeia de Oliveira.
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

#ifndef MATH_SLICER_HPP
#define MATH_SLICER_HPP

#include <vector>
#include "element.hpp"

namespace math{

namespace slicer{

template<typename T>
inline std::vector<T> from_node_id(const Node* const* nodes, size_t number_of_nodes, size_t dof){
    std::vector<size_t> ids(number_of_nodes*dof);
    for(size_t i = 0; i < number_of_nodes; ++i){
        const auto id = nodes[i]->id;
        for(size_t j = 0; j < dof; ++j){
            ids[i*dof + j] = id*dof + j;
        }
    }

    return ids;
}
inline std::vector<long> from_node_upos(const Node* const* nodes, size_t number_of_nodes, size_t dof){
    std::vector<long> ids(number_of_nodes*dof);
    for(size_t i = 0; i < number_of_nodes; ++i){
        const auto n = nodes[i];
        for(size_t j = 0; j < dof; ++j){
            ids[i*dof + j] = n->u_pos[j];
        }
    }

    return ids;
}

template<typename T>
void from_node_id(const Node* const* nodes, size_t number_of_nodes, size_t dof, std::vector<T>& pos){
    for(size_t i = 0; i < number_of_nodes; ++i){
        const auto id = nodes[i]->id;
        for(size_t j = 0; j < dof; ++j){
            pos[i*dof + j] = id*dof + j;
        }
    }
}
inline void from_node_upos(const Node* const* nodes, size_t number_of_nodes, size_t dof, std::vector<long>& upos){
    for(size_t i = 0; i < number_of_nodes; ++i){
        const auto n = nodes[i];
        for(size_t j = 0; j < dof; ++j){
            upos[i*dof + j] = n->u_pos[j];
        }
    }
}
template<typename T>
inline std::vector<T> sequence(const T start, const T stop){
    std::vector<T> seq(stop-start);
    std::iota(seq.begin(), seq.end(), start);
    return seq;
}


}

}

#endif
