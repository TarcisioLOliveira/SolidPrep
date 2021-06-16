/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   Foobar is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Foobar is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with SolidPrep.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef BEAM_GRAPH_HPP
#define BEAM_GRAPH_HPP

#include "project_data.hpp"
#include "element.hpp"
#include "element_factory.hpp"
#include "lapacke.h"

// It's not really a mesh, so...
class BeamGraph{
    public:
    BeamGraph(ProjectData* data, BeamElementFactory::BeamElementType t):data(data), type(t){}
    ~BeamGraph();

    void run();
    
    inline BeamNode* get(size_t id){
        return this->nodes[id];
    }
    inline size_t size(){
        return this->nodes.size();
    }
    private:
    ProjectData * const data;
    std::vector<BeamNode*> nodes;
    BeamElementFactory::BeamElementType type;
    void insert_element_matrix(std::vector<float>& K, const std::vector<float>& k, const std::vector<long>& pos, int w, int& n) const;
    void clear_nodes();
};

#endif
