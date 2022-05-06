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

#ifndef BEAM_GRAPH_HPP
#define BEAM_GRAPH_HPP

#include "project_data.hpp"
#include "element.hpp"
#include "element_factory.hpp"
#include "lapacke.h"

/**
 * Generates and stores the "mesh" of a beam or group of beams, being also 
 * responsible for generating the beams themselves by using the pathfinding 
 * algorithm chosen by the project. It is not exactly a mesh, as it uses only
 * 2D elements, hence the name "graph".
 *
 * Used by the BeamSizing class only.
 * @see BeamSizing
 * @see BeamLinear2D
 *
 * @deprecated doing not support multiple beams nor any kind of intersections
 * of beams.
 */
class BeamGraph{
    public:

    /**
     * Constructs the object by storing the parameters to be used.
     *
     * @param data Pointer to the project's specifications.
     * @param t Type of beam element to be used.
     */
    BeamGraph(ProjectData* data, BeamElementFactory::BeamElementType t):data(data), type(t){}
    ~BeamGraph();

    /**
     * Generates the graph.
     */
    void run();
   
    /**
     * Gets a node based on its id.
     *
     * @param id Node's id.
     *
     * @return A pointer to the node.
     */ 
    inline BeamNode* get(size_t id){
        return this->nodes[id];
    }
    
    /**
     * Returns the amount of nodes generated.
     *
     * @return Number of nodes.
     */
    inline size_t size(){
        return this->nodes.size();
    }
    private:
    ProjectData * const data;
    std::vector<BeamNode*> nodes;
    BeamElementFactory::BeamElementType type;

    /**
     * Inserts an elemental stiffness matrix into the global stiffness matrix.
     * Lower band matrices only (both K and k). Updates height of K (var n) if
     * K needed to be resized.
     *
     * @param K Global stiffness matrix.
     * @param k elemental stiffness matrix.
     * @param pos Positioning of each element (considering they are symmetric
     * matrices).
     * @param w Length of K
     * @param n Current height of K (as band matrix).
     *
     */
    void insert_element_matrix(std::vector<double>& K, const std::vector<double>& k, const std::vector<long>& pos, int w, int& n) const;

    /**
     * Erases all nodes.
     */
    void clear_nodes();
};

#endif
