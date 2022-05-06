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

#ifndef BEAM_MESHING_HPP
#define BEAM_MESHING_HPP

#include "meshing.hpp"

/**
 * Interface class for meshers designed to be used on continuum topologies by
 * assuming that they are composed of multiple intersecting beams, in order to
 * be used by beam sizing algorithms.
 *
 * @see Meshing
 * @see BeamSizing
 */
class BeamMeshing : public Meshing{
    public:
    /**
     * Represents a node situated on the topology's boundary.
     */
    struct BoundaryNode{
        BoundaryNode(MeshNode* n):node(n){}
        MeshNode* node; ///< Associated mesh node
        gp_Dir normal; ///< Normal to the node, pointing outward from the topology
    };

    /**
     * Constructs an object by predetermining element size.
     *
     * @param size Element size.
     */
    BeamMeshing(double size):Meshing(size){}

    std::vector<BoundaryNode> boundary_nodes;
};

#endif
