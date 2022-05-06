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

#include "element.hpp"


BeamNode* BeamNodeFactory::make_node(gp_Pnt p, size_t id, double dim, gp_Dir n, BeamNodeType t){
    if(t == BEAM_NODE_2D){
        return new BeamNode2D(p, id, dim, n);
    }
    return nullptr;
}

MeshNode* MeshNodeFactory::make_node(gp_Pnt p, size_t id, MeshNodeType t){
    if(t == MESH_NODE_2D){
        return new MeshNode2D(p, id);
    }
    return nullptr;
}
