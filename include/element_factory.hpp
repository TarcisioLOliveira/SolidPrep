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

#ifndef ELEMENT_FACTORY_HPP
#define ELEMENT_FACTORY_HPP

#include "element.hpp"
#include "element/beam_linear_2D.hpp"

class BeamElementFactory{
    public:
    enum BeamElementType{
        NONE,
        BEAM_LINEAR_2D
    };
    template<typename ... Args>
    static BeamElement* make_element(BeamElementType t, Args&& ... args){
        if(t == BEAM_LINEAR_2D){
            return new element::BeamLinear2D(args...);
        }

        return nullptr;
    }
    static BeamNodeFactory::BeamNodeType get_node_type(BeamElementType t){
        if(t == BEAM_LINEAR_2D){
            return BeamNodeFactory::BEAM_NODE_2D;
        }
        return BeamNodeFactory::NONE;
    }
    static size_t get_k_dimension(BeamElementType t){
        if(t == BEAM_LINEAR_2D){
            return 6;
        }
        return 0;
    }
};

#endif
