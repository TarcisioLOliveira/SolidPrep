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

#ifndef ELEMENT_FACTORY_HPP
#define ELEMENT_FACTORY_HPP

#include "element.hpp"
#include "element/beam_linear_2D.hpp"
#include "utils.hpp"
#include "spview.hpp"
#include <vector>

class ProjectData;

class BeamElementFactory{
    public:
    enum BeamElementType{
        NONE,
        BEAM_LINEAR_2D
    };
    template<typename ... Args>
    static BeamElement* make_element(BeamElementType t, Args&& ... args){
        switch(t){
            case BEAM_LINEAR_2D: return new element::BeamLinear2D(args...);
            case NONE: return nullptr;
        }

        return nullptr;
    }
    static BeamNodeFactory::BeamNodeType get_node_type(BeamElementType t){
        switch(t){
            case BEAM_LINEAR_2D: return BeamNodeFactory::BEAM_NODE_2D;
            case NONE: return BeamNodeFactory::NONE;
        }
        return BeamNodeFactory::NONE;
    }
    static size_t get_k_dimension(BeamElementType t){
        switch(t){
            case BEAM_LINEAR_2D: return 6;
            case NONE: return 0;
        }
        return 0;
    }
};

class MeshElementFactory{
    public:
    virtual ~MeshElementFactory() = default;

    inline virtual MeshElement* make_element(const ElementShape& shape) const = 0;
    inline virtual size_t get_k_dimension() const = 0;
    inline virtual size_t get_dof_per_node() const = 0;
    inline virtual size_t get_nodes_per_element() const = 0;
    /**
     * Necessary to use Gmsh's GUI. See:
     * https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
     * 
     * @return Gmsh element type number.
     */
    inline virtual size_t get_gmsh_element_type() const = 0;
    inline virtual size_t get_element_order() const = 0;
    inline virtual utils::ProblemType get_problem_type() const = 0;
    inline virtual Element::Shape get_shape_type() const = 0;
    inline virtual spview::defs::ElementType get_spview_code() const = 0;

    // Boundary information
    inline virtual size_t get_boundary_nodes_per_element() const = 0;
    inline virtual size_t get_boundary_gmsh_element_type() const = 0;
};

template<class T>
class MeshElementFactoryImpl : public MeshElementFactory{
    public:
    inline MeshElement* make_element(const ElementShape& shape) const override{
        return new T(shape);
    }
    inline size_t get_k_dimension() const override{
        return T::K_DIM;
    }
    inline size_t get_dof_per_node() const override{
        return T::NODE_DOF;
    }
    inline size_t get_nodes_per_element() const override{
        return T::NODES_PER_ELEM;
    }
    /**
     * Necessary to use Gmsh's GUI. See:
     * https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
     * 
     * @return Gmsh element type number.
     */
    inline size_t get_gmsh_element_type() const override{
        return T::GMSH_TYPE;
    }
    inline utils::ProblemType get_problem_type() const override{
        return T::PROBLEM_TYPE;
    }
    inline size_t get_element_order() const override{
        return T::ORDER;
    }
    inline Element::Shape get_shape_type() const override{
        return T::SHAPE_TYPE;
    }
    inline size_t get_boundary_nodes_per_element() const override{
        return T::BOUNDARY_NODES_PER_ELEM;
    }
    inline size_t get_boundary_gmsh_element_type() const override{
        return T::BOUNDARY_GMSH_TYPE;
    }
    inline spview::defs::ElementType get_spview_code() const override{
        return T::SPVIEW_CODE;
    }
};


#endif
