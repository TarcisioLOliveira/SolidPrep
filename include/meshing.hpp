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

#ifndef MESHING_HPP
#define MESHING_HPP

#include <TopoDS_Shape.hxx>
#include "element.hpp"
#include <vector>
#include "support.hpp"
#include "force.hpp"
#include "element_factory.hpp"

class ProjectData;

class Meshing{
    public:

    /**
     * Just generates a mesh from a shape according to the child class's input
     * parameters set in its construction. Return MeshNodes are just pointers
     * to the ones stored in the class object.
     *
     * @param s OCCT shape to be meshed
     * @return collection of ElementShapes (which are just a collection of
     *         MeshNodes
     */
    virtual std::vector<ElementShape> mesh(TopoDS_Shape s) = 0;
    /**
     * Takes a mesh and boundary conditions and returns a collection of
     * elements and its load vector, to be used later in a FiniteElement
     * object.
     *
     * Assumes you are using the same instace which generated base_mesh (for
     * simplicity).
     *
     * As it returns two vector instances, both are returned as function
     * arguments (final_mesh and load_vector).
     *
     * The object also stores the mesh elements created.
     *
     * @param base_mesh Mesh obtained from mesh()
     * @param element_type Type of element to be used
     * @param supports List of supports
     * @param force List of forces
     * @param final_mesh Output element list
     * @param load_vector Output external force vector
     */
    virtual void prepare_for_FEM(const std::vector<ElementShape>& base_mesh,
                                 MeshElementFactory::MeshElementType element_type,
                                 ProjectData* data,
                                 std::vector<MeshElement*>& final_mesh,
                                 std::vector<double>& load_vector);

    protected:
    std::vector<std::unique_ptr<MeshNode>> node_list;
    std::vector<std::unique_ptr<MeshElement>> element_list;

    std::vector<long> get_support_dof(size_t& offset, size_t id, const Support& support, MeshElementFactory::MeshElementType type) const;
    std::vector<double> get_force_dof(const Force& force, MeshElementFactory::MeshElementType type) const;
};

#endif
