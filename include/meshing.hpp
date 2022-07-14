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
    Meshing(double size):size(size){}

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
     * The object stores the results.
     *
     * @param base_mesh Mesh obtained from mesh()
     * @param element_type Type of element to be used
     * @param supports List of supports
     * @param force List of forces
     */
    virtual void prepare_for_FEM(const std::vector<ElementShape>& base_mesh,
                                 MeshElementFactory::MeshElementType element_type,
                                 ProjectData* data, bool force_only = false);

    /**
     * Removes elements below a certain density threshold. Useful for
     * validating and saving results.
     *
     * @param rho Density vector
     * @param threshold Density threshold
     */
    virtual std::vector<ElementShape> prune(const std::vector<double>& rho, double threshold);

    std::vector<std::unique_ptr<MeshNode>> node_list;
    std::vector<std::unique_ptr<MeshElement>> element_list;
    std::vector<double> load_vector;
    TopoDS_Shape shape;

    protected:
    double size;
    MeshElementFactory::MeshElementType type;
    std::vector<long> get_support_dof(size_t& offset, size_t id, const Support& support, MeshElementFactory::MeshElementType type) const;
    std::vector<double> get_force_dof(const Force& force, MeshElementFactory::MeshElementType type) const;
    void reverse_cuthill_mckee(const std::vector<ElementShape>& elem_list);
    // Inside and not on boundary
    bool is_strictly_inside2D(gp_Pnt p, TopoDS_Shape s) const;
};

#endif
