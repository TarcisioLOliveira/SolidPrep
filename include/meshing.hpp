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
#include <unordered_map>
#include "element.hpp"
#include <vector>
#include "support.hpp"
#include "force.hpp"
#include "element_factory.hpp"
#include "geometry.hpp"

class ProjectData;

class Meshing{
    public:
    Meshing(const std::vector<std::unique_ptr<Geometry>>& geometries,
            const MeshElementFactory* const elem_type,
            const double thickness):
        elem_info(elem_type), geometries(utils::extract_pointers(geometries)),
        thickness(thickness){}

    /**
     * Just generates a mesh from a group of geometries, depending on analysis
     * parameters set in their construction. Returns ElementShape instances,
     * not elements. Elements should be made using prepare_for_FEM().
     *
     * @param geoemtries The geometries to be meshed
     * @param elem_type Element type
     *
     * @return collection of ElementShapes (which are just a collection of
     *         MeshNodes
     */
    virtual void mesh(const std::vector<Force>& forces, 
                      const std::vector<Support>& supports) = 0;

    /**
     * Removes elements below a certain density threshold. Useful for
     * validating and saving results.
     *
     * @param rho Density vector
     * @param threshold Density threshold
     */
    virtual void prune(const std::vector<Force>& forces, 
                       const std::vector<Support>& supports,
                       const std::vector<double>& rho, double threshold);

    const MeshElementFactory * const elem_info;
    const std::vector<Geometry*> geometries;
    const double thickness;

    std::vector<std::unique_ptr<MeshNode>> node_list;
    std::vector<double> load_vector;

    protected:
    std::vector<long> get_support_dof(size_t& offset, size_t id, const Support& support, const MeshElementFactory* elem_maker) const;
    std::vector<double> get_force_dof(const Force& force, const MeshElementFactory* elem_maker) const;
    void reverse_cuthill_mckee(const std::vector<ElementShape>& elem_list);
    // Inside and not on boundary
    bool is_strictly_inside2D(gp_Pnt p, TopoDS_Shape s) const;

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
    virtual void prepare_for_FEM(const TopoDS_Shape& shape,
                                 const std::vector<ElementShape>& base_mesh,
                                 const std::vector<Force>& forces, 
                                 const std::vector<Support>& supports);

    /**
     * Prunes and optimizes list of nodes after generation.
     *
     * @param list List of element shapes.
     */
    virtual void prune(const std::vector<ElementShape>& list);

    /**
     * Generates the list of element shapes from mesh information.
     *
     * @param elem_tags Element tags
     * @param elem_node_tags Node tags per element
     * @param nodes_per_elem Number of nodes in each element
     * @param duplicate_map Map of duplicate nodes (optional)
     *
     * @return Vector of element shapes
     */
    virtual std::vector<ElementShape> generate_element_shapes(const std::vector<size_t>& elem_tags, const std::vector<size_t>& elem_node_tags, size_t nodes_per_elem,const std::unordered_map<size_t, size_t>& duplicate_map = std::unordered_map<size_t, size_t>());

    /**
     * Search for duplicate nodes (different tag but same position)
     * That may happen also for single geometries that get cut by a boundary
     * condition, but they behave correctly without this workaround.
     *
     * @return Mapping of duplicate nodes, redirecting the duplicates to single
     * one.
     */
    virtual std::unordered_map<size_t, size_t> find_duplicates();

    /**
     * Turns a vector of geometries into a single compound geometry.
     *
     * @param geometries Vector of Geometry instances
     *
     * @return OCCT compound as TopoDS_Shape
     */
    virtual TopoDS_Shape make_compound(const std::vector<Geometry*>& geometries) const;

    virtual bool adapt_for_boundary_condition_inside(TopoDS_Shape& shape, const std::vector<Force>& forces, const std::vector<Support>& supports);
};

#endif
