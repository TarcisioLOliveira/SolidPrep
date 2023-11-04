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
#include "spring.hpp"
#include "element_factory.hpp"
#include "geometry.hpp"

class ProjectData;

class BoundaryElement{
    public:
    const MeshNode** const nodes;
    const MeshElement* const parent;
    const gp_Dir normal;
    
    ~BoundaryElement(){
        delete[] nodes;
    }
    BoundaryElement(const std::vector<MeshNode*>& n, const MeshElement* const parent, gp_Dir normal):
        nodes(allocate_nodes(n)), parent(parent), normal(std::move(normal)){}

    inline gp_Pnt get_centroid(const size_t N) const{
        double x = 0, y = 0, z = 0;
        for(size_t i = 0; i < N; ++i){
            x += this->nodes[i]->point.X();
            y += this->nodes[i]->point.Y();
            z += this->nodes[i]->point.Z();
        }
        return gp_Pnt(x/N, y/N, z/N);
    }

    private:
    /**
     * Sets the elements from the element vector.
     *
     * @param n Nodes
     *
     * @return Newly allocated matrix containing pointer to node array.
     */
    const MeshNode** allocate_nodes(const std::vector<MeshNode*>& n){
        const MeshNode** nodes(new const MeshNode*[n.size()]);
        std::copy(n.begin(), n.end(), nodes);
        return nodes;
    }
};

class Meshing{
    public:
    Meshing(const std::vector<std::unique_ptr<Geometry>>& geometries,
            const MeshElementFactory* const elem_type,
            const double thickness):
        elem_info(elem_type), geometries(utils::extract_pointers(geometries)),
        thickness(thickness){}

    virtual ~Meshing() = default;

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
                      const std::vector<Support>& supports,
                      std::vector<Spring>& springs) = 0;

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
    std::unordered_multimap<size_t, MeshElement*> inverse_mesh;
    std::vector<BoundaryElement> boundary_elements;
    std::vector<MeshNode*> boundary_node_list;
    std::vector<Spring>* springs;

    protected:
    std::vector<bool> get_support_dof(const Support& support, const MeshElementFactory* elem_maker) const;
    std::vector<double> get_force_dof(const Force& force, const MeshElementFactory* elem_maker) const;
    void reverse_cuthill_mckee(const std::vector<ElementShape>& elem_list);
    // Inside and not on boundary
    bool is_strictly_inside2D(gp_Pnt p, TopoDS_Shape s) const;
    bool is_strictly_inside3D(gp_Pnt p, TopoDS_Shape s) const;

    /**
     * Main function for populating each geometry's element list. Uses the other
     * protected functions in this class to do so. Can be overriden if desired,
     * as can the other functions.
     *
     * @param shape Shape being meshed
     * @param geom_elem_mapping Maps element ranges to the geometries which should contain them
     * @param elem_node_tags Map each node to its elements
     * @param bound_elem_node_tags Maps each boundary node to boundary elements
     * @param id_map Maps the original node tags (e.g. from Gmsh) to their MeshNode instance
     * @param forces Forces to be used
     * @param supports Supports to be used
     * @param springs Springs to be used
     * @param deduplicate Whether to deduplicate nodes
     * @param boundary_condition_inside Whether there is one or more boundary conditions inside the geometry
     */
    virtual void generate_elements(const TopoDS_Shape& shape,
                                   const std::vector<size_t>& geom_elem_mapping, 
                                   const std::vector<size_t>& elem_node_tags, 
                                   const std::vector<size_t>& bound_elem_node_tags,
                                   std::unordered_map<size_t, MeshNode*>& id_map,
                                   const std::vector<Force>& forces, 
                                   const std::vector<Support>& supports,
                                   std::vector<Spring>& springs,
                                   const bool deduplicate,
                                   const bool boundary_condition_inside);

    /**
     * Creates a list of elements from the list of element shapes.
     *
     * @param base_mesh List of ElementShape
     * @param elem_info Element type
     *
     * @return Vector of MeshElement instances
     */
    virtual std::vector<std::unique_ptr<MeshElement>> create_element_list(
                                   const std::vector<ElementShape>& base_mesh, 
                                   const MeshElementFactory * const elem_info) const;

    /**
     * Populates the inverse_mesh member variable.
     *
     * @param element_list List of MeshElement instances (created e.g. using create_element_list())
     */
    virtual void populate_inverse_mesh(const std::vector<std::unique_ptr<MeshElement>>& element_list);

    /**
     * Populates the boundary_elements member variable.
     *
     * @param boundary_base_mesh List of shapes of boundary elements
     * @param boundary_condition_inside Whether there is one or more boundary conditions inside the geometry
     */
    virtual void populate_boundary_elements(const std::vector<ElementShape>& boundary_base_mesh,
                                            const bool boundary_condition_inside);

    /**
     * Generate positioning information for node degrees of freedom, including
     * which DOFs should be removed.
     *
     * @param supports List of supports to be applied.
     */
    void apply_supports(const std::vector<Support>& supports);

    /**
     * Generate positioning information for spring boundary condition.
     *
     * @param springs List of springs to be applied.
     */
    void apply_springs(std::vector<Spring>& springs);

    /**
     * Generates the load vector member variable based on the mesh, shape and
     * loads applied.
     *
     * @param shape Shape being meshed
     * @param forces Loads being applied
     */
    void generate_load_vector(const TopoDS_Shape& shape,
                              const std::vector<Force>& forces);


    /**
     * Distributes elements to their respective geometries.
     *
     * @param element_list List of elements
     */
    void distribute_elements(const std::vector<size_t>& geom_elem_mapping, 
                             std::vector<std::unique_ptr<MeshElement>>& element_list);

    /**
     * DEPRECATED. Only used by StandardBeamMesher (which is also being
     * deprecated).
     *
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
                                 const std::vector<size_t>& geom_elem_mapping,
                                 const std::vector<ElementShape>& base_mesh,
                                 const std::vector<Force>& forces, 
                                 const std::vector<Support>& supports);

    /**
     * Prunes and optimizes list of nodes after generation.
     *
     * @param list List of element shapes.
     */
    virtual void optimize(std::vector<ElementShape>& list, std::unordered_map<size_t, MeshNode*>* id_map = nullptr);

    /**
     * Generates the list of element shapes from mesh information.
     *
     * @param elem_node_tags Node tags per element
     * @param nodes_per_elem Number of nodes in each element
     * @param id_map Map of original ids to nodes
     * @param calc_normals Whether to calculate element normals (for boundary elements)
     *
     * @return Vector of element shapes
     */
    virtual std::vector<ElementShape> generate_element_shapes(
            const std::vector<size_t>& elem_node_tags, 
            size_t bound_nodes_per_elem,
            const std::unordered_map<size_t, MeshNode*>& id_map,
            bool calc_normals = false);

    /**
     * Search for duplicate nodes (different tag but same position)
     * That may happen also for single geometries that get cut by a boundary
     * condition, but they behave correctly without this workaround.
     *
     * @param id_map Mapping of node tags to node_list indices.
     *
     * @return Indices of duplicated nodes, to be used late in remove_duplicates().
     */
    virtual std::vector<size_t> find_duplicates(std::unordered_map<size_t, MeshNode*>& id_map);

    /**
     * Turns a vector of geometries into a single compound geometry.
     *
     * @param geometries Vector of Geometry instances
     * @param scale Scaling constant
     *
     * @return OCCT compound as TopoDS_Shape
     */
    virtual TopoDS_Shape make_compound(const std::vector<Geometry*>& geometries, const double scale = 1.0) const;

    virtual bool adapt_for_boundary_condition_inside(TopoDS_Shape& shape, const std::vector<Force>& forces, const std::vector<Support>& supports);
};

#endif
