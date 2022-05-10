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

#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <cmath>
#include <gp_Pnt.hxx>
#include <vector>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Edge.hxx>

/**
 * Basic node superclass.
 */
class Node{
    public:
    const gp_Pnt point; ///< Node position
    size_t id;
    double * const results; ///< Stores results, such as stress or reactions
    /**
     * Stores the position of the values of local displacement in the global
     * displacement vector.
     */
    long * u_pos;
    ~Node(){ delete[] results; }

    protected:
    /**
     * Constructs the node instance.
     *
     * @param p Node position.
     * @param id Node id.
     * @param dim Size of results array.
     */
    Node(gp_Pnt p, size_t id, size_t dim):point(p), id(id), results(new double[dim]()), u_pos(nullptr){}
};

/**
 * Basic superclass for 1D beam nodes.
 * Used only in BeamGraph.
 */
class BeamNode : public Node{
    public:
    const double dim;
    const gp_Dir normal;
    virtual ~BeamNode() = default;
    protected:
    /**
     * Constructs the beam node instance. 
     *
     * @param p Position.
     * @param id Node id.
     * @param res_n Size of results array.
     * @param dim Height of the cross-section at point p.
     * @param normal Normal of the cross-section at point p.
     */
    BeamNode(gp_Pnt p, size_t id, size_t res_n, double dim, gp_Dir n):Node(p, id, res_n),dim(dim), normal(n){}
};

/**
 * Basic superclass for 1D beam nodes in 2D environments.
 * Used only in BeamGraph.
 */
class BeamNode2D : public BeamNode{
    public:
    /**
     * Constructs the beam node instance.
     *
     * @param p Position.
     * @param id Node id.
     * @param dim Height of the cross-section at point p.
     * @param n Normal of the cross-section at point p.
     */
    BeamNode2D(gp_Pnt p, size_t id, double dim, gp_Dir n):BeamNode(p, id, 3, dim, n){}
};

/**
 * Class used to create instances of the desired BeamNode subclass more easily.
 */
class BeamNodeFactory{
    public:
    enum BeamNodeType{
        NONE,
        BEAM_NODE_2D
    };
    /**
     * Constructs the beam node instance of the desired element type.
     *
     * @param p Position.
     * @param id Node id.
     * @param dim Height of the cross-section at point p.
     * @param n Normal of the cross-section at point p.
     * @param t Element type.
     */
    static BeamNode* make_node(gp_Pnt p, size_t id, double dim, gp_Dir n, BeamNodeType t);
};

/**
 * Basic superclass for 2D and 3D elements.
 */
class MeshNode : public Node{
    public:
    virtual ~MeshNode() = default;
    /**
     * Gets the size of the results array.
     *
     * @return Size of results array.
     */
    virtual size_t get_result_size() const = 0;
    /**
     * Gets Von Mises stress based on results array.
     *
     * @return Von Mises stress.
     */
    virtual double get_Von_Mises() const = 0;
    protected:
    /**
     * Constructs the mesh node instance.
     *
     * @param p Position.
     * @param id Node id.
     * @param res_n Size of the result array.
     */
    MeshNode(gp_Pnt p, size_t id, size_t res_n): Node(p, id, res_n){}
};

/**
 * Superclass for 2D elements.
 *
 * @see MeshNode
 */
class MeshNode2D : public MeshNode{
    public:
    /**
     * Constructs the 2D mesh node instance.
     *
     * @param p Position.
     * @param id Node id.
     */
    MeshNode2D(gp_Pnt p, size_t id):MeshNode(p, id, 3){}
    virtual size_t get_result_size() const override{ return 3; }
    virtual double get_Von_Mises() const override{
        return std::sqrt(std::pow(results[0], 2) - results[0]*results[1] + std::pow(results[1], 2) + 3*std::pow(results[2], 2));
    }
};

/**
 * Class used to create instances of the desired MeshNode subclass more easily.
 */
class MeshNodeFactory{
    public:
    enum MeshNodeType{
        NONE,
        MESH_NODE_2D
    };
    /**
     * Constructs the 2D mesh node instance with the specified type.
     *
     * @param p Position.
     * @param id Node id.
     * @param t Type of node.
     */
    static MeshNode* make_node(gp_Pnt p, size_t id, MeshNodeType t);
};

/**
 * Struct used to store the nodes that compose an element. Used between the
 * meshing process and the element generation.
 *
 * @see Meshing
 */
struct ElementShape{
    std::vector<MeshNode*> nodes;
};

/**
 * Basic superclass for finite elements.
 */
class Element{
    public:

    std::vector<Node*> nodes;

    virtual ~Element() = default;
    /**
     * Creates and returns the elemental stiffness matrix.
     *
     * @return The element stiffness matrix.
     */
    virtual std::vector<double> get_k() const = 0;

    protected:
    /**
     * Creates the element based on a predefined set of nodes.
     *
     * @param n Nodes.
     */
    Element(std::vector<Node*> n):nodes(std::move(n)){}
};

/**
 * Basic superclass for beam elements. Used by BeamSizing and BeamGraph only,
 * inherited by the BeamLinear2D element.
 *
 * @see BeamSizing
 * @see BeamGraph
 * @see BeamLinear2D
 */
class BeamElement : public Element{
    public:

    virtual ~BeamElement() = default;
    /**
     * Creates and returns the elemental stiffness matrix.
     *
     * @return The element stiffness matrix.
     */
    virtual std::vector<double> get_k() const override = 0;
    /**
     * Calculates the internal loads at a certain element node.
     *
     * @param node Number of the node within the element's list of nodes.
     * @param u Displacement vector.
     *
     * @return A pointer to the node.
     */
    virtual BeamNode* get_internal_loads(size_t node, const std::vector<double>& u) const = 0;
    /**
     * Returns an element node.
     *
     * @param node Number of the node within the element's list of nodes.
     *
     * @return A pointer to the node.
     */
    virtual inline BeamNode* get_node(size_t node) const{ return static_cast<BeamNode*>(this->nodes[node]);}

    protected:
    /**
     * Creates a beam element using two beam nodes.
     *
     * @param p1 Beam node 1.
     * @param p2 Beam node 2.
     */
    BeamElement(BeamNode* p1, BeamNode* p2):Element({p1,p2}){}
};

/**
 * Basic superclass to define 2D and 3D elements.
 */
class MeshElement : public Element{
    public:

    virtual ~MeshElement() = default;
    /**
     * Creates and returns the elemental stiffness matrix.
     *
     * @return The element stiffness matrix.
     */
    virtual std::vector<double> get_k() const override = 0;
    /**
     * Calculates the internal load vector at a node.
     *
     * @param node Number of the node within the element's list of nodes.
     * @param u Displacement vector.
     *
     * @return A pointer to the node.
     */
    virtual MeshNode* get_internal_loads(size_t node, const std::vector<double>& u) const = 0;
    /** 
     * Calculates the Von Mises stresses at a point within the element.
     *
     * @param p The point.
     * @param u Displacement vector.
     *
     * @return Von Mises stress at point p.
     */
    virtual double get_stress_at(gp_Pnt p, const std::vector<double>& u) const = 0;
    /**
     * Calculates the Cauchy stress tensor at a point within the element.
     *
     * @param p The point.
     * @param u Displacement vector.
     *
     * @return Stress tensor at point p.
     */
    virtual std::vector<double> get_stress_tensor(gp_Pnt p, const std::vector<double>& u) const = 0;
    /** 
     * Calculates the internal loads at a point within the element. Does not
     * consider contribution from other elements.
     *
     * @param p The point.
     * @param u Displacement vector.
     *
     * @return Internal loads.
     */
    virtual std::vector<double> get_loads_at(gp_Pnt p, const std::vector<double>& u) const = 0;
    // Edge for 2D, face for 3D
    /** 
     * Calculates the intersection points between a shape (edge for 2D, face
     * for 3D) and the boundaries of the element.
     *
     * @param crosssection Intersecting shape.
     *
     * @return List of intersected points.
     */
    virtual std::vector<gp_Pnt> get_intersection_points(const TopoDS_Shape& crosssection) const = 0;
    /**
     * Returns the volume or area of the element.
     *
     * @return volume or area.
     */
    virtual double get_volume() const = 0;
    /**
     * Calculates the element's compliance.
     *
     * @param u Displacement vector.
     * @param l Virtual displacement vector (optional).
     *
     * @return The compliance.
     */
    virtual double get_compliance(const std::vector<double>& u, const std::vector<double>& l = std::vector<double>()) const = 0;
    /**
     * Calculates the elemental contribution to the virtual load, for a p-norm
     * global stress aggregation function.
     *
     * @param mult Value to be multiplied to the result. For p-norm stress 
     * based on LE et al, 2009, it is:
     * P*v*std::pow(data->new_x[i], pt*P)*std::pow(S, P-2)
     * @param u Displacement vector.
     * @param l Global virtual load vector.
     */
    virtual void get_virtual_load(double mult, gp_Pnt point, const std::vector<double>& u, std::vector<double>& l) const = 0;
    /**
     * Calculates the force vector for a distributed load.
     *
     * @param dir Load direction.
     * @param norm Distributed load intensity over area (that is, pressure).
     * @param points Delimiting points.
     * @return Force vector.
     */
    virtual std::vector<double> get_f(gp_Dir dir, double norm, std::vector<gp_Pnt> points) const = 0;
    /**
     * Returns the geometry of the element.
     *
     * @param disp Displacement of nodes (optional).
     *
     * @return Shape of the element.
     */
    virtual TopoDS_Shape get_shape(std::vector<gp_Vec> disp = std::vector<gp_Vec>()) const = 0;
    /**
     * Calculates the centroid of the element.
     *
     * @return The centroid.
     */
    virtual gp_Pnt get_centroid() const = 0;

    /**
     * Necessary to use Gmsh's GUI. See:
     * https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
     * 
     * @return Gmsh element type number.
     */
    virtual size_t get_gmsh_element_type() const = 0;

    /**
     * Returns an element node.
     *
     * @param node Number of the node within the element's list of nodes.
     *
     * @return A pointer to the node.
     */
    virtual inline MeshNode* get_node(size_t node) const{ return static_cast<MeshNode*>(this->nodes[node]);}

    protected:
    /**
     * Creates an element with the specified nodes.
     *
     * @param nodes List of nodes.
     */
    MeshElement(std::vector<MeshNode*> nodes):Element(std::vector<Node*>(nodes.begin(), nodes.end())){}
};

#endif
