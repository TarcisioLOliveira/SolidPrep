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

#include "material.hpp"
#include "math/matrix.hpp"
#include <gp_Pnt.hxx>
#include <vector>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Edge.hxx>
#include <memory>

class BoundaryElement;

/**
 * Basic node superclass.
 */
class Node{
    public:
    gp_Pnt point; ///< Node position
    size_t id;
    /**
     * Stores the position of the values of local displacement in the global
     * displacement vector.
     */
    long* u_pos;

    virtual ~Node(){
        delete[] u_pos;
    }

    Node(Node&) = delete;

    Node& operator=(Node&) = delete;

    protected:
    /**
     * Constructs the node instance.
     *
     * @param p Node position.
     * @param id Node id.
     */
    Node(gp_Pnt p, size_t id, size_t u_size):point(std::move(p)), id(id), u_pos(new long[u_size]){}
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
    BeamNode(gp_Pnt p, size_t id, double dim, gp_Dir n, size_t u_size):Node(p, id, u_size),dim(dim), normal(n){}
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
    BeamNode2D(gp_Pnt p, size_t id, double dim, gp_Dir n):BeamNode(p, id, dim, n, 3){}
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
     * Constructs the mesh node instance.
     *
     * @param p Position.
     * @param id Node id.
     */
    MeshNode(gp_Pnt p, size_t id, size_t u_size): Node(std::move(p), id, u_size){}
};

/**
 * Struct used to store the nodes that compose an element. Used between the
 * meshing process and the element generation.
 *
 * @see Meshing
 */
struct ElementShape{
    std::vector<MeshNode*> nodes;
    gp_Dir normal = gp_Dir(0,0,1);
};

/**
 * Basic superclass for finite elements.
 */
class Element{
    public:
    enum class Shape{
        TRI,
        QUAD
    };

    // Nodes sorted according to Gmsh's style.
    // Used for rendering element and other functions which require knowing
    // the geometry and its edges.
    Node** const nodes;
    size_t id = 0;

    virtual ~Element(){
        delete[] nodes;
    }

    /**
     * Creates and returns the elemental stiffness matrix.
     *
     * @return The element stiffness matrix.
     */
    // Removed just to avoid having to deal with BeamElement
    // virtual std::vector<double> get_k() const = 0;

    protected:
    /**
     * Creates the element based on a predefined set of nodes.
     *
     * @param n Nodes.
     */
    Element(const std::vector<Node*>& n):
        nodes(allocate_nodes(n)){}
    /**
     * Sets the elements from the element vector.
     *
     * @param n Nodes
     *
     * @return Newly allocated matrix containing pointer to node array.
     */
    Node** allocate_nodes(const std::vector<Node*>& n){
        Node** nodes(new Node*[n.size()]);
        std::copy(n.begin(), n.end(), nodes);
        return nodes;
    }
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
    virtual std::vector<double> get_k() const = 0;
    /**
     * Calculates the internal loads at a certain element node.
     *
     * @param node Number of the node within the element's list of nodes.
     * @param u Displacement vector.
     *
     * @return Internal loads per node.
     */
    virtual std::vector<double> get_internal_loads(size_t node, const std::vector<double>& u) const = 0;
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

class MeshElementFactory;
/**
 * Basic superclass to define 2D and 3D elements.
 */
class MeshElement : public Element{
    public:

    virtual ~MeshElement() = default;

    virtual math::Matrix get_uu(const MeshElement* const e2, const std::vector<gp_Pnt>& bounds, const gp_Dir n) const = 0;

    virtual math::Matrix get_MnMn(const MeshElement* const e2, const std::vector<double>& u_ext, const std::vector<gp_Pnt>& bounds, const gp_Dir n) const = 0;
    virtual math::Matrix get_MnMn_log(const MeshElement* const e2, const std::vector<double>& u_ext, const std::vector<gp_Pnt>& bounds, const gp_Dir n, const double C, const double K) const = 0;
    virtual double get_log_integ(const MeshElement* const e2, const std::vector<double>& u_ext, const std::vector<gp_Pnt>& bounds, const gp_Dir n, const double C, const double K) const = 0;
    virtual double get_log_integ_deriv(const MeshElement* const e2, const std::vector<long>& node_positions, const std::vector<double>& u, const std::vector<double>& du, const std::vector<gp_Pnt>& bounds, const gp_Dir n, const double C, const double K) const = 0;
    virtual math::Matrix Kue_log(const MeshElement* const e2, const std::vector<double>& u_ext, const std::vector<gp_Pnt>& bounds, const gp_Dir n, const double C, const double K) const = 0;
    virtual void Ku_log(const double mult, const MeshElement* const e2, const std::vector<long>& node_positions, const std::vector<double>& u_ext, const std::vector<gp_Pnt>& bounds, const gp_Dir n, std::vector<double>& Ku, const double C, const double K) const = 0;
    virtual void dKu_log(const double mult, const MeshElement* const e2, const std::vector<long>& node_positions, const std::vector<double>& u, const std::vector<double>& du, const std::vector<gp_Pnt>& bounds, const gp_Dir n, std::vector<double>& dKu, const double C, const double K) const = 0;

    virtual math::Matrix get_Ni(const gp_Pnt& p) const = 0;
    virtual math::Vector get_Ni_1dof(const gp_Pnt& p) const{
        (void) p;
        return math::Vector();
    }

    /**
     * Creates and returns the elemental stiffness matrix.
     *
     * @param D Constitutive matrix.
     * @param t Geometry thickness.
     *
     * @return The element stiffness matrix.
     */
    virtual math::Matrix get_k(const math::Matrix& D, const double t) const = 0;
    /**
     * Creates and returns the elemental Robin matrix.
     *
     * @param K Stiffness matrix.
     * @param t Geometry thickness.
     * @param points Points that define the boundary.
     *
     * @return The element stiffness matrix.
     */
    virtual math::Matrix get_R(const math::Matrix& K, const double t, const std::vector<gp_Pnt>& points) const = 0;
    /**
     * Calculates the internal load vector of the element.
     *
     * @param D Constitutive matrix.
     * @param t Geometry thickness.
     * @param u Displacement vector.
     *
     * @return Internal loads at node.
     */
    virtual math::Vector get_internal_loads(const math::Matrix& D, const double t, const std::vector<double>& u) const = 0;
    /** 
     * Calculates the Von Mises stresses at a point within the element.
     *
     * @param D Constitutive matrix.
     * @param p The point.
     * @param u Displacement vector.
     *
     * @return Von Mises stress at point p.
     */
    virtual double get_stress_at(const math::Matrix& D, const gp_Pnt& p, const std::vector<double>& u, const double eps = 0) const = 0;
    /** 
     * Calculates the Von Mises strain at a point within the element.
     *
     * @param p The point.
     * @param u Displacement vector.
     *
     * @return Von Mises stress at point p.
     */
    virtual double get_strain_VM(const gp_Pnt& p, const std::vector<double>& u, const double eps = 0) const = 0;
    virtual math::Vector get_principal_strains(const gp_Pnt& p, const std::vector<double>& u) const = 0;
    /**
     * Calculates the derivative of von Mises stress with D(rho).
     *
     * @param D Constitutive matrix.
     * @param dD Derivative of the constitutive matrix.
     * @param mult Value to be multiplied to the result.
     * @param point Point to measure stress.
     * @param u Displacement vector.
     */
    virtual double von_Mises_derivative(const math::Matrix& D, const math::Matrix& dD, double mult, const gp_Pnt& point, const std::vector<double>& u) const = 0;
    /**
     * Calculates the Cauchy stress tensor at a point within the element.
     *
     * @param D Constitutive matrix.
     * @param p The point.
     * @param u Displacement vector.
     *
     * @return Stress tensor at point p.
     */
    virtual math::Matrix get_stress_tensor(const math::Matrix& D, const gp_Pnt& p, const std::vector<double>& u) const = 0;

    /**
     * Calculates the Cauchy strain tensor at a point within the element.
     *
     * @param p The point.
     * @param u Displacement vector.
     *
     * @return Strain tensor at point p.
     */
    virtual math::Matrix get_strain_tensor(const gp_Pnt& p, const std::vector<double>& u) const = 0;
    /**
     * Calculates the strain vector at a point within the element.
     *
     * @param p The point.
     * @param u Displacement vector.
     *
     * @return Strain vector at point p.
     */
    virtual math::Vector get_strain_vector(const gp_Pnt& p, const std::vector<double>& u) const = 0;
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
     * @param t Geometry thickness.
     *
     * @return volume or area.
     */
    virtual double get_volume(const double t) const = 0;
    /**
     * Calculates the element's compliance.
     *
     * @param D Constitutive matrix.
     * @param t Geometry thickness.
     * @param u Displacement vector.
     *
     * @return The compliance.
     */
    virtual double get_compliance(const math::Matrix& D, const double t, const std::vector<double>& u) const = 0;
    /**
     * Calculates the element's compliance with a virtual displacement vector.
     *
     * @param D Constitutive matrix.
     * @param t Geometry thickness.
     * @param u Displacement vector.
     * @param l Virtual displacement vector.
     *
     * @return The compliance.
     */
    virtual double get_compliance(const math::Matrix& D, const double t, const std::vector<double>& u, const std::vector<double>& l) const = 0;
    /**
     * Calculates the elemental contribution to the virtual load, for a p-norm
     * global stress aggregation function.
     *
     * @param D Constitutive matrix.
     * @param mult Value to be multiplied to the result. For p-norm stress 
     * based on LE et al, 2009, it is:
     * P*v*std::pow(data->new_x[i], pt*P)*std::pow(S, P-2)
     * @param u Displacement vector.
     * @param l Global virtual load vector.
     */
    virtual void get_virtual_load(const math::Matrix& D, double mult, const gp_Pnt& point, const std::vector<double>& u, std::vector<double>& l) const = 0;
    /**
     * Calculates the force vector for a distributed load.
     *
     * @param t Geometry thickness.
     * @param vec Load vector.
     * @param points Delimiting points.
     * @return Force vector.
     */
    virtual math::Vector get_f(const double t, const math::Vector& vec, const std::vector<gp_Pnt>& points) const = 0;
    /**
     * Returns the geometry of the element.
     *
     * @return Shape of the element.
     */
    virtual TopoDS_Shape get_shape() const = 0;
    /**
     * Returns the geometry of the element, considering nodal displacements.
     *
     * @param disp Displacement of nodes.
     *
     * @return Shape of the element.
     */
    virtual TopoDS_Shape get_shape(const std::vector<gp_Vec>& disp) const = 0;
    /**
     * Calculates the centroid of the element.
     *
     * @return The centroid.
     */
    virtual gp_Pnt get_centroid() const = 0;
    /**
     * Returns a factory pointer, from which it'll possible to obtain
     * information on the element being used.
     *
     * It returns a unique_ptr because it requires a pointer (due to poly-
     * morphism and template stuff), but it's generated on the fly, so that
     * the element doesn't have to store a pointer just to return it.
     *
     * Dirty and kind of slow, but it works. Just avoid using it. It's better
     * to plan ahead and use the `elem_info`variable from Meshing or the
     * 'topopt_element` variable from ProjectData.
     *
     * @return Unique pointer to MeshElementFactory instance.
     */
    virtual inline std::unique_ptr<MeshElementFactory> get_element_info() const = 0;

    /**
     * Returns the gradient of the shape nodal density shape functions at a 
     * point (commonly the centroid).
     *
     * @param p Point where the gradient is measured.
     *
     * @return Gradient matrix (2xN for 2D, 3xN for 3D)
     */
    virtual math::Matrix get_nodal_density_gradient(gp_Pnt p) const = 0;

    /**
     * Gets the linear displacement matrix (B).
     *
     * @param point Point where it's measured at.
     *
     * @return B matrix.
     */
    virtual math::Matrix get_B(const gp_Pnt& point) const = 0;

    virtual math::Matrix get_dk_sh(const math::Matrix& D, const double t, const size_t n, const size_t dof) const{
        (void) D;
        (void) t;
        (void) n;
        (void) dof;
        return math::Matrix();
    }
    virtual double get_dV_sh(const double t, const size_t n, const size_t dof) const{
        (void) t;
        (void) n;
        (void) dof;

        return 0;
    }
    virtual math::Matrix get_dB_sh(const gp_Pnt& p, const size_t n, const size_t dof) const{
        (void) p;
        (void) n;
        (void) dof;

        return math::Matrix();
    }
    virtual void calculate_coefficients(){

    }
    /**
     * Calculates the shape derivative of von Mises stress.
     *
     * @param D Constitutive matrix.
     * @param mult Value to be multiplied to the result.
     * @param point Point to measure stress.
     * @param u Displacement vector.
     * @param n Element's node id.
     * @param dof Analyzed degree of freedom.
     */
    virtual double von_Mises_derivative_sh(const math::Matrix& D, double mult, const gp_Pnt& point, const std::vector<double>& u, const size_t n, const size_t dof) const = 0;

    virtual math::Matrix diffusion_dim_dof(const double t, const math::Matrix& A) const = 0;
    virtual math::Matrix get_dNi(const gp_Pnt& p) const{
        (void)p;
        return math::Matrix();
    }


    /**
     * Returns a 1 degree of freedom diffusion matrix.
     * 
     * @param t Geometry thickness.
     * @param A Diffusivity matrix.
     *
     * @return Diffusion matrix.
     */
    virtual math::Matrix diffusion_1dof(const double t, const math::Matrix& A) const = 0;

    /**
     * Returns a 1 degree of freedom advection matrix.
     *
     * @param t Geometry thickness.
     * @param v Direction of advection.
     *
     * @return Advection matrix.
     */
    virtual math::Matrix advection_1dof(const double t, const math::Vector& v) const = 0;

    /**
     * Returns a 1 degree of freedom absorption matrix.
     *
     * @param t Geometry thickness.
     *
     * @return Absorption matrix.
     */
    virtual math::Matrix absorption_1dof(const double t) const = 0;

    virtual math::Matrix robin_1dof(const double t, const std::vector<gp_Pnt>& points) const{
        (void)t; (void)points; return math::Matrix();
    }

    /**
     * Returns a 1 degree of freedom source vector.
     *
     * @param t Geometry thickness.
     *
     * @return source matrix.
     */
    virtual math::Vector source_1dof(const double t) const = 0;

    /**
     * Returns a 1 degree of freedom flow vector.
     *
     * @param t Geometry thickness.
     * @param nodes Delimiting boundary nodes
     *
     * @return Flow vector.
     */
    virtual math::Vector flow_1dof(const double t, const MeshNode** nodes) const = 0;

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
    MeshElement(const std::vector<MeshNode*>& nodes):
        Element(std::vector<Node*>(nodes.begin(), nodes.end()))
        {}

    inline double H(const double x, const double C, const double K) const{
        if(std::abs(K*x) < 345){
            const double ln = std::log(1 + std::exp(K*x));
            const double h = C*ln*(ln/(K*K));
            return h;
        } else if(x > 0){
            return C*x;
        } else {
            return 0;
        }
    }
    inline double dH(const double x, const double C, const double K) const{
        if(std::abs(K*x) < 345){
            const double ekx1 = 1 + std::exp(K*x);
            const double ekx_1 = 1 + std::exp(-K*x);
            const double ln = std::log(ekx1);
            const double dh = C*(ln/(K*ekx_1));
            return dh;
        } else if(x > 0){
            return C;
        } else {
            return 0;
        }
    }
    inline double ddH(const double x, const double C, const double K) const{
        if(std::abs(K*x) < 345){
            const double ekx1 = 1 + std::exp(K*x);
            const double e_kx = std::exp(-K*x);
            const double ekx_1 = 1 + e_kx;
            const double ln = std::log(ekx1);
            const double ddh = C*((ln*e_kx + 1)/(ekx_1*ekx_1));
            return ddh;
        } else {
            return 0;
        }
    }

    /**
     * Gets the matrix used to calculate the nodal force vector. Calculated
     * from the boundary integral of the interpolation matrix.
     *
     * @param t Geometry thickness.
     * @param points Points that define the boundary.
     *
     * @return The interpolation matrix for nodal forces, Nf.
     */
    virtual math::Matrix get_Nf(const double t, const std::vector<gp_Pnt>& points) const = 0;
};

namespace utils{
    template<size_t Y, size_t Z>
    class BoundaryNullifier;
}

class BoundaryMeshElement : public Element{
    public:
    const MeshElement* const parent;

    virtual ~BoundaryMeshElement() = default;

    virtual math::Matrix get_K_ext(const math::Matrix& D, const gp_Pnt& center) const = 0;
    virtual math::Matrix get_stress_integrals(const math::Matrix& D, const gp_Pnt& center) const = 0;
    virtual math::Matrix get_equilibrium_partial(const math::Matrix& D, const gp_Pnt& center, const std::vector<size_t>& stresses) const = 0;
    virtual math::Vector get_dz_vector(const math::Matrix& S, const math::Matrix& D, const double Az, const double Bz, const gp_Pnt& center) const = 0;
    virtual math::Matrix get_dz_vector_matrix(const math::Matrix& S, const math::Matrix& D, const gp_Pnt& center) const = 0;
    virtual math::Matrix get_dz_vector_matrix_1d(const math::Matrix& S, const math::Matrix& D, const math::Matrix& dS, const math::Matrix& dD, const gp_Pnt& center) const = 0;
    virtual math::Vector get_force_vector(const math::Matrix& D, const std::vector<double>& u, const gp_Pnt& center, const math::Matrix& rot, const bool transpose = false) const = 0;
    virtual math::Matrix get_force_vector_matrix(const math::Matrix& D, const gp_Pnt& center, const math::Matrix& rot) const = 0;

    /**
     * Returns a 1 degree of freedom diffusion matrix.
     * 
     * @param A Diffusivity matrix.
     *
     * @return Diffusion matrix.
     */
    virtual math::Matrix diffusion_1dof(const math::Matrix& A) const = 0;

    /**
     * Returns a 1 degree of freedom advection matrix.
     *
     * @param v Direction of advection.
     *
     * @return Advection matrix.
     */
    virtual math::Matrix advection_1dof(const math::Vector& v) const = 0;

    /**
     * Returns a 1 degree of freedom absorption matrix.
     *
     * @return Absorption matrix.
     */
    virtual math::Matrix absorption_1dof() const = 0;

    /**
     * Returns a 1 degree of freedom source vector.
     *
     * @return source matrix.
     */
    virtual math::Vector source_1dof() const = 0;
    virtual math::Vector flow_1dof(const std::array<const Node*, 2>& nodes) const = 0;

    virtual double get_area() const = 0;

    /**
     * Calculates the centroid of the element.
     *
     * @return The centroid.
     */
    virtual gp_Pnt get_centroid() const = 0;
    /**
     * Calculates the normal of the element.
     *
     * @return The normal.
     */
    virtual gp_Dir get_normal() const = 0;

    protected:
    /**
     * Creates an element with the specified nodes.
     *
     * @param nodes List of nodes.
     */
    BoundaryMeshElement(const std::vector<MeshNode*>& nodes, const MeshElement* const parent):
        Element(std::vector<Node*>(nodes.begin(), nodes.end())), parent(parent)
        {}
};

class ContactMeshElement : public Element{
    public:
    const MeshElement* const e1;
    const MeshElement* const e2;

    virtual ~ContactMeshElement() = default;

    virtual math::Matrix get_frictionless_Ge() const = 0;

    virtual math::Matrix fl2_uu(const math::Vector& l_e, const math::Vector& u1, const math::Vector& u2) const = 0;
    virtual math::Matrix fl2_uL(const math::Vector& u1, const math::Vector& u2) const = 0;

    virtual void fl2_Ku_lambda(const double EPS, const std::vector<long> u1_pos, const std::vector<long> u2_pos, const std::vector<long>& lu_pos, const std::vector<double>& u, std::vector<double>& Ku) const = 0;
    virtual void fl2_dKu_lambda(const double EPS, const std::vector<long> u1_pos, const std::vector<long> u2_pos, const std::vector<long>& lu_pos, const std::vector<double>& u, const std::vector<double>& du, std::vector<double>& Ku) const = 0;


    virtual double get_area() const = 0;

    /**
     * Calculates the centroid of the element.
     *
     * @return The centroid.
     */
    virtual gp_Pnt get_centroid() const = 0;
    /**
     * Calculates the normal of the element.
     *
     * @return The normal.
     */
    virtual gp_Dir get_normal() const = 0;
    virtual const math::Matrix& get_R() const = 0;

    protected:
    /**
     * Creates an element with the specified nodes.
     *
     * @param nodes List of nodes.
     */
    ContactMeshElement(const std::vector<MeshNode*>& nodes, const MeshElement* const e1, const MeshElement* const e2):
        Element(std::vector<Node*>(nodes.begin(), nodes.end())), e1(e1), e2(e2)
        {}
};

#endif
