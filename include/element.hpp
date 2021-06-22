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

#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <cmath>
#include <gp_Pnt.hxx>
#include <vector>
#include <TopoDS_Shape.hxx>

class Node{
    public:
    const gp_Pnt point;
    size_t id;
    float * const results;
    long * u_pos;
    ~Node(){ delete[] results; }

    protected:
    Node(gp_Pnt p, size_t id, size_t dim):point(p), id(id), results(new float[dim]()), u_pos(nullptr){}
};

class BeamNode : public Node{
    public:
    const float dim;
    const gp_Dir normal;
    virtual ~BeamNode() = default;
    protected:
    BeamNode(gp_Pnt p, size_t id, size_t res_n, float dim, gp_Dir n):Node(p, id, res_n),dim(dim), normal(n){}
};

class BeamNode2D : public BeamNode{
    public:
    BeamNode2D(gp_Pnt p, size_t id, float dim, gp_Dir n):BeamNode(p, id, 3, dim, n){}
};

class BeamNodeFactory{
    public:
    enum BeamNodeType{
        NONE,
        BEAM_NODE_2D
    };
    static BeamNode* make_node(gp_Pnt p, size_t id, float dim, gp_Dir n, BeamNodeType t);
};

class MeshNode : public Node{
    public:
    virtual ~MeshNode() = default;
    virtual size_t get_result_size() const = 0;
    virtual float get_Von_Mises() const = 0;
    protected:
    MeshNode(gp_Pnt p, size_t id, size_t res_n): Node(p, id, res_n){}
};

class MeshNode2D : public MeshNode{
    public:
    MeshNode2D(gp_Pnt p, size_t id):MeshNode(p, id, 3){}
    virtual size_t get_result_size() const override{ return 3; }
    virtual float get_Von_Mises() const override{
        return std::sqrt(std::pow(results[0], 2) - results[0]*results[1] + std::pow(results[1], 2) + 3*std::pow(results[2], 2));
    }
};

class MeshNodeFactory{
    public:
    enum MeshNodeType{
        NONE,
        MESH_NODE_2D
    };
    static MeshNode* make_node(gp_Pnt p, size_t id, MeshNodeType t);
};


struct ElementShape{
    std::vector<MeshNode*> nodes;
};

class Element{
    public:

    const std::vector<Node*> nodes;

    virtual ~Element() = default;
    virtual std::vector<float> get_k() const = 0;

    protected:
    Element(std::vector<Node*> n):nodes(std::move(n)){}
};

class BeamElement : public Element{
    public:

    virtual ~BeamElement() = default;
    virtual std::vector<float> get_k() const override = 0;
    virtual BeamNode* get_internal_loads(size_t node, const std::vector<float>& u) const = 0;
    virtual inline BeamNode* get_node(size_t node) const{ return static_cast<BeamNode*>(this->nodes[node]);}

    protected:
    BeamElement(BeamNode* p1, BeamNode* p2):Element({p1,p2}){}
};

class MeshElement : public Element{
    public:

    virtual ~MeshElement() = default;
    virtual std::vector<float> get_k() const override = 0;
    virtual MeshNode* get_stresses(size_t node, const std::vector<float>& u) const = 0;
    virtual MeshNode* get_internal_loads(size_t node, const std::vector<float>& u) const = 0;
    virtual double get_volume() const = 0;
    virtual TopoDS_Shape get_shape() const = 0;
    virtual gp_Pnt get_centroid() const = 0;

    /**
     * Necessary to use Gmsh's GUI. See:
     * https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
     * 
     * @return Gmsh element type number.
     */
    virtual size_t get_gmsh_element_type() const = 0;

    virtual inline MeshNode* get_node(size_t node) const{ return static_cast<MeshNode*>(this->nodes[node]);}

    protected:
    MeshElement(std::vector<MeshNode*> nodes):Element(std::vector<Node*>(nodes.begin(), nodes.end())){}
};

#endif
