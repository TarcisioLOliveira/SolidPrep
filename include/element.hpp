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

#include <gp_Pnt.hxx>
#include <vector>

class Node{
    public:
    const gp_Pnt point;
    const long id;
    double * const results;
    ~Node(){ delete[] results; }

    protected:
    Node(gp_Pnt p, long id, size_t dim):point(p), id(id), results(new double[dim]()){}
};

class BeamNode : public Node{
    public:
    const double dim;
    const gp_Dir normal;
    virtual ~BeamNode() = default;
    protected:
    BeamNode(gp_Pnt p, long id, size_t res_n, double dim, gp_Dir n):Node(p, id, res_n),dim(dim), normal(n){}
};

class BeamNode2D : public BeamNode{
    public:
    BeamNode2D(gp_Pnt p, long id, double dim, gp_Dir n):BeamNode(p, id, 3, dim, n){}
};

class BeamNodeFactory{
    public:
    enum BeamNodeType{
        NONE,
        BEAM_NODE_2D
    };
    static BeamNode* make_node(gp_Pnt p, long id, double dim, gp_Dir n, BeamNodeType t);
};

class MeshNode : public Node{
    public:
    virtual ~MeshNode() = default;
    protected:
    MeshNode(gp_Pnt p, long id, size_t res_n):Node(p, id, res_n){}
};

class MeshNode2D : public MeshNode{
    public:
    MeshNode2D(gp_Pnt p, long id):MeshNode(p, id, 3){}
};

class MeshNodeFactory{
    public:
    enum MeshNodeType{
        NONE,
        MESH_NODE_2D
    };
    static MeshNode* make_node(gp_Pnt p, long id, MeshNodeType t);
};

class Element{
    public:

    const std::vector<Node*> nodes;
    const std::vector<long> u_pos;

    virtual ~Element() = default;
    virtual std::vector<double> get_k() const = 0;

    protected:
    Element(std::vector<Node*> n, std::vector<long> u_pos):nodes(std::move(n)), u_pos(std::move(u_pos)){}
};

class BeamElement : public Element{
    public:

    virtual ~BeamElement() = default;
    virtual std::vector<double> get_k() const override = 0;
    virtual BeamNode* get_internal_loads(size_t node, const std::vector<double>& u) const = 0;
    virtual inline BeamNode* get_node(size_t node) const{ return static_cast<BeamNode*>(this->nodes[node]);}

    protected:
    BeamElement(BeamNode* p1, BeamNode* p2, std::vector<long> u_pos):Element({p1,p2}, std::move(u_pos)){}
};

class MeshElement : public Element{
    public:

    virtual ~MeshElement() = default;
    virtual std::vector<double> get_k() const override = 0;
    virtual MeshNode* get_stresses(size_t node, const std::vector<double>& u) const = 0;

    protected:
    MeshElement(std::vector<MeshNode*> nodes, std::vector<long> u_pos):Element(std::vector<Node*>(nodes.begin(), nodes.end()), std::move(u_pos)){}
};

#endif
