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

#ifndef ELEMENT_COMMON_HPP
#define ELEMENT_COMMON_HPP

#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS.hxx>
#include <IntTools_EdgeEdge.hxx>
#include <Standard_Handle.hxx>
#include <vector>
#include <BRepOffsetAPI_Sewing.hxx>

#include "element.hpp"
#include "logger.hpp"
#include "utils.hpp"

template<class T>
class MeshElementCommon : public MeshElement{
    public:
    virtual ~MeshElementCommon() = default;

    virtual std::vector<double> get_internal_loads(const std::vector<double>& D, const double t, const std::vector<double>& u) const override{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        const std::vector<double> k = this->get_k(D, t);

        std::vector<double> results(K_DIM, 0);
        std::vector<double> U(K_DIM, 0);
        for(size_t l = 0; l < N; ++l){
            for(size_t j = 0; j < NODE_DOF; ++j){
                if(this->nodes[l]->u_pos[j] > -1){
                    U[l*NODE_DOF + j] = u[this->nodes[l]->u_pos[j]];
                }
            }
        }
        for(size_t i = 0; i < K_DIM; ++i){
            for(size_t j = 0; j < K_DIM; ++j){
                results[i] += k[i*K_DIM + j]*U[j];
            }
        }

        return results;
    }

    virtual double get_compliance(const std::vector<double>& D, const double t, const std::vector<double>& u) const override{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        const auto k = this->get_k(D, t);

        std::vector<double> u_vec(K_DIM, 0);
        for(size_t i = 0; i < N; ++i){
            for(size_t j = 0; j < NODE_DOF; ++j){
                if(this->nodes[i]->u_pos[j] > -1){
                    u_vec[i*NODE_DOF+j] = u[this->nodes[i]->u_pos[j]];
                }
            }
        }

        std::vector<double> f_vec(K_DIM, 0);

        for(size_t i = 0; i < K_DIM; ++i){
            for(size_t j = 0; j < K_DIM; ++j){
                f_vec[i] += k[i*K_DIM + j]*u_vec[j];
            }
        }

        double result = 0;
        for(size_t i = 0; i < K_DIM; ++i){
            result += f_vec[i]*u_vec[i];
        }
        // return cblas_ddot(K_DIM, u_vec.data(), 1, f_vec.data(), 1);
        
        return result;
    }
    virtual double get_compliance(const std::vector<double>& D, const double t, const std::vector<double>& u, const std::vector<double>& l) const override{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        const auto k = this->get_k(D, t);

        std::vector<double> u_vec(K_DIM, 0);
        for(size_t i = 0; i < N; ++i){
            for(size_t j = 0; j < NODE_DOF; ++j){
                if(this->nodes[i]->u_pos[j] > -1){
                    u_vec[i*NODE_DOF+j] = u[this->nodes[i]->u_pos[j]];
                }
            }
        }

        std::vector<double> f_vec(K_DIM, 0);

        for(size_t i = 0; i < K_DIM; ++i){
            for(size_t j = 0; j < K_DIM; ++j){
                f_vec[i] += k[i*K_DIM + j]*u_vec[j];
            }
        }

        double result = 0;
        std::vector<double> l_vec(K_DIM, 0);
        for(size_t i = 0; i < N; ++i){
            for(size_t j = 0; j < NODE_DOF; ++j){
                if(this->nodes[i]->u_pos[j] > -1){
                    l_vec[i*NODE_DOF+j] = l[this->nodes[i]->u_pos[j]];
                }
            }
        }
        for(size_t i = 0; i < K_DIM; ++i){
            result += f_vec[i]*l_vec[i];
        }
        // return cblas_ddot(K_DIM, l_vec.data(), 1, f_vec.data(), 1);
        return result;
    }

    virtual gp_Pnt get_centroid() const override{
        const size_t N = T::NODES_PER_ELEM;

        double x = 0;
        double y = 0;
        double z = 0;
        for(size_t i = 0; i < N; ++i){
            const auto& n = this->nodes[i];
            x += n->point.X();
            y += n->point.Y();
            z += n->point.Z();
        }

        return gp_Pnt(x/N, y/N, z/N);
    }

    protected:
    MeshElementCommon(const std::vector<MeshNode*>& nodes):
            MeshElement(nodes)
            {}

    virtual double _von_Mises_derivative(const std::vector<double>& D, const std::vector<double>& dD, double mult, const gp_Pnt& point, const std::vector<double>& u, const std::vector<double>& V, const size_t S_SIZE) const{
        auto s = this->_get_stress_vector(D, point, u, S_SIZE);
        auto ds = this->_get_stress_vector(dD, point, u, S_SIZE);
        std::vector<double> vec1(s.size(), 0);
        for(size_t i = 0; i < S_SIZE; ++i){
            for(size_t j = 0; j < S_SIZE; ++j){
                vec1[i] += V[i*S_SIZE + j]*ds[j];
            }
        }

        double result = 0;
        for(size_t i = 0; i < S_SIZE; ++i){
            result += s[i]*ds[i];
        }
        result *= mult;

        return result;
    }

    virtual void _get_virtual_load(const std::vector<double>& D, double mult, const gp_Pnt& point, const std::vector<double>& u, std::vector<double>& l, const std::vector<double>& V, const size_t S_SIZE) const{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        const std::vector<double> DB = this->get_DB(D, point);

        std::vector<double> u_vec(K_DIM, 0);
        for(size_t i = 0; i < N; ++i){
            for(size_t j = 0; j < NODE_DOF; ++j){
                if(this->nodes[i]->u_pos[j] > -1){
                    u_vec[i*NODE_DOF+j] = u[this->nodes[i]->u_pos[j]];
                }
            }
        }

        std::vector<double> vec(S_SIZE, 0);

        for(size_t i = 0; i < S_SIZE; ++i){
            for(size_t j = 0; j < K_DIM; ++j){
                vec[i] += DB[i*K_DIM + j]*u_vec[j];
            }
        }
        std::vector<double> vec2(S_SIZE, 0);
        for(size_t i = 0; i < S_SIZE; ++i){
            for(size_t j = 0; j < S_SIZE; ++j){
                vec2[i] += V[i*S_SIZE + j]*vec[j];
            }
        }
        std::vector<double> vec3(K_DIM, 0);
        for(size_t i = 0; i < K_DIM; ++i){
            for(size_t j = 0; j < S_SIZE; ++j){
                vec3[i] += DB[j*K_DIM + i]*vec2[j];
            }
        }

        for(size_t k = 0; k < N; ++k){
            for(size_t j = 0; j < NODE_DOF; ++j){
                if(this->nodes[k]->u_pos[j] > -1){
                    l[this->nodes[k]->u_pos[j]] += mult*vec3[k*NODE_DOF+j];
                }
            }
        }
    }

    virtual std::vector<double> _get_stress_vector(const std::vector<double>& D, const gp_Pnt& p, const std::vector<double>& u, const size_t S_SIZE) const{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        const std::vector<double> DB = this->get_DB(D, p);

        std::vector<double> results(S_SIZE, 0);
        std::vector<double> U(K_DIM, 0);
        for(size_t l = 0; l < N; ++l){
            for(size_t j = 0; j < NODE_DOF; ++j){
                if(this->nodes[l]->u_pos[j] > -1){
                    U[l*NODE_DOF + j] = u[this->nodes[l]->u_pos[j]];
                }
            }
        }
        for(size_t i = 0; i < S_SIZE; ++i){
            for(size_t j = 0; j < K_DIM; ++j){
                results[i] += DB[i*K_DIM + j]*U[j];
            }
        }

        return results;
    }


    virtual std::vector<double> _get_stress_tensor(const std::vector<double>& D, const gp_Pnt& p, const std::vector<double>& u, const std::vector<size_t> indices, const size_t S_SIZE) const{

        auto results = this->_get_stress_vector(D, p, u, S_SIZE);

        std::vector<double> S(indices.size());
        for(size_t i = 0; i < indices.size(); ++i){
            S[i] = results[indices[i]];
        }

        return S;
    }
};


template<class T>
class MeshElementCommon2D : public MeshElementCommon<T>{
    public:
    static const size_t S_SIZE = 3; // Size of the stress and strain vectors
    static const size_t DIM    = 2; // Number of dimensions

    static const utils::ProblemType PROBLEM_TYPE = utils::PROBLEM_TYPE_2D;

    virtual ~MeshElementCommon2D() = default;

    virtual double get_stress_at(const std::vector<double>& D, const gp_Pnt& point, const std::vector<double>& u) const override{
        auto results = this->_get_stress_vector(D, point, u, S_SIZE);

        return std::sqrt(results[0]*results[0] - results[0]*results[1] + results[1]*results[1] + 3*results[2]*results[2]);
    }

    virtual std::vector<double> get_stress_tensor(const std::vector<double>& D, const gp_Pnt& p, const std::vector<double>& u) const override{
        const std::vector<size_t> indices{0, 2, 2, 1};
        return this->_get_stress_tensor(D, p, u, indices, S_SIZE);
    }

    virtual void get_virtual_load(const std::vector<double>& D, double mult, const gp_Pnt& point, const std::vector<double>& u, std::vector<double>& l) const override{
        const std::vector<double> V{1, -0.5, 0,
                                   -0.5, 1, 0,
                                   0,   0, 3};

        this->_get_virtual_load(D, mult, point, u, l, V, S_SIZE);
    }

    virtual double von_Mises_derivative(const std::vector<double>& D, const std::vector<double>& dD, double mult, const gp_Pnt& point, const std::vector<double>& u) const override{
        const std::vector<double> V{1, -0.5, 0,
                                   -0.5, 1, 0,
                                   0,   0, 3};

        return this->_von_Mises_derivative(D, dD, mult, point, u, V, S_SIZE);
    }

    virtual std::vector<gp_Pnt> get_intersection_points(const TopoDS_Shape& crosssection) const override{
        const size_t N = T::NODES_PER_ELEM;
        std::vector<gp_Pnt> points;

        const TopoDS_Edge line_edge = TopoDS::Edge(crosssection);

        for(size_t i = 0; i < N; ++i){
            size_t j = (i + 1) % N;

            const TopoDS_Edge e = BRepBuilderAPI_MakeEdge(this->nodes[i]->point, this->nodes[j]->point);
            const gp_Dir dir(gp_Vec(this->nodes[i]->point, this->nodes[j]->point));

            IntTools_EdgeEdge tool(line_edge, e);
            tool.Perform();
            auto common = tool.CommonParts();
            if(common.Size() > 0){
                for(auto& c:common){
                    double dist = std::abs(c.VertexParameter2()); // n.normal is a unit vector, and the line starts at 0
                    gp_Pnt p = this->nodes[i]->point.Translated(dist*dir);
                    bool contained = false;
                    for(auto& pp:points){
                        if(p.IsEqual(pp, Precision::Confusion())){
                            contained = true;
                            break;
                        }
                    }
                    if(!contained){
                        points.push_back(p);
                    }
                }
            }
        }

        return points;
    }

    virtual std::vector<double> get_f(const double t, const gp_Dir& dir, double norm, const std::vector<gp_Pnt>& points) const override{
        const size_t K_DIM = T::K_DIM;

        const double px = dir.X()*norm;
        const double py = dir.Y()*norm;

        const auto Nf = this->get_Nf(t, points);

        std::vector<double> f(K_DIM, 0);
        for(size_t i = 0; i < K_DIM; ++i){
            f[i] = Nf[T::DIM*i]*px + Nf[T::DIM*i+1]*py;
        }

        return f;
    }

    protected:
    MeshElementCommon2D(const std::vector<MeshNode*>& nodes):
            MeshElementCommon<T>(nodes)
            {}
};

template<class T>
class MeshElementCommon2DTri : public MeshElementCommon2D<T>{
    public:
    virtual ~MeshElementCommon2DTri() = default;

    static const Element::Shape SHAPE_TYPE = Element::Shape::TRI;

    virtual double get_volume(const double t) const override{
        const gp_Mat deltaM(1, this->nodes[0]->point.X(), this->nodes[0]->point.Y(),
                            1, this->nodes[1]->point.X(), this->nodes[1]->point.Y(),
                            1, this->nodes[2]->point.X(), this->nodes[2]->point.Y());

        return 0.5*std::abs(deltaM.Determinant())*t;
    }

    virtual TopoDS_Shape get_shape() const override{
        const gp_Pnt p1 = this->nodes[0]->point;
        const gp_Pnt p2 = this->nodes[1]->point;
        const gp_Pnt p3 = this->nodes[2]->point;

        return this->generate_geometry(p1, p2, p3);
    }

    virtual TopoDS_Shape get_shape(const std::vector<gp_Vec>& disp) const override{
        const gp_Pnt p1 = this->nodes[0]->point.Translated(disp[0]);
        const gp_Pnt p2 = this->nodes[1]->point.Translated(disp[1]);
        const gp_Pnt p3 = this->nodes[2]->point.Translated(disp[2]);

        return this->generate_geometry(p1, p2, p3);
    }

    protected:
    MeshElementCommon2DTri(const std::vector<MeshNode*>& nodes):
            MeshElementCommon2D<T>(nodes)
            {}

    inline TopoDS_Face generate_geometry(const gp_Pnt& p1, const gp_Pnt& p2, const gp_Pnt& p3) const{
        const TopoDS_Vertex v1 = BRepBuilderAPI_MakeVertex(p1);
        const TopoDS_Vertex v2 = BRepBuilderAPI_MakeVertex(p2);
        const TopoDS_Vertex v3 = BRepBuilderAPI_MakeVertex(p3);

        const TopoDS_Edge e1 = BRepBuilderAPI_MakeEdge(v1, v2);
        const TopoDS_Edge e2 = BRepBuilderAPI_MakeEdge(v2, v3);
        const TopoDS_Edge e3 = BRepBuilderAPI_MakeEdge(v3, v1);

        const TopoDS_Wire w = BRepBuilderAPI_MakeWire(e1, e2, e3);

        const TopoDS_Face f = BRepBuilderAPI_MakeFace(w);

        return f;
    }
};

template<class T>
class MeshElementCommon2DQuad : public MeshElementCommon2D<T>{
    public:
    virtual ~MeshElementCommon2DQuad() = default;

    static const Element::Shape SHAPE_TYPE = Element::Shape::QUAD;

    virtual double get_volume(const double t) const override{
        const gp_Mat deltaM(1, this->nodes[0]->point.X(), this->nodes[0]->point.Y(),
                            1, this->nodes[1]->point.X(), this->nodes[1]->point.Y(),
                            1, this->nodes[2]->point.X(), this->nodes[2]->point.Y());

        const gp_Mat deltaM2(1, this->nodes[2]->point.X(), this->nodes[2]->point.Y(),
                             1, this->nodes[3]->point.X(), this->nodes[3]->point.Y(),
                             1, this->nodes[0]->point.X(), this->nodes[0]->point.Y());

        return 0.5*(std::abs(deltaM.Determinant()) + std::abs(deltaM2.Determinant()))*t;
    }

    virtual TopoDS_Shape get_shape() const override{
        const gp_Pnt p1 = this->nodes[0]->point;
        const gp_Pnt p2 = this->nodes[1]->point;
        const gp_Pnt p3 = this->nodes[2]->point;
        const gp_Pnt p4 = this->nodes[3]->point;

        return this->generate_geometry(p1, p2, p3, p4);
    }

    virtual TopoDS_Shape get_shape(const std::vector<gp_Vec>& disp) const override{
        const gp_Pnt p1 = this->nodes[0]->point.Translated(disp[0]);
        const gp_Pnt p2 = this->nodes[1]->point.Translated(disp[1]);
        const gp_Pnt p3 = this->nodes[2]->point.Translated(disp[2]);
        const gp_Pnt p4 = this->nodes[3]->point.Translated(disp[3]);

        return this->generate_geometry(p1, p2, p3, p4);
    }

    protected:
    MeshElementCommon2DQuad(const std::vector<MeshNode*>& nodes):
            MeshElementCommon2D<T>(nodes)
            {}

    inline TopoDS_Face generate_geometry(const gp_Pnt& p1, const gp_Pnt& p2, const gp_Pnt& p3, const gp_Pnt& p4) const{
        const TopoDS_Vertex v1 = BRepBuilderAPI_MakeVertex(p1);
        const TopoDS_Vertex v2 = BRepBuilderAPI_MakeVertex(p2);
        const TopoDS_Vertex v3 = BRepBuilderAPI_MakeVertex(p3);
        const TopoDS_Vertex v4 = BRepBuilderAPI_MakeVertex(p4);

        const TopoDS_Edge e1 = BRepBuilderAPI_MakeEdge(v1, v2);
        const TopoDS_Edge e2 = BRepBuilderAPI_MakeEdge(v2, v3);
        const TopoDS_Edge e3 = BRepBuilderAPI_MakeEdge(v3, v4);
        const TopoDS_Edge e4 = BRepBuilderAPI_MakeEdge(v4, v1);

        const TopoDS_Wire w = BRepBuilderAPI_MakeWire(e1, e2, e3, e4);

        const TopoDS_Face f = BRepBuilderAPI_MakeFace(w);

        return f;
    }
};


template<class T>
class MeshElementCommon3D : public MeshElementCommon<T>{
    public:
    static const size_t S_SIZE = 6; // Size of the stress and strain vectors
    static const size_t DIM    = 3; // Number of dimensions

    static const utils::ProblemType PROBLEM_TYPE = utils::PROBLEM_TYPE_3D;

    virtual ~MeshElementCommon3D() = default;

    virtual double get_stress_at(const std::vector<double>& D, const gp_Pnt& point, const std::vector<double>& u) const override{
        auto s = this->_get_stress_vector(D, point, u, S_SIZE);

        const double ds1 = s[0] - s[1];
        const double ds2 = s[1] - s[2];
        const double ds3 = s[2] - s[0];

        return std::sqrt(0.5*(ds1*ds1 + ds2*ds2 + ds3*ds3 + 6*(s[3]*s[3] + s[4]*s[4] + s[5]*s[5])));
    }

    virtual std::vector<double> get_stress_tensor(const std::vector<double>& D, const gp_Pnt& p, const std::vector<double>& u) const override{
        const std::vector<size_t> indices{0, 3, 4, 3, 1, 5, 4, 5, 3};
        return this->_get_stress_tensor(D, p, u, indices, S_SIZE);
    }

    virtual void get_virtual_load(const std::vector<double>& D, double mult, const gp_Pnt& point, const std::vector<double>& u, std::vector<double>& l) const override{
        const std::vector<double> V{   1, -0.5, -0.5, 0, 0, 0,
                                    -0.5,    1, -0.5, 0, 0, 0,
                                    -0.5, -0.5,    1, 0, 0, 0,
                                       0,    0,    0, 3, 0, 0,
                                       0,    0,    0, 0, 3, 0,
                                       0,    0,    0, 0, 0, 3};

        this->_get_virtual_load(D, mult, point, u, l, V, S_SIZE);
    }

    virtual double von_Mises_derivative(const std::vector<double>& D, const std::vector<double>& dD, double mult, const gp_Pnt& point, const std::vector<double>& u) const override{
        const std::vector<double> V{   1, -0.5, -0.5, 0, 0, 0,
                                    -0.5,    1, -0.5, 0, 0, 0,
                                    -0.5, -0.5,    1, 0, 0, 0,
                                       0,    0,    0, 3, 0, 0,
                                       0,    0,    0, 0, 3, 0,
                                       0,    0,    0, 0, 0, 3};

        return this->_von_Mises_derivative(D, dD, mult, point, u, V, S_SIZE);
    }

    virtual std::vector<double> get_f(const double t, const gp_Dir& dir, double norm, const std::vector<gp_Pnt>& points) const override{
        const size_t K_DIM = T::K_DIM;

        const double px = dir.X()*norm;
        const double py = dir.Y()*norm;
        const double pz = dir.Z()*norm;

        const auto Nf = this->get_Nf(t, points);

        std::vector<double> f(K_DIM, 0);
        for(size_t i = 0; i < K_DIM; ++i){
            f[i] = Nf[T::DIM*i]*px + Nf[T::DIM*i+1]*py + Nf[T::DIM*i+2]*pz;
        }

        return f;
    }


    protected:
    MeshElementCommon3D(const std::vector<MeshNode*>& nodes):
            MeshElementCommon<T>(nodes)
            {}
};


template<class T>
class MeshElementCommon3DTet : public MeshElementCommon3D<T>{
    public:
    virtual ~MeshElementCommon3DTet() = default;

    static const Element::Shape SHAPE_TYPE = Element::Shape::TRI;

    virtual double get_volume(const double t) const override{
        (void)t;
        const gp_Pnt p0 = this->nodes[0]->point;
        const gp_Pnt p1 = this->nodes[1]->point;
        const gp_Pnt p2 = this->nodes[2]->point;
        const gp_Pnt p3 = this->nodes[3]->point;

        const gp_Vec v0(p0, p1);
        const gp_Vec v1(p0, p2);
        const gp_Vec v2(p0, p3);

        return std::abs(v0.DotCross(v1, v2))/6;
    }

    virtual TopoDS_Shape get_shape() const override{
        const gp_Pnt p1 = this->nodes[0]->point;
        const gp_Pnt p2 = this->nodes[1]->point;
        const gp_Pnt p3 = this->nodes[2]->point;
        const gp_Pnt p4 = this->nodes[3]->point;

        return this->generate_geometry(p1, p2, p3, p4);
    }

    virtual TopoDS_Shape get_shape(const std::vector<gp_Vec>& disp) const override{
        const gp_Pnt p1 = this->nodes[0]->point.Translated(disp[0]);
        const gp_Pnt p2 = this->nodes[1]->point.Translated(disp[1]);
        const gp_Pnt p3 = this->nodes[2]->point.Translated(disp[2]);
        const gp_Pnt p4 = this->nodes[3]->point.Translated(disp[3]);

        return this->generate_geometry(p1, p2, p3, p4);
    }

    virtual std::vector<gp_Pnt> get_intersection_points(const TopoDS_Shape& crosssection) const override{
        logger::log_assert(false, logger::ERROR, "get_intersection_points is not implemented for 3D elements.");

        const size_t N = T::NODES_PER_ELEM;
        std::vector<gp_Pnt> points;

        const TopoDS_Edge line_edge = TopoDS::Edge(crosssection);

        for(size_t i = 0; i < N; ++i){
            size_t j = (i + 1) % N;

            const TopoDS_Edge e = BRepBuilderAPI_MakeEdge(this->nodes[i]->point, this->nodes[j]->point);
            const gp_Dir dir(gp_Vec(this->nodes[i]->point, this->nodes[j]->point));

            IntTools_EdgeEdge tool(line_edge, e);
            tool.Perform();
            auto common = tool.CommonParts();
            if(common.Size() > 0){
                for(auto& c:common){
                    double dist = std::abs(c.VertexParameter2()); // n.normal is a unit vector, and the line starts at 0
                    gp_Pnt p = this->nodes[i]->point.Translated(dist*dir);
                    bool contained = false;
                    for(auto& pp:points){
                        if(p.IsEqual(pp, Precision::Confusion())){
                            contained = true;
                            break;
                        }
                    }
                    if(!contained){
                        points.push_back(p);
                    }
                }
            }
        }

        return points;
    }

    protected:
    MeshElementCommon3DTet(const std::vector<MeshNode*>& nodes):
            MeshElementCommon3D<T>(nodes)
            {}

    inline TopoDS_Solid generate_geometry(const gp_Pnt& p1, const gp_Pnt& p2, const gp_Pnt& p3, const gp_Pnt& p4) const{
        const TopoDS_Vertex v1 = BRepBuilderAPI_MakeVertex(p1);
        const TopoDS_Vertex v2 = BRepBuilderAPI_MakeVertex(p2);
        const TopoDS_Vertex v3 = BRepBuilderAPI_MakeVertex(p3);
        const TopoDS_Vertex v4 = BRepBuilderAPI_MakeVertex(p4);

        const TopoDS_Edge e1 = BRepBuilderAPI_MakeEdge(v1, v2);
        const TopoDS_Edge e2 = BRepBuilderAPI_MakeEdge(v1, v3);
        const TopoDS_Edge e3 = BRepBuilderAPI_MakeEdge(v1, v4);
        const TopoDS_Edge e4 = BRepBuilderAPI_MakeEdge(v2, v3);
        const TopoDS_Edge e5 = BRepBuilderAPI_MakeEdge(v3, v4);
        const TopoDS_Edge e6 = BRepBuilderAPI_MakeEdge(v4, v2);

        const TopoDS_Wire w1 = BRepBuilderAPI_MakeWire(e1, e4, e2);
        const TopoDS_Face f1 = BRepBuilderAPI_MakeFace(w1);

        const TopoDS_Wire w2 = BRepBuilderAPI_MakeWire(e2, e5, e3);
        const TopoDS_Face f2 = BRepBuilderAPI_MakeFace(w2);

        const TopoDS_Wire w3 = BRepBuilderAPI_MakeWire(e3, e6, e1);
        const TopoDS_Face f3 = BRepBuilderAPI_MakeFace(w3);

        const TopoDS_Wire w4 = BRepBuilderAPI_MakeWire(e4, e5, e6);
        const TopoDS_Face f4 = BRepBuilderAPI_MakeFace(w4);

        BRepOffsetAPI_Sewing sew;
        sew.Add(f1);
        sew.Add(f2);
        sew.Add(f3);
        sew.Add(f4);

        sew.Perform();
        const TopoDS_Shell s = TopoDS::Shell(sew.SewedShape());

        return BRepBuilderAPI_MakeSolid(s);
    }
};

#endif
