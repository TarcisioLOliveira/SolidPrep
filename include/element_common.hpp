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
#include "math/matrix.hpp"
#include "utils.hpp"
#include "utils/gauss_legendre.hpp"
#include <lapacke.h>

template<class T>
class MeshElementCommon : public MeshElement{
    public:
    virtual ~MeshElementCommon() = default;

    virtual math::Vector get_internal_loads(const math::Matrix& D, const double t, const std::vector<double>& u) const override{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        const math::Matrix k = this->get_k(D, t);

        math::Vector U(K_DIM);
        for(size_t l = 0; l < N; ++l){
            for(size_t j = 0; j < NODE_DOF; ++j){
                U[l*NODE_DOF + j] = u[this->nodes[l]->u_pos[j]];
            }
        }

        return k*U;
    }
    
    virtual double get_compliance(const math::Matrix& D, const double t, const std::vector<double>& u) const override{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        const auto k = this->get_k(D, t);

        math::Vector u_vec(K_DIM);
        for(size_t i = 0; i < N; ++i){
            for(size_t j = 0; j < NODE_DOF; ++j){
                u_vec[i*NODE_DOF+j] = u[this->nodes[i]->u_pos[j]];
            }
        }
        
        return u_vec.T()*k*u_vec;
    }
    virtual double get_compliance(const math::Matrix& D, const double t, const std::vector<double>& u, const std::vector<double>& l) const override{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        const auto k = this->get_k(D, t);

        math::Vector u_vec(K_DIM, 0);
        for(size_t i = 0; i < N; ++i){
            for(size_t j = 0; j < NODE_DOF; ++j){
                u_vec[i*NODE_DOF+j] = u[this->nodes[i]->u_pos[j]];
            }
        }
        math::Vector l_vec(K_DIM, 0);
        for(size_t i = 0; i < N; ++i){
            for(size_t j = 0; j < NODE_DOF; ++j){
                l_vec[i*NODE_DOF+j] = l[this->nodes[i]->u_pos[j]];
            }
        }
        return u_vec.T()*k*l_vec;
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

    virtual double _von_Mises_derivative(const math::Matrix& D, const math::Matrix& dD, double mult, const gp_Pnt& point, const std::vector<double>& u, const math::Matrix& V) const{
        auto s = this->_get_stress_vector(D, point, u);
        auto ds = this->_get_stress_vector(dD, point, u);

        return mult*(s.T()*V*ds);
    }

    virtual double _von_Mises_derivative_sh(const math::Matrix& D, double mult, const gp_Pnt& point, const std::vector<double>& u, const size_t n, const size_t dof, const math::Matrix& V) const{
        auto s = this->_get_stress_vector(D, point, u);
        auto ds = this->_get_stress_vector_sh(D, point, u, n, dof);

        return mult*(s.T()*V*ds);
    }

    virtual void _get_virtual_load(const math::Matrix& D, double mult, const gp_Pnt& p, const std::vector<double>& u, std::vector<double>& l, const math::Matrix& V) const{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        const math::Matrix DB = D*this->get_B(p);

        math::Vector u_vec(K_DIM, 0);
        for(size_t i = 0; i < N; ++i){
            for(size_t j = 0; j < NODE_DOF; ++j){
                if(this->nodes[i]->u_pos[j] > -1){
                    u_vec[i*NODE_DOF+j] = u[this->nodes[i]->u_pos[j]];
                }
            }
        }

        const math::Vector result = mult*(DB.T()*V*DB*u_vec);

        for(size_t k = 0; k < N; ++k){
            for(size_t j = 0; j < NODE_DOF; ++j){
                if(this->nodes[k]->u_pos[j] > -1){
                    l[this->nodes[k]->u_pos[j]] += mult*result[k*NODE_DOF+j];
                }
            }
        }
    }

    virtual math::Vector _get_stress_vector(const math::Matrix& D, const gp_Pnt& p, const std::vector<double>& u) const{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        math::Vector U(K_DIM);
        for(size_t l = 0; l < N; ++l){
            for(size_t j = 0; j < NODE_DOF; ++j){
                if(this->nodes[l]->u_pos[j] > -1){
                    U[l*NODE_DOF + j] = u[this->nodes[l]->u_pos[j]];
                }
            }
        }

        return D*this->get_B(p)*U;
    }
    virtual math::Vector _get_stress_vector_sh(const math::Matrix& D, const gp_Pnt& p, const std::vector<double>& u, const size_t n, const size_t dof) const{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        math::Vector U(K_DIM);
        for(size_t l = 0; l < N; ++l){
            for(size_t j = 0; j < NODE_DOF; ++j){
                if(this->nodes[l]->u_pos[j] > -1){
                    U[l*NODE_DOF + j] = u[this->nodes[l]->u_pos[j]];
                }
            }
        }

        return D*this->get_dB_sh(p, n, dof)*U;
    }


    virtual math::Matrix _get_stress_tensor(const math::Matrix& D, const gp_Pnt& p, const std::vector<double>& u, const std::vector<size_t> indices) const{

        const auto results = this->_get_stress_vector(D, p, u);

        const size_t N = std::sqrt(indices.size());

        math::Matrix S(N, N);
        for(size_t i = 0; i < indices.size(); ++i){
            S.data()[i] = results[indices[i]];
        }

        return S;
    }

    virtual math::Vector _get_strain_vector(const gp_Pnt& p, const std::vector<double>& u) const{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        math::Vector U(K_DIM, 0);
        for(size_t l = 0; l < N; ++l){
            for(size_t j = 0; j < NODE_DOF; ++j){
                if(this->nodes[l]->u_pos[j] > -1){
                    U[l*NODE_DOF + j] = u[this->nodes[l]->u_pos[j]];
                }
            }
        }

        return this->get_B(p)*U;
    }


    virtual math::Matrix _get_strain_tensor(const gp_Pnt& p, const std::vector<double>& u, const std::vector<size_t> indices) const{

        auto results = this->_get_strain_vector(p, u);

        const size_t N = std::sqrt(indices.size());

        math::Matrix S(N, N);
        for(size_t i = 0; i < indices.size(); ++i){
            S.data()[i] = results[indices[i]];
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

    virtual double get_stress_at(const math::Matrix& D, const gp_Pnt& point, const std::vector<double>& u, const double eps = 0) const override{
        auto results = this->_get_stress_vector(D, point, u);

        return std::sqrt(results[0]*results[0] - results[0]*results[1] + results[1]*results[1] + 3*results[2]*results[2] + eps);
    }
    virtual double get_strain_VM(const gp_Pnt& p, const std::vector<double>& u, const double eps = 0) const override{
        auto results = this->_get_strain_vector(p, u);

        return std::sqrt(results[0]*results[0] - results[0]*results[1] + results[1]*results[1] + 3*results[2]*results[2] + eps);
    }

    virtual math::Matrix get_stress_tensor(const math::Matrix& D, const gp_Pnt& p, const std::vector<double>& u) const override{
        const std::vector<size_t> indices{0, 2, 2, 1};
        return this->_get_stress_tensor(D, p, u, indices);
    }

    virtual math::Matrix get_strain_tensor(const gp_Pnt& p, const std::vector<double>& u) const override{
        const std::vector<size_t> indices{0, 2, 2, 1};
        return this->_get_strain_tensor(p, u, indices);
    }
    virtual math::Vector get_principal_strains(const gp_Pnt& p, const std::vector<double>& u) const override{
        const auto s = this->get_strain_tensor(p, u);
        math::Eigen e(s);

        return e.eigenvalues;
    }

    virtual math::Vector get_strain_vector(const gp_Pnt& p, const std::vector<double>& u) const override{
        return this->_get_strain_vector(p, u);
    }

    virtual void get_virtual_load(const math::Matrix& D, double mult, const gp_Pnt& point, const std::vector<double>& u, std::vector<double>& l) const override{
        const math::Matrix V({1, -0.5, 0,
                              -0.5, 1, 0,
                              0,   0, 3}, 3, 3);

        this->_get_virtual_load(D, mult, point, u, l, V);
    }

    virtual double von_Mises_derivative(const math::Matrix& D, const math::Matrix& dD, double mult, const gp_Pnt& point, const std::vector<double>& u) const override{
        const math::Matrix V({1, -0.5, 0,
                              -0.5, 1, 0,
                              0,   0, 3}, 3, 3);

        return this->_von_Mises_derivative(D, dD, mult, point, u, V);
    }
    virtual double von_Mises_derivative_sh(const math::Matrix& D, double mult, const gp_Pnt& point, const std::vector<double>& u, const size_t n, const size_t dof) const override{
        const math::Matrix V({1, -0.5, 0,
                              -0.5, 1, 0,
                              0,   0, 3}, 3, 3);

        return this->_von_Mises_derivative_sh(D, mult, point, u, n, dof, V);
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

    virtual math::Vector get_f(const double t, const math::Vector& vec, const std::vector<gp_Pnt>& points) const override{
        return this->get_Nf(t, points)*vec;
    }

    virtual math::Matrix get_uu(const MeshElement* const e2, const std::vector<gp_Pnt>& bounds, const gp_Dir n) const override{
        const size_t ORDER = T::ORDER;
        const size_t KW = T::K_DIM;
        const size_t DIM = this->DIM;

        math::Vector NN(2*KW, 0);
        math::Matrix MnMn(2*KW, 2*KW, 0);

        constexpr size_t GN = 2*ORDER;
        const auto& GL = utils::GaussLegendre<GN>::get();

        const gp_Pnt p1 = bounds[0];
        const gp_Pnt p2 = bounds[1];

        for(auto xi = GL.begin(); xi < GL.end(); ++xi){
            const gp_Pnt pi(
                p1.X() + (p2.X() - p1.X())*xi->x,
                p1.Y() + (p2.Y() - p1.Y())*xi->x,
                p1.Z() + (p2.Z() - p1.Z())*xi->x
            );

            const auto N1 = this->get_Ni(pi);
            const auto N2 = e2->get_Ni(pi);
            NN.fill(0);
            for(size_t i = 0; i < KW; ++i){
                for(size_t j = 0; j < DIM; ++j){
                    NN[i] += N1(j, i)*n.Coord(1+j);
                    NN[i + KW] += N2(j, i)*n.Coord(1+j);
                }
            }
            MnMn += (xi->w)*(NN*NN.T());
        }
        for(size_t i = 0; i < KW; ++i){
            for(size_t j = KW; j < 2*KW; ++j){
                MnMn(i, j) *= -1.0;
                MnMn(j, i) *= -1.0;
            }
        }
        const double drnorm = bounds[0].Distance(bounds[1]);
        MnMn *= drnorm;

        return MnMn;
    }

    virtual math::Matrix get_MnMn(const MeshElement* const e2, const std::vector<double>& u_ext, const std::vector<gp_Pnt>& bounds, const gp_Dir n) const override{
        const size_t ORDER = T::ORDER;
        const size_t NODES_PER_ELEM = T::NODES_PER_ELEM;
        const size_t DOF = T::NODE_DOF;
        const size_t KW = T::K_DIM;
        const size_t DIM = this->DIM;

        math::Vector uv1(KW), uv2(KW);
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            for(size_t j = 0; j < DOF; ++j){
                uv1[i*DOF + j] = u_ext[this->nodes[i]->u_pos[j]];
                uv2[i*DOF + j] = u_ext[e2->nodes[i]->u_pos[j]];
            }
        }
        math::Vector NN(2*KW, 0);
        math::Matrix MnMn(2*KW, 2*KW, 0);

        constexpr size_t GN = 2*ORDER + 2;
        const auto& GL = utils::GaussLegendre<GN>::get();

        const gp_Pnt p1 = bounds[0];
        const gp_Pnt p2 = bounds[1];

        for(auto xi = GL.begin(); xi < GL.end(); ++xi){

            const gp_Pnt pi(
                p1.X() + (p2.X() - p1.X())*xi->x,
                p1.Y() + (p2.Y() - p1.Y())*xi->x,
                p1.Z() + (p2.Z() - p1.Z())*xi->x
            );

            const auto N1 = this->get_Ni(pi);
            const auto N2 = e2->get_Ni(pi);
            double gp = 0;
            math::Vector up1(N1*uv1);
            math::Vector up2(N2*uv2);
            for(size_t j = 0; j < DIM; ++j){
                gp += (up2[j] - up1[j])*n.Coord(1+j);
            }
            if(gp < 1e-7){
                NN.fill(0);
                for(size_t i = 0; i < KW; ++i){
                    for(size_t j = 0; j < DIM; ++j){
                        NN[i] += N1(j, i)*n.Coord(1+j);
                        NN[i + KW] += N2(j, i)*n.Coord(1+j);
                    }
                }
                MnMn += (xi->w)*(NN*NN.T());
            }
        }
        for(size_t i = 0; i < KW; ++i){
            for(size_t j = KW; j < 2*KW; ++j){
                MnMn(i, j) *= -1.0;
                MnMn(j, i) *= -1.0;
            }
        }
        const double drnorm = bounds[0].Distance(bounds[1]);
        MnMn *= drnorm;

        return MnMn;
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

    virtual double get_stress_at(const math::Matrix& D, const gp_Pnt& point, const std::vector<double>& u, const double eps = 0) const override{
        auto s = this->_get_stress_vector(D, point, u);

        const double ds1 = s[0] - s[1];
        const double ds2 = s[1] - s[2];
        const double ds3 = s[2] - s[0];

        return std::sqrt(0.5*(ds1*ds1 + ds2*ds2 + ds3*ds3 + 6*(s[3]*s[3] + s[4]*s[4] + s[5]*s[5])) + eps);
    }
    virtual double get_strain_VM(const gp_Pnt& p, const std::vector<double>& u, const double eps = 0) const override{
        auto s = this->_get_strain_vector(p, u);

        const double ds1 = s[0] - s[1];
        const double ds2 = s[1] - s[2];
        const double ds3 = s[2] - s[0];

        return std::sqrt(0.5*(ds1*ds1 + ds2*ds2 + ds3*ds3 + 6*(s[3]*s[3] + s[4]*s[4] + s[5]*s[5])) + eps);
    }
    virtual math::Vector get_principal_strains(const gp_Pnt& p, const std::vector<double>& u) const override{
        auto s = this->get_strain_tensor(p, u);
        math::Eigen e(s);

        return e.eigenvalues;
    }

    virtual math::Matrix get_stress_tensor(const math::Matrix& D, const gp_Pnt& p, const std::vector<double>& u) const override{
        const std::vector<size_t> indices{0, 3, 4, 3, 1, 5, 4, 5, 2};
        return this->_get_stress_tensor(D, p, u, indices);
    }

    virtual math::Matrix get_strain_tensor(const gp_Pnt& p, const std::vector<double>& u) const override{
        const std::vector<size_t> indices{0, 3, 4, 3, 1, 5, 4, 5, 2};
        return this->_get_strain_tensor(p, u, indices);
    }

    virtual math::Vector get_strain_vector(const gp_Pnt& p, const std::vector<double>& u) const override{
        return this->_get_strain_vector(p, u);
    }

    virtual void get_virtual_load(const math::Matrix& D, double mult, const gp_Pnt& point, const std::vector<double>& u, std::vector<double>& l) const override{
        const math::Matrix V({   1, -0.5, -0.5, 0, 0, 0,
                              -0.5,    1, -0.5, 0, 0, 0,
                              -0.5, -0.5,    1, 0, 0, 0,
                                 0,    0,    0, 3, 0, 0,
                                 0,    0,    0, 0, 3, 0,
                                 0,    0,    0, 0, 0, 3}, 6, 6);

        this->_get_virtual_load(D, mult, point, u, l, V);
    }

    virtual double von_Mises_derivative(const math::Matrix& D, const math::Matrix& dD, double mult, const gp_Pnt& point, const std::vector<double>& u) const override{
        const math::Matrix V({   1, -0.5, -0.5, 0, 0, 0,
                              -0.5,    1, -0.5, 0, 0, 0,
                              -0.5, -0.5,    1, 0, 0, 0,
                                 0,    0,    0, 3, 0, 0,
                                 0,    0,    0, 0, 3, 0,
                                 0,    0,    0, 0, 0, 3}, 6, 6);

        return this->_von_Mises_derivative(D, dD, mult, point, u, V);
    }
    virtual double von_Mises_derivative_sh(const math::Matrix& D, double mult, const gp_Pnt& point, const std::vector<double>& u, const size_t n, const size_t dof) const override{
        const math::Matrix V({   1, -0.5, -0.5, 0, 0, 0,
                              -0.5,    1, -0.5, 0, 0, 0,
                              -0.5, -0.5,    1, 0, 0, 0,
                                 0,    0,    0, 3, 0, 0,
                                 0,    0,    0, 0, 3, 0,
                                 0,    0,    0, 0, 0, 3}, 6, 6);

        return this->_von_Mises_derivative_sh(D, mult, point, u, n, dof, V);
    }

    virtual math::Vector get_f(const double t, const math::Vector& vec, const std::vector<gp_Pnt>& points) const override{
        return this->get_Nf(t, points)*vec;
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

    virtual double get_dV_sh(const double t, const size_t n, const size_t dof) const override{
        (void)t;
        const gp_Pnt p0 = this->nodes[0]->point;
        const gp_Pnt p1 = this->nodes[1]->point;
        const gp_Pnt p2 = this->nodes[2]->point;
        const gp_Pnt p3 = this->nodes[3]->point;

        const gp_Vec v0(p0, p1);
        const gp_Vec v1(p0, p2);
        const gp_Vec v2(p0, p3);

        std::vector<gp_Pnt> dp(4, gp_Pnt(0,0,0));
        dp[n].SetCoord(1+dof, 1);

        const double mult = (v0.DotCross(v1, v2) > 0) ? 1 : -1;

        const gp_Vec dv0(dp[0], dp[1]);
        const gp_Vec dv1(dp[0], dp[2]);
        const gp_Vec dv2(dp[0], dp[3]);

        return mult*(dv0.DotCross(v1, v2) + v0.Dot(dv1.Crossed(v2) + v1.Crossed(dv2)))/6;
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

    virtual math::Matrix get_uu(const MeshElement* const e2, const std::vector<gp_Pnt>& bounds, const gp_Dir n) const override{
        const size_t ORDER = T::ORDER;
        const size_t KW = T::K_DIM;
        const size_t DIM = this->DIM;

        const auto& gli = utils::GaussLegendreTri<2*ORDER>::get();
        math::Vector NN(2*KW, 0);
        math::Matrix MnMn(2*KW, 2*KW, 0);

        for(auto it = gli.begin(); it != gli.end(); ++it){
            const gp_Pnt pi = MeshElementCommon3DTet<T>::ECTRI_GL_to_point(*it, bounds);
            const auto N1 = this->get_Ni(pi);
            const auto N2 = e2->get_Ni(pi);
            NN.fill(0);
            for(size_t i = 0; i < KW; ++i){
                for(size_t j = 0; j < DIM; ++j){
                    NN[i] += N1(j, i)*n.Coord(1+j);
                    NN[i + KW] += N2(j, i)*n.Coord(1+j);
                }
            }
            MnMn += (it->w)*(NN*NN.T());
        }
        for(size_t i = 0; i < KW; ++i){
            for(size_t j = KW; j < 2*KW; ++j){
                MnMn(i, j) *= -1.0;
                MnMn(j, i) *= -1.0;
            }
        }
        const gp_Vec v1(bounds[0], bounds[1]);
        const gp_Vec v2(bounds[0], bounds[2]);
        const double A = 0.5*v1.Crossed(v2).Magnitude();
        MnMn *= A;

        return MnMn;
    }
    virtual math::Matrix get_MnMn(const MeshElement* const e2, const std::vector<double>& u_ext, const std::vector<gp_Pnt>& bounds, const gp_Dir n) const override{
        const size_t ORDER = T::ORDER;
        const size_t NODES_PER_ELEM = T::NODES_PER_ELEM;
        const size_t DOF = T::NODE_DOF;
        const size_t KW = T::K_DIM;
        const size_t DIM = this->DIM;

        math::Vector uv1(KW), uv2(KW);
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            for(size_t j = 0; j < DOF; ++j){
                uv1[i*DOF + j] = u_ext[this->nodes[i]->u_pos[j]];
                uv2[i*DOF + j] = u_ext[e2->nodes[i]->u_pos[j]];
            }
        }
        const auto& gli = utils::GaussLegendreTri<2*ORDER + 2>::get();
        math::Vector NN(2*KW, 0);
        math::Matrix MnMn(2*KW, 2*KW, 0);

        for(auto it = gli.begin(); it != gli.end(); ++it){
            const gp_Pnt pi = MeshElementCommon3DTet<T>::ECTRI_GL_to_point(*it, bounds);
            const auto N1 = this->get_Ni(pi);
            const auto N2 = e2->get_Ni(pi);
            double gp = 0;
            math::Vector up1(N1*uv1);
            math::Vector up2(N2*uv2);
            for(size_t j = 0; j < DIM; ++j){
                gp += (up2[j] - up1[j])*n.Coord(1+j);
            }
            if(gp < 1e-7){
                NN.fill(0);
                for(size_t i = 0; i < KW; ++i){
                    for(size_t j = 0; j < DIM; ++j){
                        NN[i] += N1(j, i)*n.Coord(1+j);
                        NN[i + KW] += N2(j, i)*n.Coord(1+j);
                    }
                }
                MnMn += (it->w)*(NN*NN.T());
            }
        }
        for(size_t i = 0; i < KW; ++i){
            for(size_t j = KW; j < 2*KW; ++j){
                MnMn(i, j) *= -1.0;
                MnMn(j, i) *= -1.0;
            }
        }
        const gp_Vec v1(bounds[0], bounds[1]);
        const gp_Vec v2(bounds[0], bounds[2]);
        const double A = 0.5*v1.Crossed(v2).Magnitude();
        MnMn *= A;

        return MnMn;
    }

    protected:
    MeshElementCommon3DTet(const std::vector<MeshNode*>& nodes):
            MeshElementCommon3D<T>(nodes)
            {}

    inline static gp_Pnt ECTRI_GL_to_point(const utils::GLPointTri& gp, const std::vector<gp_Pnt>& bounds){
        return gp_Pnt
            (gp.a*bounds[0].X() + gp.b*bounds[1].X() + gp.c*bounds[2].X(),
             gp.a*bounds[0].Y() + gp.b*bounds[1].Y() + gp.c*bounds[2].Y(),
             gp.a*bounds[0].Z() + gp.b*bounds[1].Z() + gp.c*bounds[2].Z());
    } 

    inline TopoDS_Solid generate_geometry(const gp_Pnt& p1, const gp_Pnt& p2, const gp_Pnt& p3, const gp_Pnt& p4) const{
        const TopoDS_Vertex v1 = BRepBuilderAPI_MakeVertex(p1);
        const TopoDS_Vertex v2 = BRepBuilderAPI_MakeVertex(p2);
        const TopoDS_Vertex v3 = BRepBuilderAPI_MakeVertex(p3);
        const TopoDS_Vertex v4 = BRepBuilderAPI_MakeVertex(p4);

        const TopoDS_Edge e1 = BRepBuilderAPI_MakeEdge(v1, v2);
        const TopoDS_Edge e2 = BRepBuilderAPI_MakeEdge(v1, v3);
        const TopoDS_Edge e3 = BRepBuilderAPI_MakeEdge(v1, v4);

        const TopoDS_Edge e12 = BRepBuilderAPI_MakeEdge(v2, v1);
        const TopoDS_Edge e22 = BRepBuilderAPI_MakeEdge(v3, v1);
        const TopoDS_Edge e32 = BRepBuilderAPI_MakeEdge(v4, v1);

        const TopoDS_Edge e4 = BRepBuilderAPI_MakeEdge(v2, v3);
        const TopoDS_Edge e5 = BRepBuilderAPI_MakeEdge(v3, v4);
        const TopoDS_Edge e6 = BRepBuilderAPI_MakeEdge(v4, v2);

        const TopoDS_Wire w1 = BRepBuilderAPI_MakeWire(e1, e4, e22);
        const TopoDS_Face f1 = BRepBuilderAPI_MakeFace(w1);

        const TopoDS_Wire w2 = BRepBuilderAPI_MakeWire(e2, e5, e32);
        const TopoDS_Face f2 = BRepBuilderAPI_MakeFace(w2);

        const TopoDS_Wire w3 = BRepBuilderAPI_MakeWire(e3, e6, e12);
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

        TopoDS_Solid so = BRepBuilderAPI_MakeSolid(s);

        return so;
    }
};

template<class T>
class MeshElementCommon3DHex : public MeshElementCommon3D<T>{
    public:
    virtual ~MeshElementCommon3DHex() = default;

    static const Element::Shape SHAPE_TYPE = Element::Shape::QUAD;

    virtual double get_volume(const double t) const override{
        (void)t;
        // 0 1 3 5
        // 1 3 2 5
        // 2 3 6 5
        // 3 7 6 5
        // 3 7 4 5
        // 0 3 4 5

        double V = 0;
        V += this->get_volume_tet(
            this->nodes[0]->point,
            this->nodes[1]->point,
            this->nodes[3]->point,
            this->nodes[5]->point
        );
        V += this->get_volume_tet(
            this->nodes[1]->point,
            this->nodes[3]->point,
            this->nodes[2]->point,
            this->nodes[5]->point
        );
        V += this->get_volume_tet(
            this->nodes[2]->point,
            this->nodes[3]->point,
            this->nodes[6]->point,
            this->nodes[5]->point
        );
        V += this->get_volume_tet(
            this->nodes[3]->point,
            this->nodes[7]->point,
            this->nodes[6]->point,
            this->nodes[5]->point
        );
        V += this->get_volume_tet(
            this->nodes[3]->point,
            this->nodes[7]->point,
            this->nodes[4]->point,
            this->nodes[5]->point
        );
        V += this->get_volume_tet(
            this->nodes[0]->point,
            this->nodes[3]->point,
            this->nodes[4]->point,
            this->nodes[5]->point
        );

        return V;
    }

    virtual TopoDS_Shape get_shape() const override{
        std::array<gp_Pnt, 8> p;
        for(size_t i = 0; i < 8; ++i){
            p[i] = this->nodes[i]->point;
        }

        return this->generate_geometry(p);
    }

    virtual TopoDS_Shape get_shape(const std::vector<gp_Vec>& disp) const override{
        std::array<gp_Pnt, 8> p;
        for(size_t i = 0; i < 8; ++i){
            p[i] = this->nodes[i]->point.Translated(disp[i]);
        }

        return this->generate_geometry(p);
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
    MeshElementCommon3DHex(const std::vector<MeshNode*>& nodes):
            MeshElementCommon3D<T>(nodes)
            {}

    inline double get_volume_tet(const gp_Pnt& p0, const gp_Pnt& p1, const gp_Pnt& p2, const gp_Pnt& p3) const{
        const gp_Vec v0(p0, p1);
        const gp_Vec v1(p0, p2);
        const gp_Vec v2(p0, p3);

        return std::abs(v0.DotCross(v1, v2))/6;
    }

    inline TopoDS_Face get_quad(const TopoDS_Vertex& v1, const TopoDS_Vertex& v2, const TopoDS_Vertex& v3, const TopoDS_Vertex& v4) const{
        const TopoDS_Edge e1 = BRepBuilderAPI_MakeEdge(v1, v2);
        const TopoDS_Edge e2 = BRepBuilderAPI_MakeEdge(v2, v3);
        const TopoDS_Edge e3 = BRepBuilderAPI_MakeEdge(v3, v4);
        const TopoDS_Edge e4 = BRepBuilderAPI_MakeEdge(v4, v1);

        const TopoDS_Wire w = BRepBuilderAPI_MakeWire(e1, e2, e3, e4);

        const TopoDS_Face f = BRepBuilderAPI_MakeFace(w);

        return f;
    }

    inline TopoDS_Solid generate_geometry(const std::array<gp_Pnt, 8>& p) const{
        std::array<TopoDS_Vertex, 8> v;
        for(size_t i = 0; i < 8; ++i){
            v[i] = BRepBuilderAPI_MakeVertex(p[i]);
        }

        BRepOffsetAPI_Sewing sew;
        TopoDS_Face f = this->get_quad(v[0], v[1], v[2], v[3]);
        sew.Add(f);
        f = this->get_quad(v[1], v[2], v[6], v[5]);
        sew.Add(f);
        f = this->get_quad(v[4], v[5], v[6], v[7]);
        sew.Add(f);
        f = this->get_quad(v[0], v[4], v[7], v[3]);
        sew.Add(f);
        f = this->get_quad(v[0], v[4], v[5], v[1]);
        sew.Add(f);
        f = this->get_quad(v[3], v[7], v[6], v[2]);
        sew.Add(f);

        sew.Perform();
        const TopoDS_Shell s = TopoDS::Shell(sew.SewedShape());

        TopoDS_Solid so = BRepBuilderAPI_MakeSolid(s);

        return so;
    }
};
#endif
