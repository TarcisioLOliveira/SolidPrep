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
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS.hxx>
#include <IntTools_EdgeEdge.hxx>
#include <Standard_Handle.hxx>
#include <vector>

#include "element.hpp"
#include "utils.hpp"

template<class T>
class MeshElementCommon : public MeshElement{
    public:
    virtual double get_stress_at(gp_Pnt point, const std::vector<double>& u) const override{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        const std::vector<double> DB = this->get_DB(point);

        std::vector<double> results(N, 0);
        std::vector<double> u_vec(K_DIM, 0);
        for(size_t l = 0; l < N; ++l){
            for(size_t j = 0; j < NODE_DOF; ++j){
                if(this->nodes[l]->u_pos[j] > -1){
                    u_vec[l*NODE_DOF + j] = u[this->nodes[l]->u_pos[j]];
                }
            }
        }
        for(size_t i = 0; i < N; ++i){
            for(size_t j = 0; j < K_DIM; ++j){
                results[i] += DB[i*K_DIM + j]*u_vec[j];
            }
        }

        return std::sqrt(std::pow(results[0], 2) - results[0]*results[1] + std::pow(results[1], 2) + 3*std::pow(results[2], 2));
    }

    virtual std::vector<double> get_internal_loads(const std::vector<double>& u) const override{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        const std::vector<double> k = this->get_k();

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

    virtual double get_compliance(const std::vector<double>& u) const override{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        const auto k = this->get_k();

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
    virtual double get_compliance(const std::vector<double>& u, const std::vector<double>& l) const override{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        const auto k = this->get_k();

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
        std::vector<double> l_vec(9, 0);
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


    protected:
    MeshElementCommon(const std::vector<MeshNode*>& nodes, Material* m):
            MeshElement(nodes, m)
            {}

};


template<class T>
class MeshElementCommon2D : public MeshElementCommon<T>{
    public:
    static const size_t S_SIZE = 3; // Size of the stress and strain vectors
    static const size_t DIM    = 2; // Number of dimensions

    static const utils::ProblemType PROBLEM_TYPE = utils::PROBLEM_TYPE_2D;

    virtual std::vector<double> get_stress_tensor(const gp_Pnt& p, const std::vector<double>& u) const override{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        const std::vector<double> DB = this->get_DB(p);

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

        const std::vector<double> S{results[0], results[2],
                                    results[2], results[1]};

        return S;
    }

    virtual void get_virtual_load(double mult, const gp_Pnt& point, const std::vector<double>& u, std::vector<double>& l) const override{
        const size_t N = T::NODES_PER_ELEM;
        const size_t K_DIM = T::K_DIM;
        const size_t NODE_DOF = T::NODE_DOF;

        const std::vector<double> DB = this->get_DB(point);
        const std::vector<double> V{1, -0.5, 0,
                                   -0.5, 1, 0,
                                   0,   0, 3};

        std::vector<double> u_vec(K_DIM, 0);
        for(size_t i = 0; i < N; ++i){
            for(size_t j = 0; j < NODE_DOF; ++j){
                if(this->nodes[i]->u_pos[j] > -1){
                    u_vec[i*NODE_DOF+j] = u[this->nodes[i]->u_pos[j]];
                }
            }
        }

        std::vector<double> vec(N, 0);

        for(size_t i = 0; i < S_SIZE; ++i){
            for(size_t j = 0; j < K_DIM; ++j){
                vec[i] += DB[i*K_DIM + j]*u_vec[j];
            }
        }
        std::vector<double> vec2(N, 0);
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

    virtual std::vector<double> get_f(const gp_Dir& dir, double norm, const std::vector<gp_Pnt>& points) const override{
        const size_t K_DIM = T::K_DIM;

        double px = dir.X()*norm;
        double py = dir.Y()*norm;

        auto Nf = this->get_Nf(points);

        std::vector<double> f(K_DIM, 0);
        for(size_t i = 0; i < K_DIM; ++i){
            f[i] = Nf[DIM*i]*px + Nf[DIM*i+1]*py;
        }

        return f;
    }

    protected:
    const double t;
    MeshElementCommon2D(const std::vector<MeshNode*>& nodes, Material* m, double thickness):
            MeshElementCommon<T>(nodes, m), t(thickness)
            {}
};

template<class T>
class MeshElementCommon2DTri : public MeshElementCommon2D<T>{
    public:
    virtual std::vector<gp_Pnt> get_intersection_points(const TopoDS_Shape& crosssection) const override{
        const size_t N = T::NODES_PER_ELEM;
        std::vector<gp_Pnt> points;

        const TopoDS_Edge line_edge = TopoDS::Edge(crosssection);

        for(size_t i = 0; i < N; ++i){
            size_t j = (i + 1) % 3;

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

    virtual double get_volume() const override{
        const gp_Mat deltaM(1, this->nodes[0]->point.X(), this->nodes[0]->point.Y(),
                            1, this->nodes[1]->point.X(), this->nodes[1]->point.Y(),
                            1, this->nodes[2]->point.X(), this->nodes[2]->point.Y());

        return 0.5*std::abs(deltaM.Determinant())*this->t;
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
    MeshElementCommon2DTri(const std::vector<MeshNode*>& nodes, Material* m, double thickness):
            MeshElementCommon2D<T>(nodes, m, thickness)
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

#endif