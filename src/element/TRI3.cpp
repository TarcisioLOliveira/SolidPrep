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

#include "element/TRI3.hpp"
#include "cblas.h"
#include "logger.hpp"
#include "project_data.hpp"
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS.hxx>
#include <IntTools_EdgeEdge.hxx>

namespace element{

TRI3::TRI3(ElementShape s, ProjectData* data):
    MeshElement(s.nodes), mat(data->material.get()), t(data->thickness){}

std::vector<double> TRI3::get_k() const{
    size_t N = this->nodes.size();

    std::vector<double> B(3*2*N, 0);

    std::vector<gp_Pnt> p;
    for(auto n:this->nodes){
        p.push_back(n->point);
    }

    gp_Mat deltaM(1, p[0].X(), p[0].Y(), 1, p[1].X(), p[1].Y(), 1, p[2].X(), p[2].Y());

    double Delta = 0.5*std::abs(deltaM.Determinant());

    std::vector<double> a, b, c;
    for(size_t i = 0; i < N; ++i){
        size_t j = (i + 1) % 3;
        size_t k = (i + 2) % 3;

        a.push_back(p[j].X()*p[k].Y() - p[k].X()*p[j].Y());
        b.push_back(p[j].Y() - p[k].Y());
        c.push_back(p[k].X() - p[j].X());
    }

    for(size_t i = 0; i < N; ++i){
        B[i*2 + 0*2*N] = b[i]/(2*Delta);
        B[i*2 + 1*2*N] = 0;
        B[i*2 + 2*2*N] = c[i]/(2*Delta);
        B[i*2 + 0*2*N + 1] = 0;
        B[i*2 + 1*2*N + 1] = c[i]/(2*Delta);
        B[i*2 + 2*2*N + 1] = b[i]/(2*Delta);
    }

    std::vector<double> DB(3*2*N, 0);
    auto D = this->mat->stiffness_2D();

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 2*N, 3, 1, D.data(), 3, B.data(), 2*N, 0, DB.data(), 2*N);

    std::vector<double> K(2*N*2*N, 0);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 2*N, 2*N, 3, 1, B.data(), 2*N, DB.data(), 2*N, 0, K.data(), 2*N);

    cblas_dscal(K.size(), this->t*Delta, K.data(), 1);

    return K;
}

MeshNode* TRI3::get_stresses(size_t node, const std::vector<double>& u, double density) const{
    size_t N = this->nodes.size();

    auto DB = this->get_DB(this->nodes[node]->point);

    MeshNode2D* n = static_cast<MeshNode2D*>(this->nodes[node]);
    for(size_t i = 0; i < 3; ++i){
        n->results[i] = 0;
        for(size_t l = 0; l < 3; ++l){
            for(size_t j = 0; j < 2; ++j){
                if(this->nodes[l]->u_pos[j] > -1){
                    n->results[i] += DB[2*N*i + l*2 + j]*u[this->nodes[l]->u_pos[j]];
                }
            }
        }
        n->results[i] = density*std::abs(n->results[i]);
    }

    return this->get_node(node);
}

double TRI3::get_stress_at(gp_Pnt point, const std::vector<double>& u) const{
    size_t N = this->nodes.size();

    auto DB = this->get_DB(point);

    std::vector<double> results(3, 0);
    for(size_t i = 0; i < 3; ++i){
        for(size_t l = 0; l < 3; ++l){
            for(size_t j = 0; j < 2; ++j){
                if(this->nodes[l]->u_pos[j] > -1){
                    results[i] += DB[2*N*i + l*2 + j]*u[this->nodes[l]->u_pos[j]];
                }
            }
        }
    }

    return std::sqrt(std::pow(results[0], 2) - results[0]*results[1] + std::pow(results[1], 2) + 3*std::pow(results[2], 2));
}

MeshNode* TRI3::get_internal_loads(size_t node, const std::vector<double>& u) const{
    logger::log_assert(node >= 0 && node <= 3, logger::ERROR, "wrong value for BeamLinear2D node, must be either 0 or 1.");

    std::vector<double> k = this->get_k();

    MeshNode2D* n = static_cast<MeshNode2D*>(this->nodes[node]);
    for(int i = 0; i < 3; ++i){
        for(size_t l = 0; l < 3; ++l){
            for(int j = 0; j < 2; ++j){
                if(this->nodes[l]->u_pos[j] > -1){
                    n->results[i] += k[(node*2+i)*2*3+l*2+j]*u[this->nodes[l]->u_pos[j]];
                }
            }
        }
    }

    return this->get_node(node);
}

double TRI3::get_compliance(const std::vector<double>& u, const std::vector<double>& l) const{
    auto k = this->get_k();
    std::vector<double> u_vec(6, 0);
    for(size_t k = 0; k < 3; ++k){
        for(int j = 0; j < 2; ++j){
            if(this->nodes[k]->u_pos[j] > -1){
                u_vec[k*2+j] = u[this->nodes[k]->u_pos[j]];
            }
        }
    }

    std::vector<double> f_vec(6, 0);
    if(l.size() > 0){
        std::vector<double> l_vec(6, 0);
        for(size_t k = 0; k < 3; ++k){
            for(int j = 0; j < 2; ++j){
                if(this->nodes[k]->u_pos[j] > -1){
                    l_vec[k*2+j] = l[this->nodes[k]->u_pos[j]];
                }
            }
        }
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 6, 1, 6, 1, k.data(), 6, u_vec.data(), 1, 0, f_vec.data(), 1);
        return cblas_ddot(6, l_vec.data(), 1, f_vec.data(), 1);
    } else {
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 6, 1, 6, 1, k.data(), 6, u_vec.data(), 1, 0, f_vec.data(), 1);
        return cblas_ddot(6, u_vec.data(), 1, f_vec.data(), 1);
    }

}

double TRI3::get_volume() const{
    gp_Mat deltaM(1, this->nodes[0]->point.X(), this->nodes[0]->point.Y(),
                  1, this->nodes[1]->point.X(), this->nodes[1]->point.Y(),
                  1, this->nodes[2]->point.X(), this->nodes[2]->point.Y());

    return 0.5*std::abs(deltaM.Determinant())*this->t;
}

void TRI3::get_virtual_load(double mult, gp_Pnt point, const std::vector<double>& u, std::vector<double>& l) const{
    std::vector<double> DB = this->get_DB(point);
    std::vector<double> V{1, -0.5, 0,
                         -0.5, 1, 0,
                         0,   0, 3};

    std::vector<double> u_vec(6, 0);
    for(size_t k = 0; k < 3; ++k){
        for(int j = 0; j < 2; ++j){
            if(this->nodes[k]->u_pos[j] > -1){
                u_vec[k*2+j] = u[this->nodes[k]->u_pos[j]];
            }
        }
    }

    std::vector<double> f_vec(6, 0);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 1, 6, 1, DB.data(), 6, u_vec.data(), 1, 0, f_vec.data(), 1);

    std::vector<double> res(f_vec);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 1, 3, 1, V.data(), 3, res.data(), 1, 0, f_vec.data(), 1);
    res = f_vec;

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 6, 1, 3, 1, DB.data(), 6, res.data(), 1, 0, f_vec.data(), 1);
    cblas_dscal(6, mult, f_vec.data(), 1);

    for(size_t k = 0; k < 3; ++k){
        for(int j = 0; j < 2; ++j){
            if(this->nodes[k]->u_pos[j] > -1){
                l[this->nodes[k]->u_pos[j]] += f_vec[k*2+j];
            }
        }
    }
}

std::vector<double> TRI3::get_loads_at(gp_Pnt point, const std::vector<double>& u) const{
    std::vector<double> k = this->get_k();
    std::vector<double> f_vec(6,0);
    std::vector<double> u_vec(6,0);

    for(int i = 0; i < 3; ++i){
        for(size_t j = 0; j < 2; ++j){
            if(this->nodes[i]->u_pos[j] > -1){
                u_vec[3*i+j] = u[this->nodes[i]->u_pos[j]];
            }
        }
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 6, 1, 6, 1, k.data(), 6, u_vec.data(), 1, 0, f_vec.data(), 1);

    std::vector<gp_Pnt> p;
    for(auto n:this->nodes){
        p.push_back(n->point);
    }

    gp_Mat deltaM(1, p[0].X(), p[0].Y(), 1, p[1].X(), p[1].Y(), 1, p[2].X(), p[2].Y());

    double delta = 0.5*deltaM.Determinant();

    size_t N = this->nodes.size();
    std::vector<double> a, b, c;
    for(size_t i = 0; i < N; ++i){
        size_t j = (i + 1) % 3;
        size_t k = (i + 2) % 3;

        a.push_back(p[j].X()*p[k].Y() - p[k].X()*p[j].Y());
        b.push_back(p[j].Y() - p[k].Y());
        c.push_back(p[k].X() - p[j].X());
    }

    std::vector<double> L(2);
    L[0] = (a[0] + b[0]*point.X() + c[0]*point.Y())/(2*delta);
    L[1] = (a[1] + b[1]*point.X() + c[1]*point.Y())/(2*delta);
    std::vector<double> Nmat(6*2, 0);
    for(size_t i = 0; i < N; ++i){
        Nmat[2*i] = L[i];
        Nmat[2*i + 6 + 1] = L[i];
    }

    std::vector<double> res(2, 0);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 1, 6, 1, Nmat.data(), 6, f_vec.data(), 1, 0, res.data(), 1);

    return res;
}


std::vector<double> TRI3::get_average_loads(const gp_Pnt& p1, const gp_Pnt& p2, const std::vector<double>& u, const std::vector<size_t> excluded_nodes) const{
    std::vector<double> k = this->get_k();
    std::vector<double> f_vec(6,0);
    std::vector<double> u_vec(6,0);

    for(int i = 0; i < 3; ++i){
        for(size_t j = 0; j < 2; ++j){
            if(this->nodes[i]->u_pos[j] > -1){
                u_vec[3*i+j] = u[this->nodes[i]->u_pos[j]];
            }
        }
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 6, 1, 6, 1, k.data(), 6, u_vec.data(), 1, 0, f_vec.data(), 1);

    std::vector<gp_Pnt> p;
    for(auto n:this->nodes){
        p.push_back(n->point);
    }

    gp_Mat deltaM(1, p[0].X(), p[0].Y(), 1, p[1].X(), p[1].Y(), 1, p[2].X(), p[2].Y());

    double delta = 0.5*deltaM.Determinant();

    size_t N = this->nodes.size();
    std::vector<double> a, b, c;
    for(size_t i = 0; i < N; ++i){
        size_t j = (i + 1) % 3;
        size_t k = (i + 2) % 3;

        a.push_back(p[j].X()*p[k].Y() - p[k].X()*p[j].Y());
        b.push_back(p[j].Y() - p[k].Y());
        c.push_back(p[k].X() - p[j].X());
    }

    double x1 = p1.X();
    double x2 = p2.X();
    double y1 = p1.Y();
    double y2 = p2.Y();


    double dist = p1.Distance(p2);

    std::vector<double> L(3);
    L[0] = (a[0]*(x2-x1)*(y2-y1) + 0.5*b[0]*(x2*x2-x1*x1)*(y2-y1) + 0.5*c[0]*(x2-x1)*(y2*y2-y1*y1))/(2*delta);
    L[1] = (a[1]*(x2-x1)*(y2-y1) + 0.5*b[1]*(x2*x2-x1*x1)*(y2-y1) + 0.5*c[1]*(x2-x1)*(y2*y2-y1*y1))/(2*delta);
    L[2] = (a[2]*(x2-x1)*(y2-y1) + 0.5*b[2]*(x2*x2-x1*x1)*(y2-y1) + 0.5*c[2]*(x2-x1)*(y2*y2-y1*y1))/(2*delta);
    std::vector<double> Nmat(6*2, 0);
    for(size_t i = 0; i < N; ++i){
        Nmat[2*i] = L[i]/dist;
        Nmat[2*i + 6 + 1] = L[i]/dist;
    }

    std::vector<double> res(2, 0);
    std::vector<double> ss_vec(9, 0);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 9, 1, 9, 1, k.data(), 9, u_vec.data(), 1, 0, ss_vec.data(), 1);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 1, 6, 1, Nmat.data(), 6, ss_vec.data(), 1, 0, res.data(), 1);

    return res;

}

std::vector<double> TRI3::get_average_loads_and_torques(const gp_Pnt& p1, const gp_Pnt& p2, const std::vector<double>& u, const std::vector<size_t> excluded_nodes) const{
    // Missing torque. Implemented to avoid compilation errors.
    return this->get_average_loads(p1, p2, u);
}

std::vector<gp_Pnt> TRI3::get_intersection_points(const TopoDS_Shape& crosssection) const{
    size_t N = this->nodes.size();
    std::vector<gp_Pnt> points;

    TopoDS_Edge line_edge = TopoDS::Edge(crosssection);

    for(size_t i = 0; i < N; ++i){
        size_t j = (i + 1) % 3;

        TopoDS_Edge e = BRepBuilderAPI_MakeEdge(this->nodes[i]->point, this->nodes[j]->point);
        gp_Dir dir(gp_Vec(this->nodes[i]->point, this->nodes[j]->point));

        IntTools_EdgeEdge tool(line_edge, e);
        tool.Perform();
        auto common = tool.CommonParts();
        if(common.Size() > 0){
            for(auto& c:common){
                double dist = std::abs(c.VertexParameter2()); // n.normal is a unit vector, and the line starts at 0
                gp_Pnt p = this->nodes[i]->point.Translated(dist*dir);
                points.push_back(p);
            }
        }
    }

    return points;
}

TopoDS_Shape TRI3::get_shape(std::vector<gp_Vec> disp) const{
    gp_Pnt p1 = this->nodes[0]->point;
    gp_Pnt p2 = this->nodes[1]->point;
    gp_Pnt p3 = this->nodes[2]->point;

    if(disp.size() > 0){
        p1.Translate(disp[0]);
        p2.Translate(disp[1]);
        p3.Translate(disp[2]);
    }

    TopoDS_Vertex v1 = BRepBuilderAPI_MakeVertex(this->nodes[0]->point);
    TopoDS_Vertex v2 = BRepBuilderAPI_MakeVertex(this->nodes[1]->point);
    TopoDS_Vertex v3 = BRepBuilderAPI_MakeVertex(this->nodes[2]->point);

    TopoDS_Edge e1 = BRepBuilderAPI_MakeEdge(v1, v2);
    TopoDS_Edge e2 = BRepBuilderAPI_MakeEdge(v2, v3);
    TopoDS_Edge e3 = BRepBuilderAPI_MakeEdge(v3, v1);

    TopoDS_Wire w = BRepBuilderAPI_MakeWire(e1, e2, e3);

    TopoDS_Face f = BRepBuilderAPI_MakeFace(w);

    return f;
}

gp_Pnt TRI3::get_centroid() const{
    double x = 0;
    double y = 0;
    for(auto& n : this->nodes){
        x += n->point.X();
        y += n->point.Y();
    }

    return gp_Pnt(x/this->nodes.size(), y/this->nodes.size(), 0);
}
std::vector<double> TRI3::get_average_stress(const gp_Pnt& p1, const gp_Pnt& p2, const std::vector<double>& u) const{
    size_t N = this->nodes.size();

    std::vector<double> B(3*2*N, 0);

    std::vector<gp_Pnt> p;
    for(auto n:this->nodes){
        p.push_back(n->point);
    }

    gp_Mat deltaM(1, p[0].X(), p[0].Y(), 1, p[1].X(), p[1].Y(), 1, p[2].X(), p[2].Y());

    double Delta = 0.5*deltaM.Determinant();

    std::vector<double> a, b, c;
    for(size_t i = 0; i < N; ++i){
        size_t j = (i + 1) % 3;
        size_t k = (i + 2) % 3;

        a.push_back(p[j].X()*p[k].Y() - p[k].X()*p[j].Y());
        b.push_back(p[j].Y() - p[k].Y());
        c.push_back(p[k].X() - p[j].X());
    }

    double mult = (p2.X()-p1.X())*(p2.Y()-p1.Y())/p1.Distance(p2);

    for(size_t i = 0; i < N; ++i){
        B[i*2 + 0*2*N] = b[i]/(2*Delta)*mult;
        B[i*2 + 1*2*N] = 0;
        B[i*2 + 2*2*N] = c[i]/(2*Delta)*mult;
        B[i*2 + 0*2*N + 1] = 0;
        B[i*2 + 1*2*N + 1] = c[i]/(2*Delta)*mult;
        B[i*2 + 2*2*N + 1] = b[i]/(2*Delta)*mult;
    }

    std::vector<double> DB(3*2*N, 0);
    auto D = this->mat->stiffness_2D();

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 2*N, 3, 1, D.data(), 3, B.data(), 2*N, 0, DB.data(), 2*N);

    return DB;

}

std::vector<double> TRI3::get_DB(gp_Pnt point) const{
    (void)point;
    size_t N = this->nodes.size();

    std::vector<double> B(3*2*N, 0);

    std::vector<gp_Pnt> p;
    for(auto n:this->nodes){
        p.push_back(n->point);
    }

    gp_Mat deltaM(1, p[0].X(), p[0].Y(), 1, p[1].X(), p[1].Y(), 1, p[2].X(), p[2].Y());

    double Delta = 0.5*deltaM.Determinant();

    std::vector<double> a, b, c;
    for(size_t i = 0; i < N; ++i){
        size_t j = (i + 1) % 3;
        size_t k = (i + 2) % 3;

        a.push_back(p[j].X()*p[k].Y() - p[k].X()*p[j].Y());
        b.push_back(p[j].Y() - p[k].Y());
        c.push_back(p[k].X() - p[j].X());
    }

    for(size_t i = 0; i < N; ++i){
        B[i*2 + 0*2*N] = b[i]/(2*Delta);
        B[i*2 + 1*2*N] = 0;
        B[i*2 + 2*2*N] = c[i]/(2*Delta);
        B[i*2 + 0*2*N + 1] = 0;
        B[i*2 + 1*2*N + 1] = c[i]/(2*Delta);
        B[i*2 + 2*2*N + 1] = b[i]/(2*Delta);
    }

    std::vector<double> DB(3*2*N, 0);
    auto D = this->mat->stiffness_2D();

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 2*N, 3, 1, D.data(), 3, B.data(), 2*N, 0, DB.data(), 2*N);

    return DB;
}

}