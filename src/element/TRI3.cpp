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

namespace element{

TRI3::TRI3(ElementShape s, ProjectData* data):
    MeshElement(s.nodes), mat(data->material.get()), t(data->thickness){}

std::vector<float> TRI3::get_k() const{
    size_t N = this->nodes.size();

    std::vector<float> B(3*2*N, 0);

    std::vector<gp_Pnt> p;
    for(auto n:this->nodes){
        p.push_back(n->point);
    }

    gp_Mat deltaM(1, p[0].X(), p[0].Y(), 1, p[1].X(), p[1].Y(), 1, p[2].X(), p[2].Y());

    float Delta = 0.5*std::abs(deltaM.Determinant());

    std::vector<float> a, b, c;
    for(size_t i = 0; i < N; ++i){
        size_t j = (i + 1) % 3;
        size_t k = (i + 2) % 3;

        a.push_back(p[j].X()*p[k].Y() - p[k].X()*p[j].Y());
        b.push_back(p[j].Y() - p[k].Y());
        c.push_back(p[k].X() - p[j].X());
    }

    for(size_t i = 0; i < N; ++i){
        B[i*2 + 0*3*N] = b[i];
        B[i*2 + 1*3*N] = 0;
        B[i*2 + 2*3*N] = c[i];
        B[i*2 + 0*3*N + 1] = 0;
        B[i*2 + 1*3*N + 1] = c[i];
        B[i*2 + 2*3*N + 1] = b[i];
    }

    std::vector<float> DB(3*2*N, 0);
    auto D = this->mat->stiffness_2D();

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 2*N, 3, 1, D.data(), 3, B.data(), 2*N, 0, DB.data(), 2*N);

    std::vector<float> K(2*N*2*N, 0);
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 2*N, 2*N, 3, 1, B.data(), 2*N, DB.data(), 2*N, 0, K.data(), 2*N);

    cblas_sscal(K.size(), this->t/(4*Delta), K.data(), 1);

    return K;
}

MeshNode* TRI3::get_stresses(size_t node, const std::vector<float>& u, double density) const{
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
        n->results[i] = density*1e6*std::abs(n->results[i]);
    }

    return this->get_node(node);
}

double TRI3::get_stress_at(gp_Pnt point, const std::vector<float>& u) const{
    size_t N = this->nodes.size();

    auto DB = this->get_DB(point);

    std::vector<double> results(3, 0);
    for(size_t i = 0; i < 3; ++i){
        results[i] = 0;
        for(size_t l = 0; l < 3; ++l){
            for(size_t j = 0; j < 2; ++j){
                if(this->nodes[l]->u_pos[j] > -1){
                    results[i] += DB[2*N*i + l*2 + j]*u[this->nodes[l]->u_pos[j]];
                }
            }
        }
        results[i] = std::abs(results[i]);
    }

    return 1e6*std::sqrt(std::pow(results[0], 2) - results[0]*results[1] + std::pow(results[1], 2) + 3*std::pow(results[2], 2));
}

MeshNode* TRI3::get_internal_loads(size_t node, const std::vector<float>& u) const{
    logger::log_assert(node >= 0 && node <= 3, logger::ERROR, "wrong value for BeamLinear2D node, must be either 0 or 1.");

    std::vector<float> k = this->get_k();

    MeshNode2D* n = static_cast<MeshNode2D*>(this->nodes[node]);
    for(int i = 0; i < 3; ++i){
        n->results[i] = 0;
        for(size_t l = 0; l < 3; ++l){
            for(int j = 0; j < 2; ++j){
                if(this->nodes[l]->u_pos[j] > -1){
                    n->results[i] += k[node*2*3+l*2+j]*u[this->nodes[l]->u_pos[j]];
                }
            }
        }
        n->results[i] = std::abs(n->results[i]);
    }

    return this->get_node(node);
}

double TRI3::get_compliance(const std::vector<float>& u, const std::vector<float>& l) const{
    auto k = this->get_k();
    std::vector<float> u_vec(6, 0);
    for(size_t k = 0; k < 3; ++k){
        for(int j = 0; j < 2; ++j){
            u_vec[k*2+j] = u[this->nodes[k]->u_pos[j]];
        }
    }

    std::vector<float> f_vec(6, 0);
    if(l.size() > 0){
        std::vector<float> l_vec(6, 0);
        for(size_t k = 0; k < 3; ++k){
            for(int j = 0; j < 2; ++j){
                l_vec[k*2+j] = l[this->nodes[k]->u_pos[j]];
            }
        }
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 6, 1, 6, 1, k.data(), 6, l_vec.data(), 6, 0, f_vec.data(), 6);
    } else {
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 6, 1, 6, 1, k.data(), 6, u_vec.data(), 6, 0, f_vec.data(), 6);
    }

    return 1e-3*cblas_sdot(6, u_vec.data(), 1, f_vec.data(), 1);
}

double TRI3::get_volume() const{
    gp_Mat deltaM(1, this->nodes[0]->point.X(), this->nodes[0]->point.Y(),
                  1, this->nodes[1]->point.X(), this->nodes[1]->point.Y(),
                  1, this->nodes[2]->point.X(), this->nodes[2]->point.Y());

    return 0.5*std::abs(deltaM.Determinant())*1e-6;
}

void TRI3::get_virtual_load(double P, gp_Pnt point, std::vector<float>& u, std::vector<float>& l) const{
    std::vector<float> DB = this->get_DB(point);
    double stress = this->get_stress_at(point, u);
    std::vector<float> V{1, -0.5, 0,
                         -0.5, 1, 0,
                         0,   0, 1.5};

    std::vector<float> u_vec(6, 0);
    for(size_t k = 0; k < 3; ++k){
        for(int j = 0; j < 2; ++j){
            u_vec[k*2+j] = u[this->nodes[k]->u_pos[j]];
        }
    }

    std::vector<float> f_vec(6, 0);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 1, 6, 1, DB.data(), 1, u_vec.data(), 6, 0, f_vec.data(), 6);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 1, 3, 1, V.data(), 1, f_vec.data(), 6, 0, f_vec.data(), 6);
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 3, 1, 6, 1, DB.data(), 1, u_vec.data(), 6, 0, f_vec.data(), 6);
    cblas_sscal(6, P*std::pow(stress, P-2), f_vec.data(), 1);

    for(size_t k = 0; k < 3; ++k){
        for(int j = 0; j < 2; ++j){
            u[this->nodes[k]->u_pos[j]] = f_vec[k*2+j];
        }
    }
}

TopoDS_Shape TRI3::get_shape() const{
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

std::vector<float> TRI3::get_DB(gp_Pnt point) const{
    (void)point;
    size_t N = this->nodes.size();

    std::vector<float> B(3*2*N, 0);

    std::vector<gp_Pnt> p;
    for(auto n:this->nodes){
        p.push_back(n->point);
    }

    gp_Mat deltaM(1, p[0].X(), p[0].Y(), 1, p[1].X(), p[1].Y(), 1, p[2].X(), p[2].Y());

    float Delta = 0.5*std::abs(deltaM.Determinant());

    std::vector<float> a, b, c;
    for(size_t i = 0; i < N; ++i){
        size_t j = (i + 1) % 3;
        size_t k = (i + 2) % 3;

        a.push_back(p[j].X()*p[k].Y() - p[k].X()*p[j].Y());
        b.push_back(p[j].Y() - p[k].Y());
        c.push_back(p[k].X() - p[j].X());
    }

    for(size_t i = 0; i < N; ++i){
        B[i*2 + 0*3*N] = b[i];
        B[i*2 + 1*3*N] = 0;
        B[i*2 + 2*3*N] = c[i];
        B[i*2 + 0*3*N + 1] = 0;
        B[i*2 + 1*3*N + 1] = c[i];
        B[i*2 + 2*3*N + 1] = b[i];
    }

    std::vector<float> DB(3*2*N, 0);
    auto D = this->mat->stiffness_2D();

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 2*N, 3, 1, D.data(), 3, B.data(), 2*N, 0, DB.data(), 2*N);
    cblas_sscal(DB.size(), 1/(2*Delta), DB.data(), 1);

    return DB;
}

}
