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

#include "element/GT9.hpp"
#include "cblas.h"
#include "logger.hpp"
#include "project_data.hpp"

namespace element{

GT9::GT9(std::vector<long> u_pos, ElementShape s, ProjectData* data):
    MeshElement(s.nodes, std::move(u_pos)), mat(data->material.get()), t(data->thickness){}

std::vector<double> GT9::get_k() const{
    size_t N = this->nodes.size();

    std::vector<double> B(3*3*N, 0);

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
        size_t j = (i + 1) % 3;
        size_t k = (i + 2) % 3;

        double Ai = b[i]*b[k];
        double Bi = b[i]*b[j];
        double Ci = c[i]*c[k];
        double Di = c[i]*c[j];
        double Ei = c[i]*b[k] + b[i]*c[k];
        double Fi = c[i]*b[j] + b[i]*c[j];
        
        B[i*3 + 0*3*N] = (1/(4*Delta))*2*b[i];
        B[i*3 + 1*3*N] = 0;
        B[i*3 + 2*3*N] = (1/(4*Delta))*2*c[i];
        B[i*3 + 0*3*N + 1] = 0;
        B[i*3 + 1*3*N + 1] = (1/(4*Delta))*2*c[i];
        B[i*3 + 2*3*N + 1] = (1/(4*Delta))*2*b[i];
        B[i*3 + 0*3*N + 2] = (1/(4*Delta))*(Ai - Bi);
        B[i*3 + 1*3*N + 2] = (1/(4*Delta))*(Ci - Di);
        B[i*3 + 2*3*N + 2] = (1/(4*Delta))*(Ei - Fi);
    }

    std::vector<double> DB(3*3*N, 0);
    auto D = this->mat->stiffness_2D();

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3*N, 3, 1, D.data(), 3, B.data(), 3*N, 0, DB.data(), 3*N);

    std::vector<double> K(3*N*3*N, 0);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 3*N, 3*N, 3, 1, B.data(), 3*N, DB.data(), 3*N, 0, K.data(), 3*N);

    cblas_dscal(K.size(), this->t*Delta, K.data(), 1);
    for(size_t i = 2; i < 3*N; i += 3){
        for(size_t j = 0; j < 3*N; ++j){
            K[i + 3*N*j] *= 1.0/3;
        }
    }
    for(size_t i = 2; i < 3*N; i += 3){
        for(size_t j = 0; j < 3*N; ++j){
            K[j + 3*N*i] *= 1.0/3;
        }
    }
    K[2 + 3*N*2] *= 3.0/2;
    K[5 + 3*N*2] *= 3.0/4;
    K[8 + 3*N*2] *= 3.0/4;

    K[2 + 3*N*5] *= 3.0/4;
    K[5 + 3*N*5] *= 3.0/2;
    K[8 + 3*N*5] *= 3.0/4;

    K[2 + 3*N*8] *= 3.0/4;
    K[5 + 3*N*8] *= 3.0/4;
    K[8 + 3*N*8] *= 3.0/2;

    return K;
}

MeshNode* GT9::get_stresses(size_t node, const std::vector<double>& u) const{
    logger::log_assert(node >= 0 && node <= 3, logger::ERROR, "wrong value for BeamLinear2D node, must be either 0 or 1.");

    size_t N = this->nodes.size();

    double x = this->get_node(node)->point.X();
    double y = this->get_node(node)->point.Y();

    std::vector<double> B(3*3*N, 0);

    std::vector<gp_Pnt> p;
    for(auto n:this->nodes){
        p.push_back(n->point);
    }

    double Delta = (p[1].X()*p[2].Y() + p[0].X()*p[1].Y() + p[2].X()*p[0].Y())
                   - (p[1].X()*p[0].Y() + p[2].X()*p[1].Y() + p[0].X()*p[2].Y());

    std::vector<double> a, b, c;
    for(size_t i = 0; i < N; ++i){
        size_t j = (i + 1) % 3;
        size_t k = (i + 2) % 3;

        a.push_back(p[j].X()*p[k].Y() - p[k].X()*p[j].Y());
        b.push_back(p[j].Y() - p[k].Y());
        c.push_back(p[k].X() - p[j].X());
    }

    for(size_t i = 0; i < N; ++i){
        size_t j = (i + 1) % 3;
        size_t k = (i + 2) % 3;

        double Lj = a[j] + b[j]*x + c[j]*y;
        double Lk = a[k] + b[k]*x + c[k]*y;

        double Ai = b[i]*b[k];
        double Bi = b[i]*b[j];
        double Ci = c[i]*c[k];
        double Di = c[i]*c[j];
        double Ei = c[i]*b[k]*Lj + b[i]*c[k]*Lk;
        double Fi = c[i]*b[j]*Lj + b[i]*c[j]*Lk;
        
        B[i*3*3] = (1/(4*Delta))*2*b[i];
        B[i*3*3 + 1] = 0;
        B[i*3*3 + 2] = (1/(4*Delta))*(Ai - Bi);
        B[i*3*3 + 3*N] = 0;
        B[i*3*3 + 3*N + 1] = (1/(4*Delta))*2*c[i];
        B[i*3*3 + 3*N + 2] = (1/(4*Delta))*(Ci - Di);
        B[i*3*3 + 3*N] = (1/(4*Delta))*2*c[i];
        B[i*3*3 + 3*N + 1] = (1/(4*Delta))*2*b[i];
        B[i*3*3 + 3*N + 2] = (1/(4*Delta))*(Ei - Fi);
    }

    std::vector<double> DB(3*3*N, 0);
    auto D = this->mat->stiffness_2D();

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3*N, 3, 1, D.data(), 3, B.data(), 3*N, 0, DB.data(), 3*N);

    std::vector<double> K(3*N*3*N, 0);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 3*N, 3*N, 3, 1, B.data(), 3*N, DB.data(), 3*N, 0, K.data(), 3*N);

    MeshNode2D* n = static_cast<MeshNode2D*>(this->nodes[node]);
    for(size_t i = 0; i < 3; ++i){
        n->results[i] = 0;
        for(size_t j = 0; j < 3*N; ++j){
            if(this->u_pos[j] > -1){
                n->results[i] += DB[3*N*i + j]*u[this->u_pos[j]];
            }
        }
        n->results[i] = std::abs(n->results[i]);
    }

    return this->get_node(node);
}

MeshNode* GT9::get_internal_loads(size_t node, const std::vector<double>& u) const{
    logger::log_assert(node >= 0 && node <= 3, logger::ERROR, "wrong value for BeamLinear2D node, must be either 0 or 1.");

    std::vector<double> k = this->get_k();

    MeshNode2D* n = static_cast<MeshNode2D*>(this->nodes[node]);
    for(int i = 0; i < 3; ++i){
        n->results[i] = 0;
        for(int j = 0; j < 6; ++j){
            if(this->u_pos[j] > -1){
                n->results[i] += k[utils::to_triangular(node*3+i, j)]*u[this->u_pos[j]];
            }
        }
        n->results[i] = std::abs(n->results[i]);
    }

    return this->get_node(node);
}

}
