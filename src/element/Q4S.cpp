/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
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

#include "element/Q4S.hpp"
#include "cblas.h"
#include "logger.hpp"
#include <vector>
#include <lapacke.h>

namespace element{

Q4S::Q4S(ElementShape s):
    MeshElementCommon2DQuad<Q4S>(s.nodes, s.nodes_sorted)
    {}

std::vector<double> Q4S::get_k(const std::vector<double>& D, const double t) const{
    double a, b;
    this->get_dimensions(a, b);

    const auto k0 = this->get_k_base(D, t, a, b);

    // Node ordering is important in this case, so we need to use the same
    // ordering that Gmsh uses, as it's the same one that is used for the
    // calculations too.
    //
    // Otherwise, the global matrix will have a ton of zeros.
    constexpr size_t N = Q4S::NODES_PER_ELEM;
    std::vector<size_t> pos;
    for(size_t i = 0; i < N; ++i){
        for(size_t j = 0; j < N; ++j){
            if(this->nodes_sorted[j]->id == this->nodes[i]->id){
                pos.push_back(j*2);
                pos.push_back(j*2+1);
                break;
            }
        }
    }

    // Afterwards, we need to get the matrix into a sorted-by-id composition,
    // which can easily be node in a way analogous to how the global matrix
    // K is formed.
    std::vector<double> k(Q4S::K_DIM*Q4S::K_DIM);
    for(size_t i = 0; i < Q4S::K_DIM; ++i){
        for(size_t j = 0; j < Q4S::K_DIM; ++j){
            const size_t p = i*Q4S::K_DIM+j;
            const size_t ps = pos[i]*Q4S::K_DIM+pos[j];
            k[ps] = k0[p];
        }
    }

    return k;
}

std::vector<double> Q4S::get_DB(const std::vector<double>& D, const gp_Pnt& point) const{
    double a, b;
    this->get_dimensions(a, b);

    const gp_Pnt p = this->normalize(point, a, b);

    const double xi = p.X();
    const double eta = p.Y();

    const auto DB0 = this->get_DB_base(D, xi, eta, a, b);

    // Node ordering is important in this case, so we need to use the same
    // ordering that Gmsh uses, as it's the same one that is used for the
    // calculations too.
    //
    // Otherwise, the global matrix will have a ton of zeros.
    constexpr size_t N = Q4S::NODES_PER_ELEM;
    std::vector<size_t> pos;
    for(size_t i = 0; i < N; ++i){
        for(size_t j = 0; j < N; ++j){
            if(this->nodes_sorted[j]->id == this->nodes[i]->id){
                pos.push_back(j*2);
                pos.push_back(j*2+1);
                break;
            }
        }
    }

    // Afterwards, we need to get the matrix into a sorted-by-id composition,
    // which can easily be node in a way analogous to how the global matrix
    // K is formed.
    std::vector<double> DB(DB0.size());
    for(size_t i = 0; i < Q4S::K_DIM; ++i){
        for(size_t j = 0; j < 3; ++j){
            const size_t p = j*Q4S::K_DIM+i;
            const size_t ps = pos[j]*Q4S::K_DIM+pos[i];
            DB[ps] = DB0[p];
        }
    }

    return DB;
}

std::vector<double> Q4S::get_Nf(const double t, const std::vector<gp_Pnt>& points) const{
    double a, b;
    this->get_dimensions(a, b);

    const gp_Pnt p0 = this->normalize(points[0], a, b);
    const gp_Pnt p1 = this->normalize(points[1], a, b);

    const std::array<double, 2> x{p0.X(), p1.X()};
    const std::array<double, 2> y{p0.Y(), p1.Y()};

    const auto Nf0 = this->get_Nf_base(t, x, y, a, b);

    // Node ordering is important in this case, so we need to use the same
    // ordering that Gmsh uses, as it's the same one that is used for the
    // calculations too.
    //
    // Otherwise, the global matrix will have a ton of zeros.
    constexpr size_t N = Q4S::NODES_PER_ELEM;
    std::vector<size_t> pos;
    for(size_t i = 0; i < N; ++i){
        for(size_t j = 0; j < N; ++j){
            if(this->nodes_sorted[j]->id == this->nodes[i]->id){
                pos.push_back(j*2);
                pos.push_back(j*2+1);
                break;
            }
        }
    }

    // Afterwards, we need to get the matrix into a sorted-by-id composition,
    // which can easily be node in a way analogous to how the global matrix
    // K is formed.
    std::vector<double> Nf(Nf0.size());
    for(size_t i = 0; i < Q4S::K_DIM; ++i){
        for(size_t j = 0; j < Q4S::NODE_DOF; ++j){
            const size_t p = i*Q4S::NODE_DOF+j;
            const size_t ps = pos[i]*Q4S::NODE_DOF+pos[j];
            Nf[ps] = Nf0[p];
        }
    }

    return Nf;
}

std::vector<double> Q4S::get_k_base(const std::vector<double>& D, const double t, const double a, const double b) const{

    std::vector<double> k{
    t*(3*D[0]*b*b + D[0] + 3*D[8]*a*a + D[8] + 3*a*b*(D[2] + D[6]))/(12*a*a*b*b)
    ,
    t*(3*D[2]*b*b + D[2] + 3*D[7]*a*a + D[7] + 3*a*b*(D[1] + D[8]))/(12*a*a*b*b)
    ,
    t*(-3*D[0]*b*b - D[0] + 3*D[8]*a*a - D[8] + 3*a*b*(D[2] - D[6]))/(12*a*a*b*b)
    ,
    t*(-3*D[2]*b*b - D[2] + 3*D[7]*a*a - D[7] + 3*a*b*(D[1] - D[8]))/(12*a*a*b*b)
    ,
    t*(-3*D[0]*b*b + D[0] - 3*D[8]*a*a + D[8] - 3*a*b*(D[2] + D[6]))/(12*a*a*b*b)
    ,
    t*(-3*D[2]*b*b + D[2] - 3*D[7]*a*a + D[7] - 3*a*b*(D[1] + D[8]))/(12*a*a*b*b)
    ,
    t*(3*D[0]*b*b - D[0] - 3*D[8]*a*a - D[8] + 3*a*b*(-D[2] + D[6]))/(12*a*a*b*b)
    ,
    t*(3*D[2]*b*b - D[2] - 3*D[7]*a*a - D[7] + 3*a*b*(-D[1] + D[8]))/(12*a*a*b*b)
    ,
    t*(3*D[5]*a*a + D[5] + 3*D[6]*b*b + D[6] + 3*a*b*(D[3] + D[8]))/(12*a*a*b*b)
    ,
    t*(3*D[4]*a*a + D[4] + 3*D[8]*b*b + D[8] + 3*a*b*(D[5] + D[7]))/(12*a*a*b*b)
    ,
    t*(3*D[5]*a*a - D[5] - 3*D[6]*b*b - D[6] + 3*a*b*(-D[3] + D[8]))/(12*a*a*b*b)
    ,
    t*(3*D[4]*a*a - D[4] - 3*D[8]*b*b - D[8] + 3*a*b*(-D[5] + D[7]))/(12*a*a*b*b)
    ,
    t*(-3*D[5]*a*a + D[5] - 3*D[6]*b*b + D[6] - 3*a*b*(D[3] + D[8]))/(12*a*a*b*b)
    ,
    t*(-3*D[4]*a*a + D[4] - 3*D[8]*b*b + D[8] - 3*a*b*(D[5] + D[7]))/(12*a*a*b*b)
    ,
    t*(-3*D[5]*a*a - D[5] + 3*D[6]*b*b - D[6] + 3*a*b*(D[3] - D[8]))/(12*a*a*b*b)
    ,
    t*(-3*D[4]*a*a - D[4] + 3*D[8]*b*b - D[8] + 3*a*b*(D[5] - D[7]))/(12*a*a*b*b)
    ,
    t*(-3*D[0]*b*b - D[0] + 3*D[8]*a*a - D[8] + 3*a*b*(-D[2] + D[6]))/(12*a*a*b*b)
    ,
    t*(-3*D[2]*b*b - D[2] + 3*D[7]*a*a - D[7] + 3*a*b*(-D[1] + D[8]))/(12*a*a*b*b)
    ,
    t*(3*D[0]*b*b + D[0] + 3*D[8]*a*a + D[8] - 3*a*b*(D[2] + D[6]))/(12*a*a*b*b)
    ,
    t*(3*D[2]*b*b + D[2] + 3*D[7]*a*a + D[7] - 3*a*b*(D[1] + D[8]))/(12*a*a*b*b)
    ,
    t*(3*D[0]*b*b - D[0] - 3*D[8]*a*a - D[8] + 3*a*b*(D[2] - D[6]))/(12*a*a*b*b)
    ,
    t*(3*D[2]*b*b - D[2] - 3*D[7]*a*a - D[7] + 3*a*b*(D[1] - D[8]))/(12*a*a*b*b)
    ,
    t*(-3*D[0]*b*b + D[0] - 3*D[8]*a*a + D[8] + 3*a*b*(D[2] + D[6]))/(12*a*a*b*b)
    ,
    t*(-3*D[2]*b*b + D[2] - 3*D[7]*a*a + D[7] + 3*a*b*(D[1] + D[8]))/(12*a*a*b*b)
    ,
    t*(3*D[5]*a*a - D[5] - 3*D[6]*b*b - D[6] + 3*a*b*(D[3] - D[8]))/(12*a*a*b*b)
    ,
    t*(3*D[4]*a*a - D[4] - 3*D[8]*b*b - D[8] + 3*a*b*(D[5] - D[7]))/(12*a*a*b*b)
    ,
    t*(3*D[5]*a*a + D[5] + 3*D[6]*b*b + D[6] - 3*a*b*(D[3] + D[8]))/(12*a*a*b*b)
    ,
    t*(3*D[4]*a*a + D[4] + 3*D[8]*b*b + D[8] - 3*a*b*(D[5] + D[7]))/(12*a*a*b*b)
    ,
    t*(-3*D[5]*a*a - D[5] + 3*D[6]*b*b - D[6] + 3*a*b*(-D[3] + D[8]))/(12*a*a*b*b)
    ,
    t*(-3*D[4]*a*a - D[4] + 3*D[8]*b*b - D[8] + 3*a*b*(-D[5] + D[7]))/(12*a*a*b*b)
    ,
    t*(-3*D[5]*a*a + D[5] - 3*D[6]*b*b + D[6] + 3*a*b*(D[3] + D[8]))/(12*a*a*b*b)
    ,
    t*(-3*D[4]*a*a + D[4] - 3*D[8]*b*b + D[8] + 3*a*b*(D[5] + D[7]))/(12*a*a*b*b)
    ,
    t*(-3*D[0]*b*b + D[0] - 3*D[8]*a*a + D[8] - 3*a*b*(D[2] + D[6]))/(12*a*a*b*b)
    ,
    t*(-3*D[2]*b*b + D[2] - 3*D[7]*a*a + D[7] - 3*a*b*(D[1] + D[8]))/(12*a*a*b*b)
    ,
    t*(3*D[0]*b*b - D[0] - 3*D[8]*a*a - D[8] + 3*a*b*(-D[2] + D[6]))/(12*a*a*b*b)
    ,
    t*(3*D[2]*b*b - D[2] - 3*D[7]*a*a - D[7] + 3*a*b*(-D[1] + D[8]))/(12*a*a*b*b)
    ,
    t*(3*D[0]*b*b + D[0] + 3*D[8]*a*a + D[8] + 3*a*b*(D[2] + D[6]))/(12*a*a*b*b)
    ,
    t*(3*D[2]*b*b + D[2] + 3*D[7]*a*a + D[7] + 3*a*b*(D[1] + D[8]))/(12*a*a*b*b)
    ,
    t*(-3*D[0]*b*b - D[0] + 3*D[8]*a*a - D[8] + 3*a*b*(D[2] - D[6]))/(12*a*a*b*b)
    ,
    t*(-3*D[2]*b*b - D[2] + 3*D[7]*a*a - D[7] + 3*a*b*(D[1] - D[8]))/(12*a*a*b*b)
    ,
    t*(-3*D[5]*a*a + D[5] - 3*D[6]*b*b + D[6] - 3*a*b*(D[3] + D[8]))/(12*a*a*b*b)
    ,
    t*(-3*D[4]*a*a + D[4] - 3*D[8]*b*b + D[8] - 3*a*b*(D[5] + D[7]))/(12*a*a*b*b)
    ,
    t*(-3*D[5]*a*a - D[5] + 3*D[6]*b*b - D[6] + 3*a*b*(D[3] - D[8]))/(12*a*a*b*b)
    ,
    t*(-3*D[4]*a*a - D[4] + 3*D[8]*b*b - D[8] + 3*a*b*(D[5] - D[7]))/(12*a*a*b*b)
    ,
    t*(3*D[5]*a*a + D[5] + 3*D[6]*b*b + D[6] + 3*a*b*(D[3] + D[8]))/(12*a*a*b*b)
    ,
    t*(3*D[4]*a*a + D[4] + 3*D[8]*b*b + D[8] + 3*a*b*(D[5] + D[7]))/(12*a*a*b*b)
    ,
    t*(3*D[5]*a*a - D[5] - 3*D[6]*b*b - D[6] + 3*a*b*(-D[3] + D[8]))/(12*a*a*b*b)
    ,
    t*(3*D[4]*a*a - D[4] - 3*D[8]*b*b - D[8] + 3*a*b*(-D[5] + D[7]))/(12*a*a*b*b)
    ,
    t*(3*D[0]*b*b - D[0] - 3*D[8]*a*a - D[8] + 3*a*b*(D[2] - D[6]))/(12*a*a*b*b)
    ,
    t*(3*D[2]*b*b - D[2] - 3*D[7]*a*a - D[7] + 3*a*b*(D[1] - D[8]))/(12*a*a*b*b)
    ,
    t*(-3*D[0]*b*b + D[0] - 3*D[8]*a*a + D[8] + 3*a*b*(D[2] + D[6]))/(12*a*a*b*b)
    ,
    t*(-3*D[2]*b*b + D[2] - 3*D[7]*a*a + D[7] + 3*a*b*(D[1] + D[8]))/(12*a*a*b*b)
    ,
    t*(-3*D[0]*b*b - D[0] + 3*D[8]*a*a - D[8] + 3*a*b*(-D[2] + D[6]))/(12*a*a*b*b)
    ,
    t*(-3*D[2]*b*b - D[2] + 3*D[7]*a*a - D[7] + 3*a*b*(-D[1] + D[8]))/(12*a*a*b*b)
    ,
    t*(3*D[0]*b*b + D[0] + 3*D[8]*a*a + D[8] - 3*a*b*(D[2] + D[6]))/(12*a*a*b*b)
    ,
    t*(3*D[2]*b*b + D[2] + 3*D[7]*a*a + D[7] - 3*a*b*(D[1] + D[8]))/(12*a*a*b*b)
    ,
    t*(-3*D[5]*a*a - D[5] + 3*D[6]*b*b - D[6] + 3*a*b*(-D[3] + D[8]))/(12*a*a*b*b)
    ,
    t*(-3*D[4]*a*a - D[4] + 3*D[8]*b*b - D[8] + 3*a*b*(-D[5] + D[7]))/(12*a*a*b*b)
    ,
    t*(-3*D[5]*a*a + D[5] - 3*D[6]*b*b + D[6] + 3*a*b*(D[3] + D[8]))/(12*a*a*b*b)
    ,
    t*(-3*D[4]*a*a + D[4] - 3*D[8]*b*b + D[8] + 3*a*b*(D[5] + D[7]))/(12*a*a*b*b)
    ,
    t*(3*D[5]*a*a - D[5] - 3*D[6]*b*b - D[6] + 3*a*b*(D[3] - D[8]))/(12*a*a*b*b)
    ,
    t*(3*D[4]*a*a - D[4] - 3*D[8]*b*b - D[8] + 3*a*b*(D[5] - D[7]))/(12*a*a*b*b)
    ,
    t*(3*D[5]*a*a + D[5] + 3*D[6]*b*b + D[6] - 3*a*b*(D[3] + D[8]))/(12*a*a*b*b)
    ,
    t*(3*D[4]*a*a + D[4] + 3*D[8]*b*b + D[8] - 3*a*b*(D[5] + D[7]))/(12*a*a*b*b)
    };

    return k;
}

std::vector<double> Q4S::get_DB_base(const std::vector<double>& D, const double xi, const double eta, const double a, const double b) const{

    std::vector<double> DB{
    (-D[0]*b + D[0]*eta - D[2]*a + D[2]*xi)/(4*a*b)
    ,
    (-D[1]*a + D[1]*xi - D[2]*b + D[2]*eta)/(4*a*b)
    ,
    (D[0]*b - D[0]*eta - D[2]*a - D[2]*xi)/(4*a*b)
    ,
    (-D[1]*a - D[1]*xi + D[2]*b - D[2]*eta)/(4*a*b)
    ,
    (D[0]*b + D[0]*eta + D[2]*a + D[2]*xi)/(4*a*b)
    ,
    (D[1]*a + D[1]*xi + D[2]*b + D[2]*eta)/(4*a*b)
    ,
    (-D[0]*b - D[0]*eta + D[2]*a - D[2]*xi)/(4*a*b)
    ,
    (D[1]*a - D[1]*xi - D[2]*b - D[2]*eta)/(4*a*b)
    ,
    (-D[3]*b + D[3]*eta - D[5]*a + D[5]*xi)/(4*a*b)
    ,
    (-D[4]*a + D[4]*xi - D[5]*b + D[5]*eta)/(4*a*b)
    ,
    (D[3]*b - D[3]*eta - D[5]*a - D[5]*xi)/(4*a*b)
    ,
    (-D[4]*a - D[4]*xi + D[5]*b - D[5]*eta)/(4*a*b)
    ,
    (D[3]*b + D[3]*eta + D[5]*a + D[5]*xi)/(4*a*b)
    ,
    (D[4]*a + D[4]*xi + D[5]*b + D[5]*eta)/(4*a*b)
    ,
    (-D[3]*b - D[3]*eta + D[5]*a - D[5]*xi)/(4*a*b)
    ,
    (D[4]*a - D[4]*xi - D[5]*b - D[5]*eta)/(4*a*b)
    ,
    (-D[6]*b + D[6]*eta - D[8]*a + D[8]*xi)/(4*a*b)
    ,
    (-D[7]*a + D[7]*xi - D[8]*b + D[8]*eta)/(4*a*b)
    ,
    (D[6]*b - D[6]*eta - D[8]*a - D[8]*xi)/(4*a*b)
    ,
    (-D[7]*a - D[7]*xi + D[8]*b - D[8]*eta)/(4*a*b)
    ,
    (D[6]*b + D[6]*eta + D[8]*a + D[8]*xi)/(4*a*b)
    ,
    (D[7]*a + D[7]*xi + D[8]*b + D[8]*eta)/(4*a*b)
    ,
    (-D[6]*b - D[6]*eta + D[8]*a - D[8]*xi)/(4*a*b)
    ,
    (D[7]*a - D[7]*xi - D[8]*b - D[8]*eta)/(4*a*b)
    };

    return DB;

}
std::vector<double> Q4S::get_Nf_base(const double t, const std::array<double, 2> x, const std::array<double, 2> y, const double a, const double b) const{

    std::vector<double> Nf{
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a*b - 3*a*(y[0] + y[1]) - 3*b*(x[0] + x[1]) + 2*x[0]*y[0] + x[0]*y[1] + x[1]*y[0] + 2*x[1]*y[1])/(24*a*b)
    ,
    0
    ,
    0
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a*b - 3*a*(y[0] + y[1]) - 3*b*(x[0] + x[1]) + 2*x[0]*y[0] + x[0]*y[1] + x[1]*y[0] + 2*x[1]*y[1])/(24*a*b)
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a*b - 3*a*(y[0] + y[1]) + 3*b*(x[0] + x[1]) - 2*x[0]*y[0] - x[0]*y[1] - x[1]*y[0] - 2*x[1]*y[1])/(24*a*b)
    ,
    0
    ,
    0
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a*b - 3*a*(y[0] + y[1]) + 3*b*(x[0] + x[1]) - 2*x[0]*y[0] - x[0]*y[1] - x[1]*y[0] - 2*x[1]*y[1])/(24*a*b)
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a*b + 3*a*(y[0] + y[1]) + 3*b*(x[0] + x[1]) + 2*x[0]*y[0] + x[0]*y[1] + x[1]*y[0] + 2*x[1]*y[1])/(24*a*b)
    ,
    0
    ,
    0
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a*b + 3*a*(y[0] + y[1]) + 3*b*(x[0] + x[1]) + 2*x[0]*y[0] + x[0]*y[1] + x[1]*y[0] + 2*x[1]*y[1])/(24*a*b)
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a*b + 3*a*(y[0] + y[1]) - 3*b*(x[0] + x[1]) - 2*x[0]*y[0] - x[0]*y[1] - x[1]*y[0] - 2*x[1]*y[1])/(24*a*b)
    ,
    0
    ,
    0
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a*b + 3*a*(y[0] + y[1]) - 3*b*(x[0] + x[1]) - 2*x[0]*y[0] - x[0]*y[1] - x[1]*y[0] - 2*x[1]*y[1])/(24*a*b)
    };

    return Nf;

}

}


