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
#include <vector>
#include <lapacke.h>

namespace element{

Q4S::Q4S(ElementShape s):
    MeshElementCommon2DQuad<Q4S>(s.nodes){

    constexpr size_t N = Q4S::NODES_PER_ELEM;

    double maxx = this->nodes[0]->point.X();
    double minx = this->nodes[0]->point.X();
    double maxy = this->nodes[0]->point.Y();
    double miny = this->nodes[0]->point.Y();
    auto c = this->get_centroid();
    for(size_t i = 1; i < N; ++i){
        const double x = this->nodes[i]->point.X();
        const double y = this->nodes[i]->point.Y();
        if(x > c.X() && y > c.Y()){
            maxx = x;
            maxy = y;
        } else if(x < c.X() && y < c.Y()){
            minx = x;
            miny = y;
        }
    }
    this->x0 = minx;
    this->y0 = miny;
    this->a = (maxx-minx)/2;
    this->b = (maxy-miny)/2;
    // std::array<double, N> xl{minx, maxx, maxx, minx};// x0, x0+2*a, x0+2*a, x0};
    // std::array<double, N> yl{miny, miny, maxy, maxy};// y0, y0, y0+2*b, y0+2*b};

    // for(size_t i = 0; i < N; ++i){
    //     for(size_t j = 0; j < N; ++j){
    //         if(this->nodes[j]->point.IsEqual({xl[i], yl[i], 0.0}, 1e-1)){
    //             std::cout << j << " ";
    //             break;
    //         }
    //     }
    // }
    // std::cout << std::endl;

    // Enforce node ordering by rotating array
    while(!this->nodes[0]->point.IsEqual({x0, y0, 0.0}, Precision::Confusion())){
        for(size_t i = 0; i < N-1; ++i){
            std::swap(this->nodes[i], this->nodes[i+1]);
        }
    }
            
    // for(size_t i = 0; i < N-1; ++i){
    //     if(!this->nodes[i]->point.IsEqual({xl[i], yl[i], 0.0}, eps)){
    //         for(size_t j = i+1; j < N; ++j){
    //             if(this->nodes[j]->point.IsEqual({xl[i], yl[i], 0.0}, eps)){
    //                 std::swap(this->nodes[i], this->nodes[j]);
    //                 break;
    //             }
    //         }
    //     }
    // }  
}

std::vector<double> Q4S::get_k(const std::vector<double>& D, const double t) const{

    std::vector<double> k{
    D[2]*t/4 + D[6]*t/4 + (D[0]*b*b*t + D[8]*a*a*t)/(3*a*b)
    ,
    D[1]*t/4 + D[8]*t/4 + (D[2]*b*b*t + D[7]*a*a*t)/(3*a*b)
    ,
    D[2]*t/4 - D[6]*t/4 + (-2*D[0]*b*b*t + D[8]*a*a*t)/(6*a*b)
    ,
    D[1]*t/4 - D[8]*t/4 + (-2*D[2]*b*b*t + D[7]*a*a*t)/(6*a*b)
    ,
    -D[2]*t/4 - D[6]*t/4 + (-D[0]*b*b*t - D[8]*a*a*t)/(6*a*b)
    ,
    -D[1]*t/4 - D[8]*t/4 + (-D[2]*b*b*t - D[7]*a*a*t)/(6*a*b)
    ,
    -D[2]*t/4 + D[6]*t/4 + (D[0]*b*b*t - 2*D[8]*a*a*t)/(6*a*b)
    ,
    -D[1]*t/4 + D[8]*t/4 + (D[2]*b*b*t - 2*D[7]*a*a*t)/(6*a*b)
    ,
    D[3]*t/4 + D[8]*t/4 + (D[5]*a*a*t + D[6]*b*b*t)/(3*a*b)
    ,
    D[5]*t/4 + D[7]*t/4 + (D[4]*a*a*t + D[8]*b*b*t)/(3*a*b)
    ,
    -D[3]*t/4 + D[8]*t/4 + (D[5]*a*a*t - 2*D[6]*b*b*t)/(6*a*b)
    ,
    -D[5]*t/4 + D[7]*t/4 + (D[4]*a*a*t - 2*D[8]*b*b*t)/(6*a*b)
    ,
    -D[3]*t/4 - D[8]*t/4 + (-D[5]*a*a*t - D[6]*b*b*t)/(6*a*b)
    ,
    -D[5]*t/4 - D[7]*t/4 + (-D[4]*a*a*t - D[8]*b*b*t)/(6*a*b)
    ,
    D[3]*t/4 - D[8]*t/4 + (-2*D[5]*a*a*t + D[6]*b*b*t)/(6*a*b)
    ,
    D[5]*t/4 - D[7]*t/4 + (-2*D[4]*a*a*t + D[8]*b*b*t)/(6*a*b)
    ,
    -D[2]*t/4 + D[6]*t/4 + (-2*D[0]*b*b*t + D[8]*a*a*t)/(6*a*b)
    ,
    -D[1]*t/4 + D[8]*t/4 + (-2*D[2]*b*b*t + D[7]*a*a*t)/(6*a*b)
    ,
    -D[2]*t/4 - D[6]*t/4 + (D[0]*b*b*t + D[8]*a*a*t)/(3*a*b)
    ,
    -D[1]*t/4 - D[8]*t/4 + (D[2]*b*b*t + D[7]*a*a*t)/(3*a*b)
    ,
    D[2]*t/4 - D[6]*t/4 + (D[0]*b*b*t - 2*D[8]*a*a*t)/(6*a*b)
    ,
    D[1]*t/4 - D[8]*t/4 + (D[2]*b*b*t - 2*D[7]*a*a*t)/(6*a*b)
    ,
    D[2]*t/4 + D[6]*t/4 + (-D[0]*b*b*t - D[8]*a*a*t)/(6*a*b)
    ,
    D[1]*t/4 + D[8]*t/4 + (-D[2]*b*b*t - D[7]*a*a*t)/(6*a*b)
    ,
    D[3]*t/4 - D[8]*t/4 + (D[5]*a*a*t - 2*D[6]*b*b*t)/(6*a*b)
    ,
    D[5]*t/4 - D[7]*t/4 + (D[4]*a*a*t - 2*D[8]*b*b*t)/(6*a*b)
    ,
    -D[3]*t/4 - D[8]*t/4 + (D[5]*a*a*t + D[6]*b*b*t)/(3*a*b)
    ,
    -D[5]*t/4 - D[7]*t/4 + (D[4]*a*a*t + D[8]*b*b*t)/(3*a*b)
    ,
    -D[3]*t/4 + D[8]*t/4 + (-2*D[5]*a*a*t + D[6]*b*b*t)/(6*a*b)
    ,
    -D[5]*t/4 + D[7]*t/4 + (-2*D[4]*a*a*t + D[8]*b*b*t)/(6*a*b)
    ,
    D[3]*t/4 + D[8]*t/4 + (-D[5]*a*a*t - D[6]*b*b*t)/(6*a*b)
    ,
    D[5]*t/4 + D[7]*t/4 + (-D[4]*a*a*t - D[8]*b*b*t)/(6*a*b)
    ,
    -D[2]*t/4 - D[6]*t/4 + (-D[0]*b*b*t - D[8]*a*a*t)/(6*a*b)
    ,
    -D[1]*t/4 - D[8]*t/4 + (-D[2]*b*b*t - D[7]*a*a*t)/(6*a*b)
    ,
    -D[2]*t/4 + D[6]*t/4 + (D[0]*b*b*t - 2*D[8]*a*a*t)/(6*a*b)
    ,
    -D[1]*t/4 + D[8]*t/4 + (D[2]*b*b*t - 2*D[7]*a*a*t)/(6*a*b)
    ,
    D[2]*t/4 + D[6]*t/4 + (D[0]*b*b*t + D[8]*a*a*t)/(3*a*b)
    ,
    D[1]*t/4 + D[8]*t/4 + (D[2]*b*b*t + D[7]*a*a*t)/(3*a*b)
    ,
    D[2]*t/4 - D[6]*t/4 + (-2*D[0]*b*b*t + D[8]*a*a*t)/(6*a*b)
    ,
    D[1]*t/4 - D[8]*t/4 + (-2*D[2]*b*b*t + D[7]*a*a*t)/(6*a*b)
    ,
    -D[3]*t/4 - D[8]*t/4 + (-D[5]*a*a*t - D[6]*b*b*t)/(6*a*b)
    ,
    -D[5]*t/4 - D[7]*t/4 + (-D[4]*a*a*t - D[8]*b*b*t)/(6*a*b)
    ,
    D[3]*t/4 - D[8]*t/4 + (-2*D[5]*a*a*t + D[6]*b*b*t)/(6*a*b)
    ,
    D[5]*t/4 - D[7]*t/4 + (-2*D[4]*a*a*t + D[8]*b*b*t)/(6*a*b)
    ,
    D[3]*t/4 + D[8]*t/4 + (D[5]*a*a*t + D[6]*b*b*t)/(3*a*b)
    ,
    D[5]*t/4 + D[7]*t/4 + (D[4]*a*a*t + D[8]*b*b*t)/(3*a*b)
    ,
    -D[3]*t/4 + D[8]*t/4 + (D[5]*a*a*t - 2*D[6]*b*b*t)/(6*a*b)
    ,
    -D[5]*t/4 + D[7]*t/4 + (D[4]*a*a*t - 2*D[8]*b*b*t)/(6*a*b)
    ,
    D[2]*t/4 - D[6]*t/4 + (D[0]*b*b*t - 2*D[8]*a*a*t)/(6*a*b)
    ,
    D[1]*t/4 - D[8]*t/4 + (D[2]*b*b*t - 2*D[7]*a*a*t)/(6*a*b)
    ,
    D[2]*t/4 + D[6]*t/4 + (-D[0]*b*b*t - D[8]*a*a*t)/(6*a*b)
    ,
    D[1]*t/4 + D[8]*t/4 + (-D[2]*b*b*t - D[7]*a*a*t)/(6*a*b)
    ,
    -D[2]*t/4 + D[6]*t/4 + (-2*D[0]*b*b*t + D[8]*a*a*t)/(6*a*b)
    ,
    -D[1]*t/4 + D[8]*t/4 + (-2*D[2]*b*b*t + D[7]*a*a*t)/(6*a*b)
    ,
    -D[2]*t/4 - D[6]*t/4 + (D[0]*b*b*t + D[8]*a*a*t)/(3*a*b)
    ,
    -D[1]*t/4 - D[8]*t/4 + (D[2]*b*b*t + D[7]*a*a*t)/(3*a*b)
    ,
    -D[3]*t/4 + D[8]*t/4 + (-2*D[5]*a*a*t + D[6]*b*b*t)/(6*a*b)
    ,
    -D[5]*t/4 + D[7]*t/4 + (-2*D[4]*a*a*t + D[8]*b*b*t)/(6*a*b)
    ,
    D[3]*t/4 + D[8]*t/4 + (-D[5]*a*a*t - D[6]*b*b*t)/(6*a*b)
    ,
    D[5]*t/4 + D[7]*t/4 + (-D[4]*a*a*t - D[8]*b*b*t)/(6*a*b)
    ,
    D[3]*t/4 - D[8]*t/4 + (D[5]*a*a*t - 2*D[6]*b*b*t)/(6*a*b)
    ,
    D[5]*t/4 - D[7]*t/4 + (D[4]*a*a*t - 2*D[8]*b*b*t)/(6*a*b)
    ,
    -D[3]*t/4 - D[8]*t/4 + (D[5]*a*a*t + D[6]*b*b*t)/(3*a*b)
    ,
    -D[5]*t/4 - D[7]*t/4 + (D[4]*a*a*t + D[8]*b*b*t)/(3*a*b)
    };

    return k;
}

std::vector<double> Q4S::get_DB(const std::vector<double>& D, const gp_Pnt& point) const{

    const gp_Pnt p = this->normalize(point);

    const double xi = p.X();
    const double eta = p.Y();

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

std::vector<double> Q4S::get_Nf(const double t, const std::vector<gp_Pnt>& points) const{

    const gp_Pnt p0 = this->normalize(points[0]);
    const gp_Pnt p1 = this->normalize(points[1]);

    const std::array<double, 2> x{p0.X(), p1.X()};
    const std::array<double, 2> y{p0.Y(), p1.Y()};

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

std::vector<double> Q4S::get_phi_radial(const double t, const double beta, const double vp, const std::vector<double>& axis, const std::vector<double>& center, const double rho) const{
    const double ax = axis[0];
    const double ay = axis[1];
    const double cx = center[0];
    const double cy = center[1];
    std::vector<double> phi{
t*(a*a*(12*a*a*ax*ax*ax*ax + 12*a*a*ax*ax*ay*ay - 24*a*a*ax*ax + 12*a*a - 30*a*ax*ax*ax*ax*cx + 30*a*ax*ax*ax*ax*x0 - 30*a*ax*ax*ax*ay*cy + 30*a*ax*ax*ax*ay*y0 - 30*a*ax*ax*ay*ay*cx + 30*a*ax*ax*ay*ay*x0 + 60*a*ax*ax*cx - 60*a*ax*ax*x0 - 30*a*ax*ay*ay*ay*cy + 30*a*ax*ay*ay*ay*y0 + 60*a*ax*ay*cy - 60*a*ax*ay*y0 - 30*a*cx + 30*a*x0 + 30*ax*ax*ax*ax*cx*cx - 60*ax*ax*ax*ax*cx*x0 + 30*ax*ax*ax*ax*x0*x0 + 60*ax*ax*ax*ay*cx*cy - 60*ax*ax*ax*ay*cx*y0 - 60*ax*ax*ax*ay*cy*x0 + 60*ax*ax*ax*ay*x0*y0 + 30*ax*ax*ay*ay*cx*cx - 60*ax*ax*ay*ay*cx*x0 + 30*ax*ax*ay*ay*cy*cy - 60*ax*ax*ay*ay*cy*y0 + 30*ax*ax*ay*ay*x0*x0 + 30*ax*ax*ay*ay*y0*y0 - 60*ax*ax*cx*cx + 120*ax*ax*cx*x0 - 60*ax*ax*x0*x0 + 60*ax*ay*ay*ay*cx*cy - 60*ax*ay*ay*ay*cx*y0 - 60*ax*ay*ay*ay*cy*x0 + 60*ax*ay*ay*ay*x0*y0 - 120*ax*ay*cx*cy + 120*ax*ay*cx*y0 + 120*ax*ay*cy*x0 - 120*ax*ay*x0*y0 + 30*ay*ay*ay*ay*cy*cy - 60*ay*ay*ay*ay*cy*y0 + 30*ay*ay*ay*ay*y0*y0 - 60*ay*ay*cy*cy + 120*ay*ay*cy*y0 - 60*ay*ay*y0*y0 + 30*cx*cx - 60*cx*x0 + 30*cy*cy - 60*cy*y0 + 30*x0*x0 + 30*y0*y0) + 5*a*b*(6*a*a*ax*ax*ax*ay + 6*a*a*ax*ay*ay*ay - 3*a*a*ax*ay*vp - 12*a*a*ax*ay + 8*a*ax*ax*ax*ax*b - 12*a*ax*ax*ax*ay*cx + 12*a*ax*ax*ax*ay*x0 + 16*a*ax*ax*ay*ay*b - 12*a*ax*ax*ay*ay*cy + 12*a*ax*ax*ay*ay*y0 + 4*a*ax*ax*b*vp - 16*a*ax*ax*b - 12*a*ax*ay*ay*ay*cx + 12*a*ax*ay*ay*ay*x0 + 6*a*ax*ay*cx*vp + 24*a*ax*ay*cx - 6*a*ax*ay*vp*x0 - 24*a*ax*ay*x0 + 8*a*ay*ay*ay*ay*b - 12*a*ay*ay*ay*ay*cy + 12*a*ay*ay*ay*ay*y0 + 4*a*ay*ay*b*vp - 16*a*ay*ay*b + 6*a*ay*ay*cy*vp + 24*a*ay*ay*cy - 6*a*ay*ay*vp*y0 - 24*a*ay*ay*y0 + 8*a*b*beta*rho - 8*a*b*vp + 16*a*b - 6*a*cy*vp - 12*a*cy + 6*a*vp*y0 + 12*a*y0 - 12*ax*ax*ax*ax*b*cx + 12*ax*ax*ax*ax*b*x0 + 6*ax*ax*ax*ay*b*b - 12*ax*ax*ax*ay*b*cy + 12*ax*ax*ax*ay*b*y0 - 12*ax*ax*ay*ay*b*cx + 12*ax*ax*ay*ay*b*x0 + 6*ax*ax*b*cx*vp + 24*ax*ax*b*cx - 6*ax*ax*b*vp*x0 - 24*ax*ax*b*x0 + 6*ax*ay*ay*ay*b*b - 12*ax*ay*ay*ay*b*cy + 12*ax*ay*ay*ay*b*y0 - 3*ax*ay*b*b*vp - 12*ax*ay*b*b + 6*ax*ay*b*cy*vp + 24*ax*ay*b*cy - 6*ax*ay*b*vp*y0 - 24*ax*ay*b*y0 - 6*b*cx*vp - 12*b*cx + 6*b*vp*x0 + 12*b*x0) + b*b*(30*ax*ax*ax*ax*cx*cx - 60*ax*ax*ax*ax*cx*x0 + 30*ax*ax*ax*ax*x0*x0 - 30*ax*ax*ax*ay*b*cx + 30*ax*ax*ax*ay*b*x0 + 60*ax*ax*ax*ay*cx*cy - 60*ax*ax*ax*ay*cx*y0 - 60*ax*ax*ax*ay*cy*x0 + 60*ax*ax*ax*ay*x0*y0 + 12*ax*ax*ay*ay*b*b - 30*ax*ax*ay*ay*b*cy + 30*ax*ax*ay*ay*b*y0 + 30*ax*ax*ay*ay*cx*cx - 60*ax*ax*ay*ay*cx*x0 + 30*ax*ax*ay*ay*cy*cy - 60*ax*ax*ay*ay*cy*y0 + 30*ax*ax*ay*ay*x0*x0 + 30*ax*ax*ay*ay*y0*y0 - 60*ax*ax*cx*cx + 120*ax*ax*cx*x0 - 60*ax*ax*x0*x0 - 30*ax*ay*ay*ay*b*cx + 30*ax*ay*ay*ay*b*x0 + 60*ax*ay*ay*ay*cx*cy - 60*ax*ay*ay*ay*cx*y0 - 60*ax*ay*ay*ay*cy*x0 + 60*ax*ay*ay*ay*x0*y0 + 60*ax*ay*b*cx - 60*ax*ay*b*x0 - 120*ax*ay*cx*cy + 120*ax*ay*cx*y0 + 120*ax*ay*cy*x0 - 120*ax*ay*x0*y0 + 12*ay*ay*ay*ay*b*b - 30*ay*ay*ay*ay*b*cy + 30*ay*ay*ay*ay*b*y0 + 30*ay*ay*ay*ay*cy*cy - 60*ay*ay*ay*ay*cy*y0 + 30*ay*ay*ay*ay*y0*y0 - 24*ay*ay*b*b + 60*ay*ay*b*cy - 60*ay*ay*b*y0 - 60*ay*ay*cy*cy + 120*ay*ay*cy*y0 - 60*ay*ay*y0*y0 + 12*b*b - 30*b*cy + 30*b*y0 + 30*cx*cx - 60*cx*x0 + 30*cy*cy - 60*cy*y0 + 30*x0*x0 + 30*y0*y0))/(90*a*b)
,
t*(a*a*(18*a*a*ax*ax*ax*ax + 18*a*a*ax*ax*ay*ay - 36*a*a*ax*ax + 18*a*a - 30*a*ax*ax*ax*ax*cx + 30*a*ax*ax*ax*ax*x0 - 30*a*ax*ax*ax*ay*cy + 30*a*ax*ax*ax*ay*y0 - 30*a*ax*ax*ay*ay*cx + 30*a*ax*ax*ay*ay*x0 + 60*a*ax*ax*cx - 60*a*ax*ax*x0 - 30*a*ax*ay*ay*ay*cy + 30*a*ax*ay*ay*ay*y0 + 60*a*ax*ay*cy - 60*a*ax*ay*y0 - 30*a*cx + 30*a*x0 + 15*ax*ax*ax*ax*cx*cx - 30*ax*ax*ax*ax*cx*x0 + 15*ax*ax*ax*ax*x0*x0 + 30*ax*ax*ax*ay*cx*cy - 30*ax*ax*ax*ay*cx*y0 - 30*ax*ax*ax*ay*cy*x0 + 30*ax*ax*ax*ay*x0*y0 + 15*ax*ax*ay*ay*cx*cx - 30*ax*ax*ay*ay*cx*x0 + 15*ax*ax*ay*ay*cy*cy - 30*ax*ax*ay*ay*cy*y0 + 15*ax*ax*ay*ay*x0*x0 + 15*ax*ax*ay*ay*y0*y0 - 30*ax*ax*cx*cx + 60*ax*ax*cx*x0 - 30*ax*ax*x0*x0 + 30*ax*ay*ay*ay*cx*cy - 30*ax*ay*ay*ay*cx*y0 - 30*ax*ay*ay*ay*cy*x0 +30*ax*ay*ay*ay*x0*y0 - 60*ax*ay*cx*cy + 60*ax*ay*cx*y0 + 60*ax*ay*cy*x0 - 60*ax*ay*x0*y0 + 15*ay*ay*ay*ay*cy*cy - 30*ay*ay*ay*ay*cy*y0 + 15*ay*ay*ay*ay*y0*y0 - 30*ay*ay*cy*cy + 60*ay*ay*cy*y0 - 30*ay*ay*y0*y0 + 15*cx*cx - 30*cx*x0 + 15*cy*cy - 30*cy*y0 + 15*x0*x0 + 15*y0*y0) + 5*a*b*(6*a*a*ax*ax*ax*ay + 6*a*a*ax*ay*ay*ay - 3*a*a*ax*ay*vp - 12*a*a*ax*ay - 8*a*ax*ax*ax*ax*b - 6*a*ax*ax*ax*ay*cx + 6*a*ax*ax*ax*ay*x0 - 4*a*ax*ax*ay*ay*b - 6*a*ax*ax*ay*ay*cy + 6*a*ax*ax*ay*ay*y0 + 8*a*ax*ax*b*vp + 16*a*ax*ax*b - 6*a*ax*ay*ay*ay*cx + 6*a*ax*ay*ay*ay*x0 + 3*a*ax*ay*cx*vp + 12*a*ax*ay*cx - 3*a*ax*ay*vp*x0 - 12*a*ax*ay*x0 + 4*a*ay*ay*ay*ay*b - 6*a*ay*ay*ay*ay*cy + 6*a*ay*ay*ay*ay*y0 + 2*a*ay*ay*b*vp - 8*a*ay*ay*b + 3*a*ay*ay*cy*vp + 12*a*ay*ay*cy - 3*a*ay*ay*vp*y0 - 12*a*ay*ay*y0 + 4*a*b*beta*rho - 10*a*b*vp - 4*a*b - 3*a*cy*vp - 6*a*cy + 3*a*vp*y0 + 6*a*y0 + 12*ax*ax*ax*ax*b*cx - 12*ax*ax*ax*ax*b*x0 - 6*ax*ax*ax*ay*b*b + 12*ax*ax*ax*ay*b*cy - 12*ax*ax*ax*ay*b*y0 + 12*ax*ax*ay*ay*b*cx - 12*ax*ax*ay*ay*b*x0 - 6*ax*ax*b*cx*vp - 24*ax*ax*b*cx + 6*ax*ax*b*vp*x0 + 24*ax*ax*b*x0 - 6*ax*ay*ay*ay*b*b + 12*ax*ay*ay*ay*b*cy - 12*ax*ay*ay*ay*b*y0 + 3*ax*ay*b*b*vp + 12*ax*ay*b*b - 6*ax*ay*b*cy*vp - 24*ax*ay*b*cy + 6*ax*ay*b*vp*y0 + 24*ax*ay*b*y0 + 6*b*cx*vp + 12*b*cx - 6*b*vp*x0 - 12*b*x0) + b*b*(-30*ax*ax*ax*ax*cx*cx + 60*ax*ax*ax*ax*cx*x0 - 30*ax*ax*ax*ax*x0*x0 + 30*ax*ax*ax*ay*b*cx - 30*ax*ax*ax*ay*b*x0 - 60*ax*ax*ax*ay*cx*cy + 60*ax*ax*ax*ay*cx*y0 + 60*ax*ax*ax*ay*cy*x0 - 60*ax*ax*ax*ay*x0*y0 - 12*ax*ax*ay*ay*b*b + 30*ax*ax*ay*ay*b*cy - 30*ax*ax*ay*ay*b*y0 - 30*ax*ax*ay*ay*cx*cx + 60*ax*ax*ay*ay*cx*x0 - 30*ax*ax*ay*ay*cy*cy + 60*ax*ax*ay*ay*cy*y0 - 30*ax*ax*ay*ay*x0*x0 - 30*ax*ax*ay*ay*y0*y0 + 60*ax*ax*cx*cx - 120*ax*ax*cx*x0 + 60*ax*ax*x0*x0 + 30*ax*ay*ay*ay*b*cx - 30*ax*ay*ay*ay*b*x0 - 60*ax*ay*ay*ay*cx*cy + 60*ax*ay*ay*ay*cx*y0 + 60*ax*ay*ay*ay*cy*x0 - 60*ax*ay*ay*ay*x0*y0 - 60*ax*ay*b*cx + 60*ax*ay*b*x0 + 120*ax*ay*cx*cy - 120*ax*ay*cx*y0 - 120*ax*ay*cy*x0 + 120*ax*ay*x0*y0 - 12*ay*ay*ay*ay*b*b + 30*ay*ay*ay*ay*b*cy - 30*ay*ay*ay*ay*b*y0 - 30*ay*ay*ay*ay*cy*cy + 60*ay*ay*ay*ay*cy*y0 - 30*ay*ay*ay*ay*y0*y0 + 24*ay*ay*b*b - 60*ay*ay*b*cy + 60*ay*ay*b*y0 + 60*ay*ay*cy*cy - 120*ay*ay*cy*y0 + 60*ay*ay*y0*y0 - 12*b*b + 30*b*cy - 30*b*y0 - 30*cx*cx + 60*cx*x0 - 30*cy*cy + 60*cy*y0 - 30*x0*x0 - 30*y0*y0))/(90*a*b)
,
t*(a*a*(-18*a*a*ax*ax*ax*ax - 18*a*a*ax*ax*ay*ay + 36*a*a*ax*ax - 18*a*a + 30*a*ax*ax*ax*ax*cx - 30*a*ax*ax*ax*ax*x0 + 30*a*ax*ax*ax*ay*cy - 30*a*ax*ax*ax*ay*y0 + 30*a*ax*ax*ay*ay*cx - 30*a*ax*ax*ay*ay*x0 - 60*a*ax*ax*cx + 60*a*ax*ax*x0 + 30*a*ax*ay*ay*ay*cy - 30*a*ax*ay*ay*ay*y0 - 60*a*ax*ay*cy + 60*a*ax*ay*y0 + 30*a*cx - 30*a*x0 - 15*ax*ax*ax*ax*cx*cx + 30*ax*ax*ax*ax*cx*x0 - 15*ax*ax*ax*ax*x0*x0 - 30*ax*ax*ax*ay*cx*cy + 30*ax*ax*ax*ay*cx*y0 + 30*ax*ax*ax*ay*cy*x0 - 30*ax*ax*ax*ay*x0*y0 - 15*ax*ax*ay*ay*cx*cx + 30*ax*ax*ay*ay*cx*x0 - 15*ax*ax*ay*ay*cy*cy +30*ax*ax*ay*ay*cy*y0 - 15*ax*ax*ay*ay*x0*x0 - 15*ax*ax*ay*ay*y0*y0 + 30*ax*ax*cx*cx - 60*ax*ax*cx*x0 + 30*ax*ax*x0*x0 - 30*ax*ay*ay*ay*cx*cy + 30*ax*ay*ay*ay*cx*y0 + 30*ax*ay*ay*ay*cy*x0 - 30*ax*ay*ay*ay*x0*y0 + 60*ax*ay*cx*cy - 60*ax*ay*cx*y0 - 60*ax*ay*cy*x0 + 60*ax*ay*x0*y0 - 15*ay*ay*ay*ay*cy*cy + 30*ay*ay*ay*ay*cy*y0 - 15*ay*ay*ay*ay*y0*y0 + 30*ay*ay*cy*cy - 60*ay*ay*cy*y0 + 30*ay*ay*y0*y0 - 15*cx*cx + 30*cx*x0 - 15*cy*cy + 30*cy*y0 - 15*x0*x0 - 15*y0*y0) + 5*a*b*(-6*a*a*ax*ax*ax*ay - 6*a*a*ax*ay*ay*ay + 3*a*a*ax*ay*vp + 12*a*a*ax*ay - 4*a*ax*ax*ax*ax*b + 6*a*ax*ax*ax*ay*cx - 6*a*ax*ax*ax*ay*x0 - 8*a*ax*ax*ay*ay*b + 6*a*ax*ax*ay*ay*cy - 6*a*ax*ax*ay*ay*y0 + 4*a*ax*ax*b*vp + 8*a*ax*ax*b + 6*a*ax*ay*ay*ay*cx - 6*a*ax*ay*ay*ay*x0 - 3*a*ax*ay*cx*vp - 12*a*ax*ay*cx + 3*a*ax*ay*vp*x0 + 12*a*ax*ay*x0 - 4*a*ay*ay*ay*ay*b + 6*a*ay*ay*ay*ay*cy - 6*a*ay*ay*ay*ay*y0 + 4*a*ay*ay*b*vp + 8*a*ay*ay*b - 3*a*ay*ay*cy*vp - 12*a*ay*ay*cy + 3*a*ay*ay*vp*y0 + 12*a*ay*ay*y0 + 2*a*b*beta*rho - 8*a*b*vp - 8*a*b + 3*a*cy*vp + 6*a*cy - 3*a*vp*y0 - 6*a*y0 + 6*ax*ax*ax*ax*b*cx - 6*ax*ax*ax*ax*b*x0 - 6*ax*ax*ax*ay*b*b + 6*ax*ax*ax*ay*b*cy - 6*ax*ax*ax*ay*b*y0 + 6*ax*ax*ay*ay*b*cx - 6*ax*ax*ay*ay*b*x0 - 3*ax*ax*b*cx*vp - 12*ax*ax*b*cx + 3*ax*ax*b*vp*x0 + 12*ax*ax*b*x0 - 6*ax*ay*ay*ay*b*b + 6*ax*ay*ay*ay*b*cy - 6*ax*ay*ay*ay*b*y0 + 3*ax*ay*b*b*vp + 12*ax*ay*b*b - 3*ax*ay*b*cy*vp - 12*ax*ay*b*cy + 3*ax*ay*b*vp*y0 + 12*ax*ay*b*y0 + 3*b*cx*vp + 6*b*cx - 3*b*vp*x0 - 6*b*x0) + b*b*(-15*ax*ax*ax*ax*cx*cx + 30*ax*ax*ax*ax*cx*x0 - 15*ax*ax*ax*ax*x0*x0 + 30*ax*ax*ax*ay*b*cx - 30*ax*ax*ax*ay*b*x0 - 30*ax*ax*ax*ay*cx*cy + 30*ax*ax*ax*ay*cx*y0 + 30*ax*ax*ax*ay*cy*x0 - 30*ax*ax*ax*ay*x0*y0 - 18*ax*ax*ay*ay*b*b + 30*ax*ax*ay*ay*b*cy - 30*ax*ax*ay*ay*b*y0 - 15*ax*ax*ay*ay*cx*cx + 30*ax*ax*ay*ay*cx*x0 - 15*ax*ax*ay*ay*cy*cy + 30*ax*ax*ay*ay*cy*y0 - 15*ax*ax*ay*ay*x0*x0 - 15*ax*ax*ay*ay*y0*y0 + 30*ax*ax*cx*cx - 60*ax*ax*cx*x0 + 30*ax*ax*x0*x0 + 30*ax*ay*ay*ay*b*cx - 30*ax*ay*ay*ay*b*x0 - 30*ax*ay*ay*ay*cx*cy + 30*ax*ay*ay*ay*cx*y0 + 30*ax*ay*ay*ay*cy*x0 - 30*ax*ay*ay*ay*x0*y0 - 60*ax*ay*b*cx + 60*ax*ay*b*x0 + 60*ax*ay*cx*cy - 60*ax*ay*cx*y0 - 60*ax*ay*cy*x0 + 60*ax*ay*x0*y0 - 18*ay*ay*ay*ay*b*b + 30*ay*ay*ay*ay*b*cy - 30*ay*ay*ay*ay*b*y0 - 15*ay*ay*ay*ay*cy*cy + 30*ay*ay*ay*ay*cy*y0 - 15*ay*ay*ay*ay*y0*y0 + 36*ay*ay*b*b - 60*ay*ay*b*cy + 60*ay*ay*b*y0 + 30*ay*ay*cy*cy - 60*ay*ay*cy*y0 + 30*ay*ay*y0*y0 - 18*b*b + 30*b*cy - 30*b*y0 - 15*cx*cx + 30*cx*x0 - 15*cy*cy + 30*cy*y0 - 15*x0*x0 - 15*y0*y0))/(90*a*b)
,
t*(a*a*(-12*a*a*ax*ax*ax*ax - 12*a*a*ax*ax*ay*ay + 24*a*a*ax*ax - 12*a*a + 30*a*ax*ax*ax*ax*cx - 30*a*ax*ax*ax*ax*x0 + 30*a*ax*ax*ax*ay*cy - 30*a*ax*ax*ax*ay*y0 + 30*a*ax*ax*ay*ay*cx - 30*a*ax*ax*ay*ay*x0 - 60*a*ax*ax*cx + 60*a*ax*ax*x0 + 30*a*ax*ay*ay*ay*cy - 30*a*ax*ay*ay*ay*y0 - 60*a*ax*ay*cy + 60*a*ax*ay*y0 + 30*a*cx - 30*a*x0 - 30*ax*ax*ax*ax*cx*cx + 60*ax*ax*ax*ax*cx*x0 - 30*ax*ax*ax*ax*x0*x0 - 60*ax*ax*ax*ay*cx*cy + 60*ax*ax*ax*ay*cx*y0 + 60*ax*ax*ax*ay*cy*x0 - 60*ax*ax*ax*ay*x0*y0 - 30*ax*ax*ay*ay*cx*cx + 60*ax*ax*ay*ay*cx*x0 - 30*ax*ax*ay*ay*cy*cy +60*ax*ax*ay*ay*cy*y0 - 30*ax*ax*ay*ay*x0*x0 - 30*ax*ax*ay*ay*y0*y0 + 60*ax*ax*cx*cx - 120*ax*ax*cx*x0 + 60*ax*ax*x0*x0 - 60*ax*ay*ay*ay*cx*cy + 60*ax*ay*ay*ay*cx*y0 + 60*ax*ay*ay*ay*cy*x0- 60*ax*ay*ay*ay*x0*y0 + 120*ax*ay*cx*cy - 120*ax*ay*cx*y0 - 120*ax*ay*cy*x0 + 120*ax*ay*x0*y0 - 30*ay*ay*ay*ay*cy*cy + 60*ay*ay*ay*ay*cy*y0 - 30*ay*ay*ay*ay*y0*y0 + 60*ay*ay*cy*cy - 120*ay*ay*cy*y0 + 60*ay*ay*y0*y0 - 30*cx*cx + 60*cx*x0 - 30*cy*cy + 60*cy*y0 - 30*x0*x0 - 30*y0*y0) + 5*a*b*(-6*a*a*ax*ax*ax*ay - 6*a*a*ax*ay*ay*ay + 3*a*a*ax*ay*vp + 12*a*a*ax*ay + 4*a*ax*ax*ax*ax*b + 12*a*ax*ax*ax*ay*cx - 12*a*ax*ax*ax*ay*x0 - 4*a*ax*ax*ay*ay*b + 12*a*ax*ax*ay*ay*cy - 12*a*ax*ax*ay*ay*y0 + 2*a*ax*ax*b*vp - 8*a*ax*ax*b + 12*a*ax*ay*ay*ay*cx - 12*a*ax*ay*ay*ay*x0 - 6*a*ax*ay*cx*vp - 24*a*ax*ay*cx + 6*a*ax*ay*vp*x0 + 24*a*ax*ay*x0 - 8*a*ay*ay*ay*ay*b + 12*a*ay*ay*ay*ay*cy - 12*a*ay*ay*ay*ay*y0 + 8*a*ay*ay*b*vp + 16*a*ay*ay*b - 6*a*ay*ay*cy*vp - 24*a*ay*ay*cy + 6*a*ay*ay*vp*y0 + 24*a*ay*ay*y0 + 4*a*b*beta*rho - 10*a*b*vp - 4*a*b + 6*a*cy*vp + 12*a*cy - 6*a*vp*y0 - 12*a*y0 - 6*ax*ax*ax*ax*b*cx + 6*ax*ax*ax*ax*b*x0 + 6*ax*ax*ax*ay*b*b - 6*ax*ax*ax*ay*b*cy + 6*ax*ax*ax*ay*b*y0 - 6*ax*ax*ay*ay*b*cx + 6*ax*ax*ay*ay*b*x0 + 3*ax*ax*b*cx*vp + 12*ax*ax*b*cx - 3*ax*ax*b*vp*x0 - 12*ax*ax*b*x0 + 6*ax*ay*ay*ay*b*b - 6*ax*ay*ay*ay*b*cy + 6*ax*ay*ay*ay*b*y0 - 3*ax*ay*b*b*vp - 12*ax*ay*b*b + 3*ax*ay*b*cy*vp + 12*ax*ay*b*cy - 3*ax*ay*b*vp*y0 - 12*ax*ay*b*y0 - 3*b*cx*vp - 6*b*cx + 3*b*vp*x0 + 6*b*x0) + b*b*(15*ax*ax*ax*ax*cx*cx - 30*ax*ax*ax*ax*cx*x0 + 15*ax*ax*ax*ax*x0*x0 - 30*ax*ax*ax*ay*b*cx + 30*ax*ax*ax*ay*b*x0 + 30*ax*ax*ax*ay*cx*cy - 30*ax*ax*ax*ay*cx*y0 - 30*ax*ax*ax*ay*cy*x0 + 30*ax*ax*ax*ay*x0*y0 + 18*ax*ax*ay*ay*b*b - 30*ax*ax*ay*ay*b*cy + 30*ax*ax*ay*ay*b*y0 + 15*ax*ax*ay*ay*cx*cx - 30*ax*ax*ay*ay*cx*x0 + 15*ax*ax*ay*ay*cy*cy - 30*ax*ax*ay*ay*cy*y0 + 15*ax*ax*ay*ay*x0*x0 + 15*ax*ax*ay*ay*y0*y0 - 30*ax*ax*cx*cx + 60*ax*ax*cx*x0 - 30*ax*ax*x0*x0 - 30*ax*ay*ay*ay*b*cx + 30*ax*ay*ay*ay*b*x0 + 30*ax*ay*ay*ay*cx*cy - 30*ax*ay*ay*ay*cx*y0 - 30*ax*ay*ay*ay*cy*x0 + 30*ax*ay*ay*ay*x0*y0 + 60*ax*ay*b*cx - 60*ax*ay*b*x0 - 60*ax*ay*cx*cy + 60*ax*ay*cx*y0 + 60*ax*ay*cy*x0 - 60*ax*ay*x0*y0 + 18*ay*ay*ay*ay*b*b - 30*ay*ay*ay*ay*b*cy + 30*ay*ay*ay*ay*b*y0 + 15*ay*ay*ay*ay*cy*cy - 30*ay*ay*ay*ay*cy*y0 + 15*ay*ay*ay*ay*y0*y0 - 36*ay*ay*b*b + 60*ay*ay*b*cy - 60*ay*ay*b*y0 - 30*ay*ay*cy*cy + 60*ay*ay*cy*y0 - 30*ay*ay*y0*y0 + 18*b*b - 30*b*cy + 30*b*y0 + 15*cx*cx - 30*cx*x0 + 15*cy*cy - 30*cy*y0 + 15*x0*x0 + 15*y0*y0))/(90*a*b)
,
t*(a*a*(18*a*a*ax*ax*ax*ax + 18*a*a*ax*ax*ay*ay - 36*a*a*ax*ax + 18*a*a - 30*a*ax*ax*ax*ax*cx + 30*a*ax*ax*ax*ax*x0 - 30*a*ax*ax*ax*ay*cy + 30*a*ax*ax*ax*ay*y0 - 30*a*ax*ax*ay*ay*cx + 30*a*ax*ax*ay*ay*x0 + 60*a*ax*ax*cx - 60*a*ax*ax*x0 - 30*a*ax*ay*ay*ay*cy + 30*a*ax*ay*ay*ay*y0 + 60*a*ax*ay*cy - 60*a*ax*ay*y0 - 30*a*cx + 30*a*x0 + 15*ax*ax*ax*ax*cx*cx - 30*ax*ax*ax*ax*cx*x0 + 15*ax*ax*ax*ax*x0*x0 + 30*ax*ax*ax*ay*cx*cy - 30*ax*ax*ax*ay*cx*y0 - 30*ax*ax*ax*ay*cy*x0 + 30*ax*ax*ax*ay*x0*y0 + 15*ax*ax*ay*ay*cx*cx - 30*ax*ax*ay*ay*cx*x0 + 15*ax*ax*ay*ay*cy*cy - 30*ax*ax*ay*ay*cy*y0 + 15*ax*ax*ay*ay*x0*x0 + 15*ax*ax*ay*ay*y0*y0 - 30*ax*ax*cx*cx + 60*ax*ax*cx*x0 - 30*ax*ax*x0*x0 + 30*ax*ay*ay*ay*cx*cy - 30*ax*ay*ay*ay*cx*y0 - 30*ax*ay*ay*ay*cy*x0 +30*ax*ay*ay*ay*x0*y0 - 60*ax*ay*cx*cy + 60*ax*ay*cx*y0 + 60*ax*ay*cy*x0 - 60*ax*ay*x0*y0 + 15*ay*ay*ay*ay*cy*cy - 30*ay*ay*ay*ay*cy*y0 + 15*ay*ay*ay*ay*y0*y0 - 30*ay*ay*cy*cy + 60*ay*ay*cy*y0 - 30*ay*ay*y0*y0 + 15*cx*cx - 30*cx*x0 + 15*cy*cy - 30*cy*y0 + 15*x0*x0 + 15*y0*y0) + 5*a*b*(6*a*a*ax*ax*ax*ay + 6*a*a*ax*ay*ay*ay - 3*a*a*ax*ay*vp - 12*a*a*ax*ay - 8*a*ax*ax*ax*ax*b - 6*a*ax*ax*ax*ay*cx + 6*a*ax*ax*ax*ay*x0 - 4*a*ax*ax*ay*ay*b - 6*a*ax*ax*ay*ay*cy + 6*a*ax*ax*ay*ay*y0 - 4*a*ax*ax*b*vp + 16*a*ax*ax*b - 6*a*ax*ay*ay*ay*cx + 6*a*ax*ay*ay*ay*x0 + 3*a*ax*ay*cx*vp + 12*a*ax*ay*cx - 3*a*ax*ay*vp*x0 - 12*a*ax*ay*x0 + 4*a*ay*ay*ay*ay*b - 6*a*ay*ay*ay*ay*cy + 6*a*ay*ay*ay*ay*y0 + 2*a*ay*ay*b*vp - 8*a*ay*ay*b + 3*a*ay*ay*cy*vp + 12*a*ay*ay*cy - 3*a*ay*ay*vp*y0 - 12*a*ay*ay*y0 + 4*a*b*beta*rho + 2*a*b*vp - 4*a*b - 3*a*cy*vp - 6*a*cy + 3*a*vp*y0 + 6*a*y0 + 12*ax*ax*ax*ax*b*cx - 12*ax*ax*ax*ax*b*x0 - 6*ax*ax*ax*ay*b*b + 12*ax*ax*ax*ay*b*cy - 12*ax*ax*ax*ay*b*y0 + 12*ax*ax*ay*ay*b*cx - 12*ax*ax*ay*ay*b*x0 + 6*ax*ax*b*cx*vp - 24*ax*ax*b*cx - 6*ax*ax*b*vp*x0 + 24*ax*ax*b*x0 - 6*ax*ay*ay*ay*b*b + 12*ax*ay*ay*ay*b*cy - 12*ax*ay*ay*ay*b*y0 - 3*ax*ay*b*b*vp + 12*ax*ay*b*b + 6*ax*ay*b*cy*vp - 24*ax*ay*b*cy - 6*ax*ay*b*vp*y0 + 24*ax*ay*b*y0 - 6*b*cx*vp + 12*b*cx + 6*b*vp*x0 - 12*b*x0) + b*b*(-30*ax*ax*ax*ax*cx*cx + 60*ax*ax*ax*ax*cx*x0 - 30*ax*ax*ax*ax*x0*x0 + 30*ax*ax*ax*ay*b*cx - 30*ax*ax*ax*ay*b*x0 - 60*ax*ax*ax*ay*cx*cy + 60*ax*ax*ax*ay*cx*y0 + 60*ax*ax*ax*ay*cy*x0 - 60*ax*ax*ax*ay*x0*y0 - 12*ax*ax*ay*ay*b*b + 30*ax*ax*ay*ay*b*cy - 30*ax*ax*ay*ay*b*y0 - 30*ax*ax*ay*ay*cx*cx + 60*ax*ax*ay*ay*cx*x0 - 30*ax*ax*ay*ay*cy*cy + 60*ax*ax*ay*ay*cy*y0 - 30*ax*ax*ay*ay*x0*x0 - 30*ax*ax*ay*ay*y0*y0 + 60*ax*ax*cx*cx - 120*ax*ax*cx*x0 + 60*ax*ax*x0*x0 + 30*ax*ay*ay*ay*b*cx - 30*ax*ay*ay*ay*b*x0 - 60*ax*ay*ay*ay*cx*cy + 60*ax*ay*ay*ay*cx*y0 + 60*ax*ay*ay*ay*cy*x0 - 60*ax*ay*ay*ay*x0*y0 - 60*ax*ay*b*cx + 60*ax*ay*b*x0 + 120*ax*ay*cx*cy - 120*ax*ay*cx*y0 - 120*ax*ay*cy*x0 + 120*ax*ay*x0*y0 - 12*ay*ay*ay*ay*b*b + 30*ay*ay*ay*ay*b*cy - 30*ay*ay*ay*ay*b*y0 - 30*ay*ay*ay*ay*cy*cy + 60*ay*ay*ay*ay*cy*y0 - 30*ay*ay*ay*ay*y0*y0 + 24*ay*ay*b*b - 60*ay*ay*b*cy + 60*ay*ay*b*y0 + 60*ay*ay*cy*cy - 120*ay*ay*cy*y0 + 60*ay*ay*y0*y0 - 12*b*b + 30*b*cy - 30*b*y0 - 30*cx*cx + 60*cx*x0 - 30*cy*cy + 60*cy*y0 - 30*x0*x0 - 30*y0*y0))/(90*a*b)
,
t*(a*a*(72*a*a*ax*ax*ax*ax + 72*a*a*ax*ax*ay*ay - 144*a*a*ax*ax + 72*a*a - 90*a*ax*ax*ax*ax*cx + 90*a*ax*ax*ax*ax*x0 - 90*a*ax*ax*ax*ay*cy + 90*a*ax*ax*ax*ay*y0 - 90*a*ax*ax*ay*ay*cx + 90*a*ax*ax*ay*ay*x0 + 180*a*ax*ax*cx - 180*a*ax*ax*x0 - 90*a*ax*ay*ay*ay*cy + 90*a*ax*ay*ay*ay*y0 + 180*a*ax*ay*cy - 180*a*ax*ay*y0 - 90*a*cx + 90*a*x0 + 30*ax*ax*ax*ax*cx*cx - 60*ax*ax*ax*ax*cx*x0 + 30*ax*ax*ax*ax*x0*x0 + 60*ax*ax*ax*ay*cx*cy - 60*ax*ax*ax*ay*cx*y0 - 60*ax*ax*ax*ay*cy*x0 + 60*ax*ax*ax*ay*x0*y0 + 30*ax*ax*ay*ay*cx*cx - 60*ax*ax*ay*ay*cx*x0 + 30*ax*ax*ay*ay*cy*cy - 60*ax*ax*ay*ay*cy*y0 + 30*ax*ax*ay*ay*x0*x0 + 30*ax*ax*ay*ay*y0*y0 - 60*ax*ax*cx*cx + 120*ax*ax*cx*x0 - 60*ax*ax*x0*x0 + 60*ax*ay*ay*ay*cx*cy - 60*ax*ay*ay*ay*cx*y0 - 60*ax*ay*ay*ay*cy*x0 + 60*ax*ay*ay*ay*x0*y0 - 120*ax*ay*cx*cy + 120*ax*ay*cx*y0 + 120*ax*ay*cy*x0 - 120*ax*ay*x0*y0 + 30*ay*ay*ay*ay*cy*cy - 60*ay*ay*ay*ay*cy*y0 + 30*ay*ay*ay*ay*y0*y0 - 60*ay*ay*cy*cy + 120*ay*ay*cy*y0 - 60*ay*ay*y0*y0 + 30*cx*cx - 60*cx*x0 + 30*cy*cy - 60*cy*y0 + 30*x0*x0 + 30*y0*y0) + 5*a*b*(18*a*a*ax*ax*ax*ay + 18*a*a*ax*ay*ay*ay - 9*a*a*ax*ay*vp - 36*a*a*ax*ay + 8*a*ax*ax*ax*ax*b - 12*a*ax*ax*ax*ay*cx + 12*a*ax*ax*ax*ay*x0 + 16*a*ax*ax*ay*ay*b - 12*a*ax*ax*ay*ay*cy + 12*a*ax*ax*ay*ay*y0 + 16*a*ax*ax*b*vp - 16*a*ax*ax*b - 12*a*ax*ay*ay*ay*cx + 12*a*ax*ay*ay*ay*x0 + 6*a*ax*ay*cx*vp + 24*a*ax*ay*cx - 6*a*ax*ay*vp*x0 - 24*a*ax*ay*x0 + 8*a*ay*ay*ay*ay*b - 12*a*ay*ay*ay*ay*cy + 12*a*ay*ay*ay*ay*y0 + 4*a*ay*ay*b*vp - 16*a*ay*ay*b + 6*a*ay*ay*cy*vp + 24*a*ay*ay*cy - 6*a*ay*ay*vp*y0 - 24*a*ay*ay*y0 + 8*a*b*beta*rho - 20*a*b*vp + 16*a*b - 6*a*cy*vp - 12*a*cy + 6*a*vp*y0 + 12*a*y0 - 12*ax*ax*ax*ax*b*cx + 12*ax*ax*ax*ax*b*x0 + 6*ax*ax*ax*ay*b*b - 12*ax*ax*ax*ay*b*cy + 12*ax*ax*ax*ay*b*y0 - 12*ax*ax*ay*ay*b*cx + 12*ax*ax*ay*ay*b*x0 - 6*ax*ax*b*cx*vp + 24*ax*ax*b*cx + 6*ax*ax*b*vp*x0 - 24*ax*ax*b*x0 + 6*ax*ay*ay*ay*b*b- 12*ax*ay*ay*ay*b*cy + 12*ax*ay*ay*ay*b*y0 + 3*ax*ay*b*b*vp - 12*ax*ay*b*b - 6*ax*ay*b*cy*vp + 24*ax*ay*b*cy + 6*ax*ay*b*vp*y0 - 24*ax*ay*b*y0 + 6*b*cx*vp - 12*b*cx - 6*b*vp*x0 + 12*b*x0) + b*b*(30*ax*ax*ax*ax*cx*cx - 60*ax*ax*ax*ax*cx*x0 + 30*ax*ax*ax*ax*x0*x0 - 30*ax*ax*ax*ay*b*cx + 30*ax*ax*ax*ay*b*x0 + 60*ax*ax*ax*ay*cx*cy - 60*ax*ax*ax*ay*cx*y0 - 60*ax*ax*ax*ay*cy*x0+ 60*ax*ax*ax*ay*x0*y0 + 12*ax*ax*ay*ay*b*b - 30*ax*ax*ay*ay*b*cy + 30*ax*ax*ay*ay*b*y0 + 30*ax*ax*ay*ay*cx*cx - 60*ax*ax*ay*ay*cx*x0 + 30*ax*ax*ay*ay*cy*cy - 60*ax*ax*ay*ay*cy*y0 + 30*ax*ax*ay*ay*x0*x0 + 30*ax*ax*ay*ay*y0*y0 - 60*ax*ax*cx*cx + 120*ax*ax*cx*x0 - 60*ax*ax*x0*x0 - 30*ax*ay*ay*ay*b*cx + 30*ax*ay*ay*ay*b*x0 + 60*ax*ay*ay*ay*cx*cy - 60*ax*ay*ay*ay*cx*y0 - 60*ax*ay*ay*ay*cy*x0 + 60*ax*ay*ay*ay*x0*y0 + 60*ax*ay*b*cx - 60*ax*ay*b*x0 - 120*ax*ay*cx*cy + 120*ax*ay*cx*y0 + 120*ax*ay*cy*x0 - 120*ax*ay*x0*y0 + 12*ay*ay*ay*ay*b*b - 30*ay*ay*ay*ay*b*cy + 30*ay*ay*ay*ay*b*y0 + 30*ay*ay*ay*ay*cy*cy - 60*ay*ay*ay*ay*cy*y0 + 30*ay*ay*ay*ay*y0*y0 - 24*ay*ay*b*b + 60*ay*ay*b*cy - 60*ay*ay*b*y0 - 60*ay*ay*cy*cy + 120*ay*ay*cy*y0 - 60*ay*ay*y0*y0 + 12*b*b - 30*b*cy + 30*b*y0 + 30*cx*cx - 60*cx*x0 + 30*cy*cy - 60*cy*y0 + 30*x0*x0 + 30*y0*y0))/(90*a*b)
,
t*(a*a*(-72*a*a*ax*ax*ax*ax - 72*a*a*ax*ax*ay*ay + 144*a*a*ax*ax - 72*a*a + 90*a*ax*ax*ax*ax*cx - 90*a*ax*ax*ax*ax*x0 + 90*a*ax*ax*ax*ay*cy - 90*a*ax*ax*ax*ay*y0 + 90*a*ax*ax*ay*ay*cx - 90*a*ax*ax*ay*ay*x0 - 180*a*ax*ax*cx + 180*a*ax*ax*x0 + 90*a*ax*ay*ay*ay*cy - 90*a*ax*ay*ay*ay*y0 - 180*a*ax*ay*cy + 180*a*ax*ay*y0 + 90*a*cx - 90*a*x0 - 30*ax*ax*ax*ax*cx*cx + 60*ax*ax*ax*ax*cx*x0 - 30*ax*ax*ax*ax*x0*x0 - 60*ax*ax*ax*ay*cx*cy + 60*ax*ax*ax*ay*cx*y0 + 60*ax*ax*ax*ay*cy*x0 - 60*ax*ax*ax*ay*x0*y0 - 30*ax*ax*ay*ay*cx*cx + 60*ax*ax*ay*ay*cx*x0 - 30*ax*ax*ay*ay*cy*cy + 60*ax*ax*ay*ay*cy*y0 - 30*ax*ax*ay*ay*x0*x0 - 30*ax*ax*ay*ay*y0*y0 + 60*ax*ax*cx*cx - 120*ax*ax*cx*x0 + 60*ax*ax*x0*x0 - 60*ax*ay*ay*ay*cx*cy + 60*ax*ay*ay*ay*cx*y0 + 60*ax*ay*ay*ay*cy*x0 - 60*ax*ay*ay*ay*x0*y0 + 120*ax*ay*cx*cy - 120*ax*ay*cx*y0 - 120*ax*ay*cy*x0 + 120*ax*ay*x0*y0 - 30*ay*ay*ay*ay*cy*cy + 60*ay*ay*ay*ay*cy*y0 - 30*ay*ay*ay*ay*y0*y0 + 60*ay*ay*cy*cy -120*ay*ay*cy*y0 + 60*ay*ay*y0*y0 - 30*cx*cx + 60*cx*x0 - 30*cy*cy + 60*cy*y0 - 30*x0*x0 - 30*y0*y0) + 5*a*b*(-18*a*a*ax*ax*ax*ay - 18*a*a*ax*ay*ay*ay + 9*a*a*ax*ay*vp + 36*a*a*ax*ay + 4*a*ax*ax*ax*ax*b + 12*a*ax*ax*ax*ay*cx - 12*a*ax*ax*ax*ay*x0 - 4*a*ax*ax*ay*ay*b + 12*a*ax*ax*ay*ay*cy - 12*a*ax*ax*ay*ay*y0 + 8*a*ax*ax*b*vp - 8*a*ax*ax*b + 12*a*ax*ay*ay*ay*cx - 12*a*ax*ay*ay*ay*x0 - 6*a*ax*ay*cx*vp - 24*a*ax*ay*cx + 6*a*ax*ay*vp*x0 + 24*a*ax*ay*x0 - 8*a*ay*ay*ay*ay*b + 12*a*ay*ay*ay*ay*cy - 12*a*ay*ay*ay*ay*y0 + 8*a*ay*ay*b*vp + 16*a*ay*ay*b - 6*a*ay*ay*cy*vp - 24*a*ay*ay*cy + 6*a*ay*ay*vp*y0 + 24*a*ay*ay*y0 + 4*a*b*beta*rho - 16*a*b*vp - 4*a*b + 6*a*cy*vp + 12*a*cy - 6*a*vp*y0 - 12*a*y0 - 6*ax*ax*ax*ax*b*cx + 6*ax*ax*ax*ax*b*x0 + 6*ax*ax*ax*ay*b*b - 6*ax*ax*ax*ay*b*cy + 6*ax*ax*ax*ay*b*y0 - 6*ax*ax*ay*ay*b*cx + 6*ax*ax*ay*ay*b*x0 - 3*ax*ax*b*cx*vp + 12*ax*ax*b*cx + 3*ax*ax*b*vp*x0 - 12*ax*ax*b*x0 + 6*ax*ay*ay*ay*b*b - 6*ax*ay*ay*ay*b*cy + 6*ax*ay*ay*ay*b*y0 + 3*ax*ay*b*b*vp - 12*ax*ay*b*b - 3*ax*ay*b*cy*vp + 12*ax*ay*b*cy + 3*ax*ay*b*vp*y0 - 12*ax*ay*b*y0 + 3*b*cx*vp - 6*b*cx - 3*b*vp*x0 + 6*b*x0) + b*b*(15*ax*ax*ax*ax*cx*cx - 30*ax*ax*ax*ax*cx*x0 + 15*ax*ax*ax*ax*x0*x0 - 30*ax*ax*ax*ay*b*cx + 30*ax*ax*ax*ay*b*x0 + 30*ax*ax*ax*ay*cx*cy - 30*ax*ax*ax*ay*cx*y0 - 30*ax*ax*ax*ay*cy*x0 + 30*ax*ax*ax*ay*x0*y0 + 18*ax*ax*ay*ay*b*b - 30*ax*ax*ay*ay*b*cy + 30*ax*ax*ay*ay*b*y0 + 15*ax*ax*ay*ay*cx*cx - 30*ax*ax*ay*ay*cx*x0 + 15*ax*ax*ay*ay*cy*cy - 30*ax*ax*ay*ay*cy*y0 + 15*ax*ax*ay*ay*x0*x0 + 15*ax*ax*ay*ay*y0*y0 - 30*ax*ax*cx*cx + 60*ax*ax*cx*x0 - 30*ax*ax*x0*x0 - 30*ax*ay*ay*ay*b*cx + 30*ax*ay*ay*ay*b*x0 + 30*ax*ay*ay*ay*cx*cy - 30*ax*ay*ay*ay*cx*y0 - 30*ax*ay*ay*ay*cy*x0 + 30*ax*ay*ay*ay*x0*y0 + 60*ax*ay*b*cx - 60*ax*ay*b*x0 - 60*ax*ay*cx*cy + 60*ax*ay*cx*y0 + 60*ax*ay*cy*x0 - 60*ax*ay*x0*y0 + 18*ay*ay*ay*ay*b*b - 30*ay*ay*ay*ay*b*cy + 30*ay*ay*ay*ay*b*y0 + 15*ay*ay*ay*ay*cy*cy - 30*ay*ay*ay*ay*cy*y0 + 15*ay*ay*ay*ay*y0*y0 - 36*ay*ay*b*b + 60*ay*ay*b*cy - 60*ay*ay*b*y0 - 30*ay*ay*cy*cy + 60*ay*ay*cy*y0 - 30*ay*ay*y0*y0 + 18*b*b - 30*b*cy+ 30*b*y0 + 15*cx*cx - 30*cx*x0 + 15*cy*cy - 30*cy*y0 + 15*x0*x0 + 15*y0*y0))/(90*a*b)
,
t*(a*a*(-18*a*a*ax*ax*ax*ax - 18*a*a*ax*ax*ay*ay + 36*a*a*ax*ax - 18*a*a + 30*a*ax*ax*ax*ax*cx - 30*a*ax*ax*ax*ax*x0 + 30*a*ax*ax*ax*ay*cy - 30*a*ax*ax*ax*ay*y0 + 30*a*ax*ax*ay*ay*cx - 30*a*ax*ax*ay*ay*x0 - 60*a*ax*ax*cx + 60*a*ax*ax*x0 + 30*a*ax*ay*ay*ay*cy - 30*a*ax*ay*ay*ay*y0 - 60*a*ax*ay*cy + 60*a*ax*ay*y0 + 30*a*cx - 30*a*x0 - 15*ax*ax*ax*ax*cx*cx + 30*ax*ax*ax*ax*cx*x0 - 15*ax*ax*ax*ax*x0*x0 - 30*ax*ax*ax*ay*cx*cy + 30*ax*ax*ax*ay*cx*y0 + 30*ax*ax*ax*ay*cy*x0 - 30*ax*ax*ax*ay*x0*y0 - 15*ax*ax*ay*ay*cx*cx + 30*ax*ax*ay*ay*cx*x0 - 15*ax*ax*ay*ay*cy*cy +30*ax*ax*ay*ay*cy*y0 - 15*ax*ax*ay*ay*x0*x0 - 15*ax*ax*ay*ay*y0*y0 + 30*ax*ax*cx*cx - 60*ax*ax*cx*x0 + 30*ax*ax*x0*x0 - 30*ax*ay*ay*ay*cx*cy + 30*ax*ay*ay*ay*cx*y0 + 30*ax*ay*ay*ay*cy*x0 - 30*ax*ay*ay*ay*x0*y0 + 60*ax*ay*cx*cy - 60*ax*ay*cx*y0 - 60*ax*ay*cy*x0 + 60*ax*ay*x0*y0 - 15*ay*ay*ay*ay*cy*cy + 30*ay*ay*ay*ay*cy*y0 - 15*ay*ay*ay*ay*y0*y0 + 30*ay*ay*cy*cy - 60*ay*ay*cy*y0 + 30*ay*ay*y0*y0 - 15*cx*cx + 30*cx*x0 - 15*cy*cy + 30*cy*y0 - 15*x0*x0 - 15*y0*y0) + 5*a*b*(-6*a*a*ax*ax*ax*ay - 6*a*a*ax*ay*ay*ay + 3*a*a*ax*ay*vp + 12*a*a*ax*ay - 4*a*ax*ax*ax*ax*b + 6*a*ax*ax*ax*ay*cx - 6*a*ax*ax*ax*ay*x0 - 8*a*ax*ax*ay*ay*b + 6*a*ax*ax*ay*ay*cy - 6*a*ax*ax*ay*ay*y0 - 2*a*ax*ax*b*vp + 8*a*ax*ax*b + 6*a*ax*ay*ay*ay*cx - 6*a*ax*ay*ay*ay*x0 - 3*a*ax*ay*cx*vp - 12*a*ax*ay*cx + 3*a*ax*ay*vp*x0 + 12*a*ax*ay*x0 - 4*a*ay*ay*ay*ay*b + 6*a*ay*ay*ay*ay*cy - 6*a*ay*ay*ay*ay*y0 + 4*a*ay*ay*b*vp + 8*a*ay*ay*b - 3*a*ay*ay*cy*vp - 12*a*ay*ay*cy + 3*a*ay*ay*vp*y0 + 12*a*ay*ay*y0 + 2*a*b*beta*rho - 2*a*b*vp - 8*a*b + 3*a*cy*vp + 6*a*cy - 3*a*vp*y0 - 6*a*y0 + 6*ax*ax*ax*ax*b*cx - 6*ax*ax*ax*ax*b*x0 - 6*ax*ax*ax*ay*b*b + 6*ax*ax*ax*ay*b*cy - 6*ax*ax*ax*ay*b*y0 + 6*ax*ax*ay*ay*b*cx - 6*ax*ax*ay*ay*b*x0 + 3*ax*ax*b*cx*vp - 12*ax*ax*b*cx - 3*ax*ax*b*vp*x0 + 12*ax*ax*b*x0 - 6*ax*ay*ay*ay*b*b + 6*ax*ay*ay*ay*b*cy - 6*ax*ay*ay*ay*b*y0 - 3*ax*ay*b*b*vp + 12*ax*ay*b*b + 3*ax*ay*b*cy*vp - 12*ax*ay*b*cy - 3*ax*ay*b*vp*y0 + 12*ax*ay*b*y0 - 3*b*cx*vp + 6*b*cx + 3*b*vp*x0 - 6*b*x0) + b*b*(-15*ax*ax*ax*ax*cx*cx + 30*ax*ax*ax*ax*cx*x0 - 15*ax*ax*ax*ax*x0*x0 + 30*ax*ax*ax*ay*b*cx - 30*ax*ax*ax*ay*b*x0 - 30*ax*ax*ax*ay*cx*cy + 30*ax*ax*ax*ay*cx*y0 + 30*ax*ax*ax*ay*cy*x0 - 30*ax*ax*ax*ay*x0*y0 - 18*ax*ax*ay*ay*b*b + 30*ax*ax*ay*ay*b*cy - 30*ax*ax*ay*ay*b*y0 - 15*ax*ax*ay*ay*cx*cx + 30*ax*ax*ay*ay*cx*x0 - 15*ax*ax*ay*ay*cy*cy + 30*ax*ax*ay*ay*cy*y0 - 15*ax*ax*ay*ay*x0*x0 - 15*ax*ax*ay*ay*y0*y0 + 30*ax*ax*cx*cx - 60*ax*ax*cx*x0 + 30*ax*ax*x0*x0 + 30*ax*ay*ay*ay*b*cx - 30*ax*ay*ay*ay*b*x0 - 30*ax*ay*ay*ay*cx*cy + 30*ax*ay*ay*ay*cx*y0 + 30*ax*ay*ay*ay*cy*x0 - 30*ax*ay*ay*ay*x0*y0 - 60*ax*ay*b*cx + 60*ax*ay*b*x0 + 60*ax*ay*cx*cy - 60*ax*ay*cx*y0 - 60*ax*ay*cy*x0 + 60*ax*ay*x0*y0 - 18*ay*ay*ay*ay*b*b + 30*ay*ay*ay*ay*b*cy - 30*ay*ay*ay*ay*b*y0 - 15*ay*ay*ay*ay*cy*cy + 30*ay*ay*ay*ay*cy*y0 - 15*ay*ay*ay*ay*y0*y0 + 36*ay*ay*b*b - 60*ay*ay*b*cy + 60*ay*ay*b*y0 + 30*ay*ay*cy*cy - 60*ay*ay*cy*y0 + 30*ay*ay*y0*y0 - 18*b*b + 30*b*cy - 30*b*y0 - 15*cx*cx + 30*cx*x0 - 15*cy*cy + 30*cy*y0 - 15*x0*x0 - 15*y0*y0))/(90*a*b)
,
t*(a*a*(-18*a*a*ax*ax*ax*ax - 18*a*a*ax*ax*ay*ay + 36*a*a*ax*ax - 18*a*a + 30*a*ax*ax*ax*ax*cx - 30*a*ax*ax*ax*ax*x0 + 30*a*ax*ax*ax*ay*cy - 30*a*ax*ax*ax*ay*y0 + 30*a*ax*ax*ay*ay*cx - 30*a*ax*ax*ay*ay*x0 - 60*a*ax*ax*cx + 60*a*ax*ax*x0 + 30*a*ax*ay*ay*ay*cy - 30*a*ax*ay*ay*ay*y0 - 60*a*ax*ay*cy + 60*a*ax*ay*y0 + 30*a*cx - 30*a*x0 - 15*ax*ax*ax*ax*cx*cx + 30*ax*ax*ax*ax*cx*x0 - 15*ax*ax*ax*ax*x0*x0 - 30*ax*ax*ax*ay*cx*cy + 30*ax*ax*ax*ay*cx*y0 + 30*ax*ax*ax*ay*cy*x0 - 30*ax*ax*ax*ay*x0*y0 - 15*ax*ax*ay*ay*cx*cx + 30*ax*ax*ay*ay*cx*x0 - 15*ax*ax*ay*ay*cy*cy +30*ax*ax*ay*ay*cy*y0 - 15*ax*ax*ay*ay*x0*x0 - 15*ax*ax*ay*ay*y0*y0 + 30*ax*ax*cx*cx - 60*ax*ax*cx*x0 + 30*ax*ax*x0*x0 - 30*ax*ay*ay*ay*cx*cy + 30*ax*ay*ay*ay*cx*y0 + 30*ax*ay*ay*ay*cy*x0 - 30*ax*ay*ay*ay*x0*y0 + 60*ax*ay*cx*cy - 60*ax*ay*cx*y0 - 60*ax*ay*cy*x0 + 60*ax*ay*x0*y0 - 15*ay*ay*ay*ay*cy*cy + 30*ay*ay*ay*ay*cy*y0 - 15*ay*ay*ay*ay*y0*y0 + 30*ay*ay*cy*cy - 60*ay*ay*cy*y0 + 30*ay*ay*y0*y0 - 15*cx*cx + 30*cx*x0 - 15*cy*cy + 30*cy*y0 - 15*x0*x0 - 15*y0*y0) + 5*a*b*(-6*a*a*ax*ax*ax*ay - 6*a*a*ax*ay*ay*ay - 3*a*a*ax*ay*vp + 12*a*a*ax*ay - 4*a*ax*ax*ax*ax*b + 6*a*ax*ax*ax*ay*cx - 6*a*ax*ax*ax*ay*x0 - 8*a*ax*ax*ay*ay*b + 6*a*ax*ax*ay*ay*cy - 6*a*ax*ax*ay*ay*y0 - 2*a*ax*ax*b*vp + 8*a*ax*ax*b + 6*a*ax*ay*ay*ay*cx - 6*a*ax*ay*ay*ay*x0 + 3*a*ax*ay*cx*vp - 12*a*ax*ay*cx - 3*a*ax*ay*vp*x0 + 12*a*ax*ay*x0 - 4*a*ay*ay*ay*ay*b + 6*a*ay*ay*ay*ay*cy - 6*a*ay*ay*ay*ay*y0 - 2*a*ay*ay*b*vp + 8*a*ay*ay*b + 3*a*ay*ay*cy*vp - 12*a*ay*ay*cy - 3*a*ay*ay*vp*y0 + 12*a*ay*ay*y0 + 2*a*b*beta*rho + 4*a*b*vp - 8*a*b - 3*a*cy*vp + 6*a*cy + 3*a*vp*y0 - 6*a*y0 + 6*ax*ax*ax*ax*b*cx - 6*ax*ax*ax*ax*b*x0 - 6*ax*ax*ax*ay*b*b + 6*ax*ax*ax*ay*b*cy - 6*ax*ax*ax*ay*b*y0 + 6*ax*ax*ay*ay*b*cx - 6*ax*ax*ay*ay*b*x0 + 3*ax*ax*b*cx*vp - 12*ax*ax*b*cx - 3*ax*ax*b*vp*x0 + 12*ax*ax*b*x0 - 6*ax*ay*ay*ay*b*b + 6*ax*ay*ay*ay*b*cy - 6*ax*ay*ay*ay*b*y0 - 3*ax*ay*b*b*vp + 12*ax*ay*b*b + 3*ax*ay*b*cy*vp - 12*ax*ay*b*cy - 3*ax*ay*b*vp*y0 + 12*ax*ay*b*y0 - 3*b*cx*vp + 6*b*cx + 3*b*vp*x0 - 6*b*x0) + b*b*(-15*ax*ax*ax*ax*cx*cx + 30*ax*ax*ax*ax*cx*x0 - 15*ax*ax*ax*ax*x0*x0 + 30*ax*ax*ax*ay*b*cx - 30*ax*ax*ax*ay*b*x0 - 30*ax*ax*ax*ay*cx*cy + 30*ax*ax*ax*ay*cx*y0 + 30*ax*ax*ax*ay*cy*x0 - 30*ax*ax*ax*ay*x0*y0 - 18*ax*ax*ay*ay*b*b + 30*ax*ax*ay*ay*b*cy - 30*ax*ax*ay*ay*b*y0 - 15*ax*ax*ay*ay*cx*cx + 30*ax*ax*ay*ay*cx*x0 - 15*ax*ax*ay*ay*cy*cy + 30*ax*ax*ay*ay*cy*y0 - 15*ax*ax*ay*ay*x0*x0 - 15*ax*ax*ay*ay*y0*y0 + 30*ax*ax*cx*cx - 60*ax*ax*cx*x0 + 30*ax*ax*x0*x0 + 30*ax*ay*ay*ay*b*cx - 30*ax*ay*ay*ay*b*x0 - 30*ax*ay*ay*ay*cx*cy + 30*ax*ay*ay*ay*cx*y0 + 30*ax*ay*ay*ay*cy*x0 - 30*ax*ay*ay*ay*x0*y0 - 60*ax*ay*b*cx + 60*ax*ay*b*x0 + 60*ax*ay*cx*cy - 60*ax*ay*cx*y0 - 60*ax*ay*cy*x0 + 60*ax*ay*x0*y0 - 18*ay*ay*ay*ay*b*b + 30*ay*ay*ay*ay*b*cy - 30*ay*ay*ay*ay*b*y0 - 15*ay*ay*ay*ay*cy*cy + 30*ay*ay*ay*ay*cy*y0 - 15*ay*ay*ay*ay*y0*y0 + 36*ay*ay*b*b - 60*ay*ay*b*cy + 60*ay*ay*b*y0 + 30*ay*ay*cy*cy - 60*ay*ay*cy*y0 + 30*ay*ay*y0*y0 - 18*b*b + 30*b*cy - 30*b*y0 - 15*cx*cx + 30*cx*x0 - 15*cy*cy + 30*cy*y0 - 15*x0*x0 - 15*y0*y0))/(90*a*b)
,
t*(a*a*(-72*a*a*ax*ax*ax*ax - 72*a*a*ax*ax*ay*ay + 144*a*a*ax*ax - 72*a*a + 90*a*ax*ax*ax*ax*cx - 90*a*ax*ax*ax*ax*x0 + 90*a*ax*ax*ax*ay*cy - 90*a*ax*ax*ax*ay*y0 + 90*a*ax*ax*ay*ay*cx - 90*a*ax*ax*ay*ay*x0 - 180*a*ax*ax*cx + 180*a*ax*ax*x0 + 90*a*ax*ay*ay*ay*cy - 90*a*ax*ay*ay*ay*y0 - 180*a*ax*ay*cy + 180*a*ax*ay*y0 + 90*a*cx - 90*a*x0 - 30*ax*ax*ax*ax*cx*cx + 60*ax*ax*ax*ax*cx*x0 - 30*ax*ax*ax*ax*x0*x0 - 60*ax*ax*ax*ay*cx*cy + 60*ax*ax*ax*ay*cx*y0 + 60*ax*ax*ax*ay*cy*x0 - 60*ax*ax*ax*ay*x0*y0 - 30*ax*ax*ay*ay*cx*cx + 60*ax*ax*ay*ay*cx*x0 - 30*ax*ax*ay*ay*cy*cy + 60*ax*ax*ay*ay*cy*y0 - 30*ax*ax*ay*ay*x0*x0 - 30*ax*ax*ay*ay*y0*y0 + 60*ax*ax*cx*cx - 120*ax*ax*cx*x0 + 60*ax*ax*x0*x0 - 60*ax*ay*ay*ay*cx*cy + 60*ax*ay*ay*ay*cx*y0 + 60*ax*ay*ay*ay*cy*x0 - 60*ax*ay*ay*ay*x0*y0 + 120*ax*ay*cx*cy - 120*ax*ay*cx*y0 - 120*ax*ay*cy*x0 + 120*ax*ay*x0*y0 - 30*ay*ay*ay*ay*cy*cy + 60*ay*ay*ay*ay*cy*y0 - 30*ay*ay*ay*ay*y0*y0 + 60*ay*ay*cy*cy -120*ay*ay*cy*y0 + 60*ay*ay*y0*y0 - 30*cx*cx + 60*cx*x0 - 30*cy*cy + 60*cy*y0 - 30*x0*x0 - 30*y0*y0) + 5*a*b*(-18*a*a*ax*ax*ax*ay - 18*a*a*ax*ay*ay*ay - 9*a*a*ax*ay*vp + 36*a*a*ax*ay + 4*a*ax*ax*ax*ax*b + 12*a*ax*ax*ax*ay*cx - 12*a*ax*ax*ax*ay*x0 - 4*a*ax*ax*ay*ay*b + 12*a*ax*ax*ay*ay*cy - 12*a*ax*ax*ay*ay*y0 + 8*a*ax*ax*b*vp - 8*a*ax*ax*b + 12*a*ax*ay*ay*ay*cx - 12*a*ax*ay*ay*ay*x0 + 6*a*ax*ay*cx*vp - 24*a*ax*ay*cx - 6*a*ax*ay*vp*x0 + 24*a*ax*ay*x0 - 8*a*ay*ay*ay*ay*b + 12*a*ay*ay*ay*ay*cy - 12*a*ay*ay*ay*ay*y0 - 4*a*ay*ay*b*vp + 16*a*ay*ay*b + 6*a*ay*ay*cy*vp - 24*a*ay*ay*cy - 6*a*ay*ay*vp*y0 + 24*a*ay*ay*y0 + 4*a*b*beta*rho - 4*a*b*vp - 4*a*b - 6*a*cy*vp + 12*a*cy + 6*a*vp*y0 - 12*a*y0 - 6*ax*ax*ax*ax*b*cx + 6*ax*ax*ax*ax*b*x0 + 6*ax*ax*ax*ay*b*b - 6*ax*ax*ax*ay*b*cy + 6*ax*ax*ax*ay*b*y0 - 6*ax*ax*ay*ay*b*cx + 6*ax*ax*ay*ay*b*x0 - 3*ax*ax*b*cx*vp + 12*ax*ax*b*cx + 3*ax*ax*b*vp*x0 - 12*ax*ax*b*x0 + 6*ax*ay*ay*ay*b*b - 6*ax*ay*ay*ay*b*cy + 6*ax*ay*ay*ay*b*y0 + 3*ax*ay*b*b*vp - 12*ax*ay*b*b - 3*ax*ay*b*cy*vp + 12*ax*ay*b*cy + 3*ax*ay*b*vp*y0 - 12*ax*ay*b*y0 + 3*b*cx*vp - 6*b*cx - 3*b*vp*x0 + 6*b*x0) + b*b*(15*ax*ax*ax*ax*cx*cx - 30*ax*ax*ax*ax*cx*x0 + 15*ax*ax*ax*ax*x0*x0 - 30*ax*ax*ax*ay*b*cx + 30*ax*ax*ax*ay*b*x0 + 30*ax*ax*ax*ay*cx*cy - 30*ax*ax*ax*ay*cx*y0 - 30*ax*ax*ax*ay*cy*x0 + 30*ax*ax*ax*ay*x0*y0 + 18*ax*ax*ay*ay*b*b - 30*ax*ax*ay*ay*b*cy + 30*ax*ax*ay*ay*b*y0 + 15*ax*ax*ay*ay*cx*cx - 30*ax*ax*ay*ay*cx*x0 + 15*ax*ax*ay*ay*cy*cy - 30*ax*ax*ay*ay*cy*y0 + 15*ax*ax*ay*ay*x0*x0 + 15*ax*ax*ay*ay*y0*y0 - 30*ax*ax*cx*cx + 60*ax*ax*cx*x0 - 30*ax*ax*x0*x0 - 30*ax*ay*ay*ay*b*cx + 30*ax*ay*ay*ay*b*x0 + 30*ax*ay*ay*ay*cx*cy - 30*ax*ay*ay*ay*cx*y0 - 30*ax*ay*ay*ay*cy*x0 + 30*ax*ay*ay*ay*x0*y0 + 60*ax*ay*b*cx - 60*ax*ay*b*x0 - 60*ax*ay*cx*cy + 60*ax*ay*cx*y0 + 60*ax*ay*cy*x0 - 60*ax*ay*x0*y0 + 18*ay*ay*ay*ay*b*b - 30*ay*ay*ay*ay*b*cy + 30*ay*ay*ay*ay*b*y0 + 15*ay*ay*ay*ay*cy*cy - 30*ay*ay*ay*ay*cy*y0 + 15*ay*ay*ay*ay*y0*y0 - 36*ay*ay*b*b + 60*ay*ay*b*cy - 60*ay*ay*b*y0 - 30*ay*ay*cy*cy + 60*ay*ay*cy*y0 - 30*ay*ay*y0*y0 + 18*b*b - 30*b*cy + 30*b*y0 + 15*cx*cx - 30*cx*x0 + 15*cy*cy - 30*cy*y0 + 15*x0*x0 + 15*y0*y0))/(90*a*b)
,
t*(a*a*(72*a*a*ax*ax*ax*ax + 72*a*a*ax*ax*ay*ay - 144*a*a*ax*ax + 72*a*a - 90*a*ax*ax*ax*ax*cx + 90*a*ax*ax*ax*ax*x0 - 90*a*ax*ax*ax*ay*cy + 90*a*ax*ax*ax*ay*y0 - 90*a*ax*ax*ay*ay*cx + 90*a*ax*ax*ay*ay*x0 + 180*a*ax*ax*cx - 180*a*ax*ax*x0 - 90*a*ax*ay*ay*ay*cy + 90*a*ax*ay*ay*ay*y0 + 180*a*ax*ay*cy - 180*a*ax*ay*y0 - 90*a*cx + 90*a*x0 + 30*ax*ax*ax*ax*cx*cx - 60*ax*ax*ax*ax*cx*x0 + 30*ax*ax*ax*ax*x0*x0 + 60*ax*ax*ax*ay*cx*cy - 60*ax*ax*ax*ay*cx*y0 - 60*ax*ax*ax*ay*cy*x0 + 60*ax*ax*ax*ay*x0*y0 + 30*ax*ax*ay*ay*cx*cx - 60*ax*ax*ay*ay*cx*x0 + 30*ax*ax*ay*ay*cy*cy - 60*ax*ax*ay*ay*cy*y0 + 30*ax*ax*ay*ay*x0*x0 + 30*ax*ax*ay*ay*y0*y0 - 60*ax*ax*cx*cx + 120*ax*ax*cx*x0 - 60*ax*ax*x0*x0 + 60*ax*ay*ay*ay*cx*cy - 60*ax*ay*ay*ay*cx*y0 - 60*ax*ay*ay*ay*cy*x0 + 60*ax*ay*ay*ay*x0*y0 - 120*ax*ay*cx*cy + 120*ax*ay*cx*y0 + 120*ax*ay*cy*x0 - 120*ax*ay*x0*y0 + 30*ay*ay*ay*ay*cy*cy - 60*ay*ay*ay*ay*cy*y0 + 30*ay*ay*ay*ay*y0*y0 - 60*ay*ay*cy*cy + 120*ay*ay*cy*y0 - 60*ay*ay*y0*y0 + 30*cx*cx - 60*cx*x0 + 30*cy*cy - 60*cy*y0 + 30*x0*x0 + 30*y0*y0) + 5*a*b*(18*a*a*ax*ax*ax*ay + 18*a*a*ax*ay*ay*ay + 9*a*a*ax*ay*vp - 36*a*a*ax*ay + 8*a*ax*ax*ax*ax*b - 12*a*ax*ax*ax*ay*cx + 12*a*ax*ax*ax*ay*x0 + 16*a*ax*ax*ay*ay*b - 12*a*ax*ax*ay*ay*cy + 12*a*ax*ax*ay*ay*y0 + 16*a*ax*ax*b*vp - 16*a*ax*ax*b - 12*a*ax*ay*ay*ay*cx + 12*a*ax*ay*ay*ay*x0 - 6*a*ax*ay*cx*vp + 24*a*ax*ay*cx + 6*a*ax*ay*vp*x0 - 24*a*ax*ay*x0 + 8*a*ay*ay*ay*ay*b - 12*a*ay*ay*ay*ay*cy + 12*a*ay*ay*ay*ay*y0 + 16*a*ay*ay*b*vp - 16*a*ay*ay*b - 6*a*ay*ay*cy*vp + 24*a*ay*ay*cy + 6*a*ay*ay*vp*y0 - 24*a*ay*ay*y0 + 8*a*b*beta*rho - 32*a*b*vp + 16*a*b + 6*a*cy*vp - 12*a*cy - 6*a*vp*y0 + 12*a*y0 - 12*ax*ax*ax*ax*b*cx + 12*ax*ax*ax*ax*b*x0 + 18*ax*ax*ax*ay*b*b - 12*ax*ax*ax*ay*b*cy + 12*ax*ax*ax*ay*b*y0 - 12*ax*ax*ay*ay*b*cx + 12*ax*ax*ay*ay*b*x0 - 6*ax*ax*b*cx*vp + 24*ax*ax*b*cx + 6*ax*ax*b*vp*x0 - 24*ax*ax*b*x0 + 18*ax*ay*ay*ay*b*b - 12*ax*ay*ay*ay*b*cy + 12*ax*ay*ay*ay*b*y0 + 9*ax*ay*b*b*vp - 36*ax*ay*b*b - 6*ax*ay*b*cy*vp + 24*ax*ay*b*cy + 6*ax*ay*b*vp*y0 - 24*ax*ay*b*y0 + 6*b*cx*vp - 12*b*cx - 6*b*vp*x0 + 12*b*x0) + b*b*(30*ax*ax*ax*ax*cx*cx - 60*ax*ax*ax*ax*cx*x0 + 30*ax*ax*ax*ax*x0*x0 - 90*ax*ax*ax*ay*b*cx + 90*ax*ax*ax*ay*b*x0 + 60*ax*ax*ax*ay*cx*cy - 60*ax*ax*ax*ay*cx*y0 - 60*ax*ax*ax*ay*cy*x0 + 60*ax*ax*ax*ay*x0*y0 + 72*ax*ax*ay*ay*b*b - 90*ax*ax*ay*ay*b*cy + 90*ax*ax*ay*ay*b*y0 + 30*ax*ax*ay*ay*cx*cx - 60*ax*ax*ay*ay*cx*x0 + 30*ax*ax*ay*ay*cy*cy - 60*ax*ax*ay*ay*cy*y0 + 30*ax*ax*ay*ay*x0*x0 + 30*ax*ax*ay*ay*y0*y0 - 60*ax*ax*cx*cx + 120*ax*ax*cx*x0 - 60*ax*ax*x0*x0 - 90*ax*ay*ay*ay*b*cx + 90*ax*ay*ay*ay*b*x0 + 60*ax*ay*ay*ay*cx*cy - 60*ax*ay*ay*ay*cx*y0 - 60*ax*ay*ay*ay*cy*x0 + 60*ax*ay*ay*ay*x0*y0 + 180*ax*ay*b*cx - 180*ax*ay*b*x0 - 120*ax*ay*cx*cy + 120*ax*ay*cx*y0 + 120*ax*ay*cy*x0 - 120*ax*ay*x0*y0 + 72*ay*ay*ay*ay*b*b - 90*ay*ay*ay*ay*b*cy + 90*ay*ay*ay*ay*b*y0 + 30*ay*ay*ay*ay*cy*cy - 60*ay*ay*ay*ay*cy*y0 + 30*ay*ay*ay*ay*y0*y0 - 144*ay*ay*b*b + 180*ay*ay*b*cy - 180*ay*ay*b*y0 - 60*ay*ay*cy*cy + 120*ay*ay*cy*y0 - 60*ay*ay*y0*y0 + 72*b*b - 90*b*cy + 90*b*y0 + 30*cx*cx - 60*cx*x0 + 30*cy*cy - 60*cy*y0 + 30*x0*x0 + 30*y0*y0))/(90*a*b)
,
t*(a*a*(18*a*a*ax*ax*ax*ax + 18*a*a*ax*ax*ay*ay - 36*a*a*ax*ax + 18*a*a - 30*a*ax*ax*ax*ax*cx + 30*a*ax*ax*ax*ax*x0 - 30*a*ax*ax*ax*ay*cy + 30*a*ax*ax*ax*ay*y0 - 30*a*ax*ax*ay*ay*cx + 30*a*ax*ax*ay*ay*x0 + 60*a*ax*ax*cx - 60*a*ax*ax*x0 - 30*a*ax*ay*ay*ay*cy + 30*a*ax*ay*ay*ay*y0 + 60*a*ax*ay*cy - 60*a*ax*ay*y0 - 30*a*cx + 30*a*x0 + 15*ax*ax*ax*ax*cx*cx - 30*ax*ax*ax*ax*cx*x0 + 15*ax*ax*ax*ax*x0*x0 + 30*ax*ax*ax*ay*cx*cy - 30*ax*ax*ax*ay*cx*y0 - 30*ax*ax*ax*ay*cy*x0 + 30*ax*ax*ax*ay*x0*y0 + 15*ax*ax*ay*ay*cx*cx - 30*ax*ax*ay*ay*cx*x0 + 15*ax*ax*ay*ay*cy*cy - 30*ax*ax*ay*ay*cy*y0 + 15*ax*ax*ay*ay*x0*x0 + 15*ax*ax*ay*ay*y0*y0 - 30*ax*ax*cx*cx + 60*ax*ax*cx*x0 - 30*ax*ax*x0*x0 + 30*ax*ay*ay*ay*cx*cy - 30*ax*ay*ay*ay*cx*y0 - 30*ax*ay*ay*ay*cy*x0 +30*ax*ay*ay*ay*x0*y0 - 60*ax*ay*cx*cy + 60*ax*ay*cx*y0 + 60*ax*ay*cy*x0 - 60*ax*ay*x0*y0 + 15*ay*ay*ay*ay*cy*cy - 30*ay*ay*ay*ay*cy*y0 + 15*ay*ay*ay*ay*y0*y0 - 30*ay*ay*cy*cy + 60*ay*ay*cy*y0 - 30*ay*ay*y0*y0 + 15*cx*cx - 30*cx*x0 + 15*cy*cy - 30*cy*y0 + 15*x0*x0 + 15*y0*y0) + 5*a*b*(6*a*a*ax*ax*ax*ay + 6*a*a*ax*ay*ay*ay + 3*a*a*ax*ay*vp - 12*a*a*ax*ay - 8*a*ax*ax*ax*ax*b - 6*a*ax*ax*ax*ay*cx + 6*a*ax*ax*ax*ay*x0 - 4*a*ax*ax*ay*ay*b - 6*a*ax*ax*ay*ay*cy + 6*a*ax*ax*ay*ay*y0 - 4*a*ax*ax*b*vp + 16*a*ax*ax*b - 6*a*ax*ay*ay*ay*cx + 6*a*ax*ay*ay*ay*x0 - 3*a*ax*ay*cx*vp + 12*a*ax*ay*cx + 3*a*ax*ay*vp*x0 - 12*a*ax*ay*x0 + 4*a*ay*ay*ay*ay*b - 6*a*ay*ay*ay*ay*cy + 6*a*ay*ay*ay*ay*y0 + 8*a*ay*ay*b*vp - 8*a*ay*ay*b - 3*a*ay*ay*cy*vp + 12*a*ay*ay*cy + 3*a*ay*ay*vp*y0 - 12*a*ay*ay*y0 + 4*a*b*beta*rho - 4*a*b*vp - 4*a*b + 3*a*cy*vp - 6*a*cy - 3*a*vp*y0 + 6*a*y0 + 12*ax*ax*ax*ax*b*cx - 12*ax*ax*ax*ax*b*x0 - 18*ax*ax*ax*ay*b*b + 12*ax*ax*ax*ay*b*cy - 12*ax*ax*ax*ay*b*y0 + 12*ax*ax*ay*ay*b*cx - 12*ax*ax*ay*ay*b*x0 + 6*ax*ax*b*cx*vp - 24*ax*ax*b*cx - 6*ax*ax*b*vp*x0 + 24*ax*ax*b*x0 - 18*ax*ay*ay*ay*b*b + 12*ax*ay*ay*ay*b*cy - 12*ax*ay*ay*ay*b*y0 - 9*ax*ay*b*b*vp + 36*ax*ay*b*b + 6*ax*ay*b*cy*vp - 24*ax*ay*b*cy - 6*ax*ay*b*vp*y0 + 24*ax*ay*b*y0 - 6*b*cx*vp + 12*b*cx + 6*b*vp*x0 - 12*b*x0) + b*b*(-30*ax*ax*ax*ax*cx*cx + 60*ax*ax*ax*ax*cx*x0 - 30*ax*ax*ax*ax*x0*x0 + 90*ax*ax*ax*ay*b*cx - 90*ax*ax*ax*ay*b*x0 - 60*ax*ax*ax*ay*cx*cy + 60*ax*ax*ax*ay*cx*y0 + 60*ax*ax*ax*ay*cy*x0 - 60*ax*ax*ax*ay*x0*y0 -72*ax*ax*ay*ay*b*b + 90*ax*ax*ay*ay*b*cy - 90*ax*ax*ay*ay*b*y0 - 30*ax*ax*ay*ay*cx*cx + 60*ax*ax*ay*ay*cx*x0 - 30*ax*ax*ay*ay*cy*cy + 60*ax*ax*ay*ay*cy*y0 - 30*ax*ax*ay*ay*x0*x0 - 30*ax*ax*ay*ay*y0*y0 + 60*ax*ax*cx*cx - 120*ax*ax*cx*x0 + 60*ax*ax*x0*x0 + 90*ax*ay*ay*ay*b*cx - 90*ax*ay*ay*ay*b*x0 - 60*ax*ay*ay*ay*cx*cy + 60*ax*ay*ay*ay*cx*y0 + 60*ax*ay*ay*ay*cy*x0 - 60*ax*ay*ay*ay*x0*y0 - 180*ax*ay*b*cx + 180*ax*ay*b*x0 + 120*ax*ay*cx*cy - 120*ax*ay*cx*y0 - 120*ax*ay*cy*x0 + 120*ax*ay*x0*y0 - 72*ay*ay*ay*ay*b*b + 90*ay*ay*ay*ay*b*cy - 90*ay*ay*ay*ay*b*y0 - 30*ay*ay*ay*ay*cy*cy + 60*ay*ay*ay*ay*cy*y0 - 30*ay*ay*ay*ay*y0*y0 + 144*ay*ay*b*b - 180*ay*ay*b*cy + 180*ay*ay*b*y0 + 60*ay*ay*cy*cy - 120*ay*ay*cy*y0 + 60*ay*ay*y0*y0 - 72*b*b + 90*b*cy -90*b*y0 - 30*cx*cx + 60*cx*x0 - 30*cy*cy + 60*cy*y0 - 30*x0*x0 - 30*y0*y0))/(90*a*b)
,
t*(a*a*(-12*a*a*ax*ax*ax*ax - 12*a*a*ax*ax*ay*ay + 24*a*a*ax*ax - 12*a*a + 30*a*ax*ax*ax*ax*cx - 30*a*ax*ax*ax*ax*x0 + 30*a*ax*ax*ax*ay*cy - 30*a*ax*ax*ax*ay*y0 + 30*a*ax*ax*ay*ay*cx - 30*a*ax*ax*ay*ay*x0 - 60*a*ax*ax*cx + 60*a*ax*ax*x0 + 30*a*ax*ay*ay*ay*cy - 30*a*ax*ay*ay*ay*y0 - 60*a*ax*ay*cy + 60*a*ax*ay*y0 + 30*a*cx - 30*a*x0 - 30*ax*ax*ax*ax*cx*cx + 60*ax*ax*ax*ax*cx*x0 - 30*ax*ax*ax*ax*x0*x0 - 60*ax*ax*ax*ay*cx*cy + 60*ax*ax*ax*ay*cx*y0 + 60*ax*ax*ax*ay*cy*x0 - 60*ax*ax*ax*ay*x0*y0 - 30*ax*ax*ay*ay*cx*cx + 60*ax*ax*ay*ay*cx*x0 - 30*ax*ax*ay*ay*cy*cy +60*ax*ax*ay*ay*cy*y0 - 30*ax*ax*ay*ay*x0*x0 - 30*ax*ax*ay*ay*y0*y0 + 60*ax*ax*cx*cx - 120*ax*ax*cx*x0 + 60*ax*ax*x0*x0 - 60*ax*ay*ay*ay*cx*cy + 60*ax*ay*ay*ay*cx*y0 + 60*ax*ay*ay*ay*cy*x0- 60*ax*ay*ay*ay*x0*y0 + 120*ax*ay*cx*cy - 120*ax*ay*cx*y0 - 120*ax*ay*cy*x0 + 120*ax*ay*x0*y0 - 30*ay*ay*ay*ay*cy*cy + 60*ay*ay*ay*ay*cy*y0 - 30*ay*ay*ay*ay*y0*y0 + 60*ay*ay*cy*cy - 120*ay*ay*cy*y0 + 60*ay*ay*y0*y0 - 30*cx*cx + 60*cx*x0 - 30*cy*cy + 60*cy*y0 - 30*x0*x0 - 30*y0*y0) + 5*a*b*(-6*a*a*ax*ax*ax*ay - 6*a*a*ax*ay*ay*ay - 3*a*a*ax*ay*vp + 12*a*a*ax*ay + 4*a*ax*ax*ax*ax*b + 12*a*ax*ax*ax*ay*cx - 12*a*ax*ax*ax*ay*x0 - 4*a*ax*ax*ay*ay*b + 12*a*ax*ax*ay*ay*cy - 12*a*ax*ax*ay*ay*y0 + 2*a*ax*ax*b*vp - 8*a*ax*ax*b + 12*a*ax*ay*ay*ay*cx - 12*a*ax*ay*ay*ay*x0 + 6*a*ax*ay*cx*vp - 24*a*ax*ay*cx - 6*a*ax*ay*vp*x0 + 24*a*ax*ay*x0 - 8*a*ay*ay*ay*ay*b + 12*a*ay*ay*ay*ay*cy - 12*a*ay*ay*ay*ay*y0 - 4*a*ay*ay*b*vp + 16*a*ay*ay*b + 6*a*ay*ay*cy*vp - 24*a*ay*ay*cy - 6*a*ay*ay*vp*y0 + 24*a*ay*ay*y0 + 4*a*b*beta*rho + 2*a*b*vp - 4*a*b - 6*a*cy*vp + 12*a*cy + 6*a*vp*y0 - 12*a*y0 - 6*ax*ax*ax*ax*b*cx + 6*ax*ax*ax*ax*b*x0 + 6*ax*ax*ax*ay*b*b- 6*ax*ax*ax*ay*b*cy + 6*ax*ax*ax*ay*b*y0 - 6*ax*ax*ay*ay*b*cx + 6*ax*ax*ay*ay*b*x0 + 3*ax*ax*b*cx*vp + 12*ax*ax*b*cx - 3*ax*ax*b*vp*x0 - 12*ax*ax*b*x0 + 6*ax*ay*ay*ay*b*b - 6*ax*ay*ay*ay*b*cy + 6*ax*ay*ay*ay*b*y0 - 3*ax*ay*b*b*vp - 12*ax*ay*b*b + 3*ax*ay*b*cy*vp + 12*ax*ay*b*cy - 3*ax*ay*b*vp*y0 - 12*ax*ay*b*y0 - 3*b*cx*vp - 6*b*cx + 3*b*vp*x0 + 6*b*x0) + b*b*(15*ax*ax*ax*ax*cx*cx - 30*ax*ax*ax*ax*cx*x0 + 15*ax*ax*ax*ax*x0*x0 - 30*ax*ax*ax*ay*b*cx + 30*ax*ax*ax*ay*b*x0 + 30*ax*ax*ax*ay*cx*cy - 30*ax*ax*ax*ay*cx*y0 - 30*ax*ax*ax*ay*cy*x0 + 30*ax*ax*ax*ay*x0*y0 + 18*ax*ax*ay*ay*b*b - 30*ax*ax*ay*ay*b*cy + 30*ax*ax*ay*ay*b*y0 + 15*ax*ax*ay*ay*cx*cx - 30*ax*ax*ay*ay*cx*x0 + 15*ax*ax*ay*ay*cy*cy - 30*ax*ax*ay*ay*cy*y0 + 15*ax*ax*ay*ay*x0*x0 + 15*ax*ax*ay*ay*y0*y0 - 30*ax*ax*cx*cx + 60*ax*ax*cx*x0 - 30*ax*ax*x0*x0 - 30*ax*ay*ay*ay*b*cx + 30*ax*ay*ay*ay*b*x0 + 30*ax*ay*ay*ay*cx*cy - 30*ax*ay*ay*ay*cx*y0 - 30*ax*ay*ay*ay*cy*x0 + 30*ax*ay*ay*ay*x0*y0 + 60*ax*ay*b*cx - 60*ax*ay*b*x0 - 60*ax*ay*cx*cy + 60*ax*ay*cx*y0 + 60*ax*ay*cy*x0 - 60*ax*ay*x0*y0 + 18*ay*ay*ay*ay*b*b - 30*ay*ay*ay*ay*b*cy + 30*ay*ay*ay*ay*b*y0 + 15*ay*ay*ay*ay*cy*cy - 30*ay*ay*ay*ay*cy*y0 + 15*ay*ay*ay*ay*y0*y0 - 36*ay*ay*b*b + 60*ay*ay*b*cy - 60*ay*ay*b*y0 - 30*ay*ay*cy*cy + 60*ay*ay*cy*y0 - 30*ay*ay*y0*y0 + 18*b*b - 30*b*cy + 30*b*y0 + 15*cx*cx - 30*cx*x0 + 15*cy*cy - 30*cy*y0 + 15*x0*x0 + 15*y0*y0))/(90*a*b)
,
t*(a*a*(-18*a*a*ax*ax*ax*ax - 18*a*a*ax*ax*ay*ay + 36*a*a*ax*ax - 18*a*a + 30*a*ax*ax*ax*ax*cx - 30*a*ax*ax*ax*ax*x0 + 30*a*ax*ax*ax*ay*cy - 30*a*ax*ax*ax*ay*y0 + 30*a*ax*ax*ay*ay*cx - 30*a*ax*ax*ay*ay*x0 - 60*a*ax*ax*cx + 60*a*ax*ax*x0 + 30*a*ax*ay*ay*ay*cy - 30*a*ax*ay*ay*ay*y0 - 60*a*ax*ay*cy + 60*a*ax*ay*y0 + 30*a*cx - 30*a*x0 - 15*ax*ax*ax*ax*cx*cx + 30*ax*ax*ax*ax*cx*x0 - 15*ax*ax*ax*ax*x0*x0 - 30*ax*ax*ax*ay*cx*cy + 30*ax*ax*ax*ay*cx*y0 + 30*ax*ax*ax*ay*cy*x0 - 30*ax*ax*ax*ay*x0*y0 - 15*ax*ax*ay*ay*cx*cx + 30*ax*ax*ay*ay*cx*x0 - 15*ax*ax*ay*ay*cy*cy +30*ax*ax*ay*ay*cy*y0 - 15*ax*ax*ay*ay*x0*x0 - 15*ax*ax*ay*ay*y0*y0 + 30*ax*ax*cx*cx - 60*ax*ax*cx*x0 + 30*ax*ax*x0*x0 - 30*ax*ay*ay*ay*cx*cy + 30*ax*ay*ay*ay*cx*y0 + 30*ax*ay*ay*ay*cy*x0 - 30*ax*ay*ay*ay*x0*y0 + 60*ax*ay*cx*cy - 60*ax*ay*cx*y0 - 60*ax*ay*cy*x0 + 60*ax*ay*x0*y0 - 15*ay*ay*ay*ay*cy*cy + 30*ay*ay*ay*ay*cy*y0 - 15*ay*ay*ay*ay*y0*y0 + 30*ay*ay*cy*cy - 60*ay*ay*cy*y0 + 30*ay*ay*y0*y0 - 15*cx*cx + 30*cx*x0 - 15*cy*cy + 30*cy*y0 - 15*x0*x0 - 15*y0*y0) + 5*a*b*(-6*a*a*ax*ax*ax*ay - 6*a*a*ax*ay*ay*ay - 3*a*a*ax*ay*vp + 12*a*a*ax*ay - 4*a*ax*ax*ax*ax*b + 6*a*ax*ax*ax*ay*cx - 6*a*ax*ax*ax*ay*x0 - 8*a*ax*ax*ay*ay*b + 6*a*ax*ax*ay*ay*cy - 6*a*ax*ax*ay*ay*y0 + 4*a*ax*ax*b*vp + 8*a*ax*ax*b + 6*a*ax*ay*ay*ay*cx - 6*a*ax*ay*ay*ay*x0 + 3*a*ax*ay*cx*vp - 12*a*ax*ay*cx - 3*a*ax*ay*vp*x0 + 12*a*ax*ay*x0 - 4*a*ay*ay*ay*ay*b + 6*a*ay*ay*ay*ay*cy - 6*a*ay*ay*ay*ay*y0 - 2*a*ay*ay*b*vp + 8*a*ay*ay*b + 3*a*ay*ay*cy*vp - 12*a*ay*ay*cy - 3*a*ay*ay*vp*y0 + 12*a*ay*ay*y0 + 2*a*b*beta*rho - 2*a*b*vp - 8*a*b - 3*a*cy*vp + 6*a*cy + 3*a*vp*y0 - 6*a*y0 + 6*ax*ax*ax*ax*b*cx - 6*ax*ax*ax*ax*b*x0 - 6*ax*ax*ax*ay*b*b + 6*ax*ax*ax*ay*b*cy - 6*ax*ax*ax*ay*b*y0 + 6*ax*ax*ay*ay*b*cx - 6*ax*ax*ay*ay*b*x0 - 3*ax*ax*b*cx*vp - 12*ax*ax*b*cx + 3*ax*ax*b*vp*x0 + 12*ax*ax*b*x0 - 6*ax*ay*ay*ay*b*b + 6*ax*ay*ay*ay*b*cy - 6*ax*ay*ay*ay*b*y0 + 3*ax*ay*b*b*vp + 12*ax*ay*b*b - 3*ax*ay*b*cy*vp - 12*ax*ay*b*cy + 3*ax*ay*b*vp*y0 + 12*ax*ay*b*y0 + 3*b*cx*vp + 6*b*cx - 3*b*vp*x0 - 6*b*x0) + b*b*(-15*ax*ax*ax*ax*cx*cx + 30*ax*ax*ax*ax*cx*x0 - 15*ax*ax*ax*ax*x0*x0 + 30*ax*ax*ax*ay*b*cx - 30*ax*ax*ax*ay*b*x0 - 30*ax*ax*ax*ay*cx*cy + 30*ax*ax*ax*ay*cx*y0 + 30*ax*ax*ax*ay*cy*x0 - 30*ax*ax*ax*ay*x0*y0 - 18*ax*ax*ay*ay*b*b + 30*ax*ax*ay*ay*b*cy - 30*ax*ax*ay*ay*b*y0 - 15*ax*ax*ay*ay*cx*cx + 30*ax*ax*ay*ay*cx*x0 - 15*ax*ax*ay*ay*cy*cy + 30*ax*ax*ay*ay*cy*y0 - 15*ax*ax*ay*ay*x0*x0 - 15*ax*ax*ay*ay*y0*y0 + 30*ax*ax*cx*cx - 60*ax*ax*cx*x0 + 30*ax*ax*x0*x0 + 30*ax*ay*ay*ay*b*cx - 30*ax*ay*ay*ay*b*x0 - 30*ax*ay*ay*ay*cx*cy + 30*ax*ay*ay*ay*cx*y0 + 30*ax*ay*ay*ay*cy*x0 - 30*ax*ay*ay*ay*x0*y0 - 60*ax*ay*b*cx + 60*ax*ay*b*x0 + 60*ax*ay*cx*cy - 60*ax*ay*cx*y0 - 60*ax*ay*cy*x0 + 60*ax*ay*x0*y0 - 18*ay*ay*ay*ay*b*b + 30*ay*ay*ay*ay*b*cy - 30*ay*ay*ay*ay*b*y0 - 15*ay*ay*ay*ay*cy*cy + 30*ay*ay*ay*ay*cy*y0 - 15*ay*ay*ay*ay*y0*y0 + 36*ay*ay*b*b - 60*ay*ay*b*cy + 60*ay*ay*b*y0 + 30*ay*ay*cy*cy - 60*ay*ay*cy*y0 + 30*ay*ay*y0*y0 - 18*b*b + 30*b*cy - 30*b*y0 - 15*cx*cx + 30*cx*x0 - 15*cy*cy + 30*cy*y0 - 15*x0*x0 - 15*y0*y0))/(90*a*b)
,
t*(a*a*(18*a*a*ax*ax*ax*ax + 18*a*a*ax*ax*ay*ay - 36*a*a*ax*ax + 18*a*a - 30*a*ax*ax*ax*ax*cx + 30*a*ax*ax*ax*ax*x0 - 30*a*ax*ax*ax*ay*cy + 30*a*ax*ax*ax*ay*y0 - 30*a*ax*ax*ay*ay*cx + 30*a*ax*ax*ay*ay*x0 + 60*a*ax*ax*cx - 60*a*ax*ax*x0 - 30*a*ax*ay*ay*ay*cy + 30*a*ax*ay*ay*ay*y0 + 60*a*ax*ay*cy - 60*a*ax*ay*y0 - 30*a*cx + 30*a*x0 + 15*ax*ax*ax*ax*cx*cx - 30*ax*ax*ax*ax*cx*x0 + 15*ax*ax*ax*ax*x0*x0 + 30*ax*ax*ax*ay*cx*cy - 30*ax*ax*ax*ay*cx*y0 - 30*ax*ax*ax*ay*cy*x0 + 30*ax*ax*ax*ay*x0*y0 + 15*ax*ax*ay*ay*cx*cx - 30*ax*ax*ay*ay*cx*x0 + 15*ax*ax*ay*ay*cy*cy - 30*ax*ax*ay*ay*cy*y0 + 15*ax*ax*ay*ay*x0*x0 + 15*ax*ax*ay*ay*y0*y0 - 30*ax*ax*cx*cx + 60*ax*ax*cx*x0 - 30*ax*ax*x0*x0 + 30*ax*ay*ay*ay*cx*cy - 30*ax*ay*ay*ay*cx*y0 - 30*ax*ay*ay*ay*cy*x0 +30*ax*ay*ay*ay*x0*y0 - 60*ax*ay*cx*cy + 60*ax*ay*cx*y0 + 60*ax*ay*cy*x0 - 60*ax*ay*x0*y0 + 15*ay*ay*ay*ay*cy*cy - 30*ay*ay*ay*ay*cy*y0 + 15*ay*ay*ay*ay*y0*y0 - 30*ay*ay*cy*cy + 60*ay*ay*cy*y0 - 30*ay*ay*y0*y0 + 15*cx*cx - 30*cx*x0 + 15*cy*cy - 30*cy*y0 + 15*x0*x0 + 15*y0*y0) + 5*a*b*(6*a*a*ax*ax*ax*ay + 6*a*a*ax*ay*ay*ay + 3*a*a*ax*ay*vp - 12*a*a*ax*ay - 8*a*ax*ax*ax*ax*b - 6*a*ax*ax*ax*ay*cx + 6*a*ax*ax*ax*ay*x0 - 4*a*ax*ax*ay*ay*b - 6*a*ax*ax*ay*ay*cy + 6*a*ax*ax*ay*ay*y0 + 8*a*ax*ax*b*vp + 16*a*ax*ax*b - 6*a*ax*ay*ay*ay*cx + 6*a*ax*ay*ay*ay*x0 - 3*a*ax*ay*cx*vp + 12*a*ax*ay*cx + 3*a*ax*ay*vp*x0 - 12*a*ax*ay*x0 + 4*a*ay*ay*ay*ay*b - 6*a*ay*ay*ay*ay*cy + 6*a*ay*ay*ay*ay*y0 + 8*a*ay*ay*b*vp - 8*a*ay*ay*b - 3*a*ay*ay*cy*vp + 12*a*ay*ay*cy + 3*a*ay*ay*vp*y0 - 12*a*ay*ay*y0 + 4*a*b*beta*rho - 16*a*b*vp - 4*a*b + 3*a*cy*vp - 6*a*cy - 3*a*vp*y0 + 6*a*y0 + 12*ax*ax*ax*ax*b*cx - 12*ax*ax*ax*ax*b*x0 - 18*ax*ax*ax*ay*b*b + 12*ax*ax*ax*ay*b*cy - 12*ax*ax*ax*ay*b*y0 + 12*ax*ax*ay*ay*b*cx - 12*ax*ax*ay*ay*b*x0 - 6*ax*ax*b*cx*vp - 24*ax*ax*b*cx + 6*ax*ax*b*vp*x0 + 24*ax*ax*b*x0 - 18*ax*ay*ay*ay*b*b + 12*ax*ay*ay*ay*b*cy - 12*ax*ay*ay*ay*b*y0 + 9*ax*ay*b*b*vp + 36*ax*ay*b*b - 6*ax*ay*b*cy*vp - 24*ax*ay*b*cy + 6*ax*ay*b*vp*y0 + 24*ax*ay*b*y0 + 6*b*cx*vp + 12*b*cx - 6*b*vp*x0 - 12*b*x0) + b*b*(-30*ax*ax*ax*ax*cx*cx + 60*ax*ax*ax*ax*cx*x0 - 30*ax*ax*ax*ax*x0*x0 + 90*ax*ax*ax*ay*b*cx - 90*ax*ax*ax*ay*b*x0 - 60*ax*ax*ax*ay*cx*cy + 60*ax*ax*ax*ay*cx*y0 + 60*ax*ax*ax*ay*cy*x0 - 60*ax*ax*ax*ay*x0*y0 - 72*ax*ax*ay*ay*b*b + 90*ax*ax*ay*ay*b*cy - 90*ax*ax*ay*ay*b*y0 - 30*ax*ax*ay*ay*cx*cx + 60*ax*ax*ay*ay*cx*x0 - 30*ax*ax*ay*ay*cy*cy + 60*ax*ax*ay*ay*cy*y0 - 30*ax*ax*ay*ay*x0*x0 - 30*ax*ax*ay*ay*y0*y0 + 60*ax*ax*cx*cx - 120*ax*ax*cx*x0 + 60*ax*ax*x0*x0 + 90*ax*ay*ay*ay*b*cx - 90*ax*ay*ay*ay*b*x0 - 60*ax*ay*ay*ay*cx*cy + 60*ax*ay*ay*ay*cx*y0 + 60*ax*ay*ay*ay*cy*x0 - 60*ax*ay*ay*ay*x0*y0 - 180*ax*ay*b*cx + 180*ax*ay*b*x0 + 120*ax*ay*cx*cy - 120*ax*ay*cx*y0 - 120*ax*ay*cy*x0 + 120*ax*ay*x0*y0 - 72*ay*ay*ay*ay*b*b + 90*ay*ay*ay*ay*b*cy - 90*ay*ay*ay*ay*b*y0 - 30*ay*ay*ay*ay*cy*cy + 60*ay*ay*ay*ay*cy*y0 - 30*ay*ay*ay*ay*y0*y0 + 144*ay*ay*b*b - 180*ay*ay*b*cy + 180*ay*ay*b*y0 + 60*ay*ay*cy*cy - 120*ay*ay*cy*y0 + 60*ay*ay*y0*y0 - 72*b*b + 90*b*cy - 90*b*y0 - 30*cx*cx + 60*cx*x0 - 30*cy*cy + 60*cy*y0 - 30*x0*x0 - 30*y0*y0))/(90*a*b)
,
t*(a*a*(12*a*a*ax*ax*ax*ax + 12*a*a*ax*ax*ay*ay - 24*a*a*ax*ax + 12*a*a - 30*a*ax*ax*ax*ax*cx + 30*a*ax*ax*ax*ax*x0 - 30*a*ax*ax*ax*ay*cy + 30*a*ax*ax*ax*ay*y0 - 30*a*ax*ax*ay*ay*cx + 30*a*ax*ax*ay*ay*x0 + 60*a*ax*ax*cx - 60*a*ax*ax*x0 - 30*a*ax*ay*ay*ay*cy + 30*a*ax*ay*ay*ay*y0 + 60*a*ax*ay*cy - 60*a*ax*ay*y0 - 30*a*cx + 30*a*x0 + 30*ax*ax*ax*ax*cx*cx - 60*ax*ax*ax*ax*cx*x0 + 30*ax*ax*ax*ax*x0*x0 + 60*ax*ax*ax*ay*cx*cy - 60*ax*ax*ax*ay*cx*y0 - 60*ax*ax*ax*ay*cy*x0 + 60*ax*ax*ax*ay*x0*y0 + 30*ax*ax*ay*ay*cx*cx - 60*ax*ax*ay*ay*cx*x0 + 30*ax*ax*ay*ay*cy*cy - 60*ax*ax*ay*ay*cy*y0 + 30*ax*ax*ay*ay*x0*x0 + 30*ax*ax*ay*ay*y0*y0 - 60*ax*ax*cx*cx + 120*ax*ax*cx*x0 - 60*ax*ax*x0*x0 + 60*ax*ay*ay*ay*cx*cy - 60*ax*ay*ay*ay*cx*y0 - 60*ax*ay*ay*ay*cy*x0 + 60*ax*ay*ay*ay*x0*y0 - 120*ax*ay*cx*cy + 120*ax*ay*cx*y0 + 120*ax*ay*cy*x0 - 120*ax*ay*x0*y0 + 30*ay*ay*ay*ay*cy*cy - 60*ay*ay*ay*ay*cy*y0 + 30*ay*ay*ay*ay*y0*y0 - 60*ay*ay*cy*cy + 120*ay*ay*cy*y0 - 60*ay*ay*y0*y0 + 30*cx*cx - 60*cx*x0 + 30*cy*cy - 60*cy*y0 + 30*x0*x0 + 30*y0*y0) + 5*a*b*(6*a*a*ax*ax*ax*ay + 6*a*a*ax*ay*ay*ay + 3*a*a*ax*ay*vp - 12*a*a*ax*ay + 8*a*ax*ax*ax*ax*b - 12*a*ax*ax*ax*ay*cx + 12*a*ax*ax*ax*ay*x0 + 16*a*ax*ax*ay*ay*b - 12*a*ax*ax*ay*ay*cy + 12*a*ax*ax*ay*ay*y0 + 4*a*ax*ax*b*vp - 16*a*ax*ax*b - 12*a*ax*ay*ay*ay*cx + 12*a*ax*ay*ay*ay*x0 - 6*a*ax*ay*cx*vp + 24*a*ax*ay*cx + 6*a*ax*ay*vp*x0 - 24*a*ax*ay*x0 + 8*a*ay*ay*ay*ay*b - 12*a*ay*ay*ay*ay*cy + 12*a*ay*ay*ay*ay*y0 + 16*a*ay*ay*b*vp - 16*a*ay*ay*b - 6*a*ay*ay*cy*vp + 24*a*ay*ay*cy + 6*a*ay*ay*vp*y0 - 24*a*ay*ay*y0 + 8*a*b*beta*rho - 20*a*b*vp + 16*a*b + 6*a*cy*vp - 12*a*cy - 6*a*vp*y0 + 12*a*y0 - 12*ax*ax*ax*ax*b*cx + 12*ax*ax*ax*ax*b*x0 + 18*ax*ax*ax*ay*b*b - 12*ax*ax*ax*ay*b*cy + 12*ax*ax*ax*ay*b*y0 - 12*ax*ax*ay*ay*b*cx + 12*ax*ax*ay*ay*b*x0 + 6*ax*ax*b*cx*vp + 24*ax*ax*b*cx - 6*ax*ax*b*vp*x0 - 24*ax*ax*b*x0 + 18*ax*ay*ay*ay*b*b - 12*ax*ay*ay*ay*b*cy + 12*ax*ay*ay*ay*b*y0 - 9*ax*ay*b*b*vp - 36*ax*ay*b*b + 6*ax*ay*b*cy*vp + 24*ax*ay*b*cy - 6*ax*ay*b*vp*y0 - 24*ax*ay*b*y0 - 6*b*cx*vp - 12*b*cx + 6*b*vp*x0 + 12*b*x0) + b*b*(30*ax*ax*ax*ax*cx*cx - 60*ax*ax*ax*ax*cx*x0 + 30*ax*ax*ax*ax*x0*x0 - 90*ax*ax*ax*ay*b*cx + 90*ax*ax*ax*ay*b*x0 + 60*ax*ax*ax*ay*cx*cy - 60*ax*ax*ax*ay*cx*y0 - 60*ax*ax*ax*ay*cy*x0 + 60*ax*ax*ax*ay*x0*y0 + 72*ax*ax*ay*ay*b*b - 90*ax*ax*ay*ay*b*cy + 90*ax*ax*ay*ay*b*y0 + 30*ax*ax*ay*ay*cx*cx - 60*ax*ax*ay*ay*cx*x0 + 30*ax*ax*ay*ay*cy*cy - 60*ax*ax*ay*ay*cy*y0 + 30*ax*ax*ay*ay*x0*x0 + 30*ax*ax*ay*ay*y0*y0 - 60*ax*ax*cx*cx + 120*ax*ax*cx*x0 - 60*ax*ax*x0*x0 - 90*ax*ay*ay*ay*b*cx + 90*ax*ay*ay*ay*b*x0 + 60*ax*ay*ay*ay*cx*cy - 60*ax*ay*ay*ay*cx*y0 - 60*ax*ay*ay*ay*cy*x0 + 60*ax*ay*ay*ay*x0*y0 + 180*ax*ay*b*cx - 180*ax*ay*b*x0 - 120*ax*ay*cx*cy + 120*ax*ay*cx*y0 + 120*ax*ay*cy*x0 - 120*ax*ay*x0*y0 + 72*ay*ay*ay*ay*b*b - 90*ay*ay*ay*ay*b*cy + 90*ay*ay*ay*ay*b*y0 + 30*ay*ay*ay*ay*cy*cy - 60*ay*ay*ay*ay*cy*y0 + 30*ay*ay*ay*ay*y0*y0 - 144*ay*ay*b*b + 180*ay*ay*b*cy - 180*ay*ay*b*y0 - 60*ay*ay*cy*cy + 120*ay*ay*cy*y0 - 60*ay*ay*y0*y0 + 72*b*b - 90*b*cy + 90*b*y0 + 30*cx*cx - 60*cx*x0 + 30*cy*cy - 60*cy*y0 + 30*x0*x0 + 30*y0*y0))/(90*a*b)
};

    return phi;
}

std::vector<double> Q4S::get_phi_grad(const double t, const double beta) const{
    std::vector<double> phi{
    4*a*b*beta*t/9
    ,
    2*a*b*beta*t/9
    ,
    a*b*beta*t/9
    ,
    2*a*b*beta*t/9
    ,
    2*a*b*beta*t/9
    ,
    4*a*b*beta*t/9
    ,
    2*a*b*beta*t/9
    ,
    a*b*beta*t/9
    ,
    a*b*beta*t/9
    ,
    2*a*b*beta*t/9
    ,
    4*a*b*beta*t/9
    ,
    2*a*b*beta*t/9
    ,
    2*a*b*beta*t/9
    ,
    a*b*beta*t/9
    ,
    2*a*b*beta*t/9
    ,
    4*a*b*beta*t/9
    };

    return phi;
}

std::vector<double> Q4S::get_phi_unidirectional(const double t, const double beta, const double l, const std::vector<double>& v, const double vn) const{
    std::vector<double> phi{
    t*(-3*a*a*l*l + a*b*(-4*a*b*beta - 3*a*l*v[1]*vn - 3*b*l*v[0]*vn) - 3*b*b*l*l)/(9*a*b)
    ,
    t*(-3*a*a*l*l + a*b*(-4*a*b*beta - 3*a*l*v[1]*vn + 6*b*l*v[0]*vn) + 6*b*b*l*l)/(18*a*b)
    ,
    -a*b*beta*t/9 + a*l*t*v[1]*vn/6 + a*l*l*t/(6*b) + b*l*t*v[0]*vn/6 + b*l*l*t/(6*a)
    ,
    t*(6*a*a*l*l + a*b*(-4*a*b*beta + 6*a*l*v[1]*vn - 3*b*l*v[0]*vn) - 3*b*b*l*l)/(18*a*b)
    ,
    t*(-3*a*a*l*l - a*b*(4*a*b*beta + 3*a*l*v[1]*vn + 6*b*l*v[0]*vn) + 6*b*b*l*l)/(18*a*b)
    ,
    t*(-3*a*a*l*l + a*b*(-4*a*b*beta - 3*a*l*v[1]*vn + 3*b*l*v[0]*vn) - 3*b*b*l*l)/(9*a*b)
    ,
    t*(6*a*a*l*l + a*b*(-4*a*b*beta + 6*a*l*v[1]*vn + 3*b*l*v[0]*vn) - 3*b*b*l*l)/(18*a*b)
    ,
    -a*b*beta*t/9 + a*l*t*v[1]*vn/6 + a*l*l*t/(6*b) - b*l*t*v[0]*vn/6 + b*l*l*t/(6*a)
    ,
    -a*b*beta*t/9 - a*l*t*v[1]*vn/6 + a*l*l*t/(6*b) - b*l*t*v[0]*vn/6 + b*l*l*t/(6*a)
    ,
    t*(6*a*a*l*l + a*b*(-4*a*b*beta - 6*a*l*v[1]*vn + 3*b*l*v[0]*vn) - 3*b*b*l*l)/(18*a*b)
    ,
    t*(-3*a*a*l*l + a*b*(-4*a*b*beta + 3*a*l*v[1]*vn + 3*b*l*v[0]*vn) - 3*b*b*l*l)/(9*a*b)
    ,
    t*(-3*a*a*l*l + a*b*(-4*a*b*beta + 3*a*l*v[1]*vn - 6*b*l*v[0]*vn) + 6*b*b*l*l)/(18*a*b)
    ,
    t*(6*a*a*l*l - a*b*(4*a*b*beta + 6*a*l*v[1]*vn + 3*b*l*v[0]*vn) - 3*b*b*l*l)/(18*a*b)
    ,
    -a*b*beta*t/9 - a*l*t*v[1]*vn/6 + a*l*l*t/(6*b) + b*l*t*v[0]*vn/6 + b*l*l*t/(6*a)
    ,
    t*(-3*a*a*l*l + a*b*(-4*a*b*beta + 3*a*l*v[1]*vn + 6*b*l*v[0]*vn) + 6*b*b*l*l)/(18*a*b)
    ,
    t*(-3*a*a*l*l + a*b*(-4*a*b*beta + 3*a*l*v[1]*vn - 3*b*l*v[0]*vn) - 3*b*b*l*l)/(9*a*b)
    };

    return phi;
}

std::vector<double> Q4S::helmholtz_tensor(const double t, const double r) const{
    std::vector<double> h{
    4*a*b*t/9 + a*r*r*t/(3*b) + b*r*r*t/(3*a)
    ,
    2*a*b*t/9 + a*r*r*t/(6*b) - b*r*r*t/(3*a)
    ,
    a*b*t/9 - a*r*r*t/(6*b) - b*r*r*t/(6*a)
    ,
    2*a*b*t/9 - a*r*r*t/(3*b) + b*r*r*t/(6*a)
    ,
    2*a*b*t/9 + a*r*r*t/(6*b) - b*r*r*t/(3*a)
    ,
    4*a*b*t/9 + a*r*r*t/(3*b) + b*r*r*t/(3*a)
    ,
    2*a*b*t/9 - a*r*r*t/(3*b) + b*r*r*t/(6*a)
    ,
    a*b*t/9 - a*r*r*t/(6*b) - b*r*r*t/(6*a)
    ,
    a*b*t/9 - a*r*r*t/(6*b) - b*r*r*t/(6*a)
    ,
    2*a*b*t/9 - a*r*r*t/(3*b) + b*r*r*t/(6*a)
    ,
    4*a*b*t/9 + a*r*r*t/(3*b) + b*r*r*t/(3*a)
    ,
    2*a*b*t/9 + a*r*r*t/(6*b) - b*r*r*t/(3*a)
    ,
    2*a*b*t/9 - a*r*r*t/(3*b) + b*r*r*t/(6*a)
    ,
    a*b*t/9 - a*r*r*t/(6*b) - b*r*r*t/(6*a)
    ,
    2*a*b*t/9 + a*r*r*t/(6*b) - b*r*r*t/(3*a)
    ,
    4*a*b*t/9 + a*r*r*t/(3*b) + b*r*r*t/(3*a)
    };

    return h;
}

std::vector<double> Q4S::helmholtz_vector(const double t) const{
    const double V = 4*a*b*t;
    const double N = V/4;

    return std::vector<double>(NODES_PER_ELEM, N);
}

std::vector<double> Q4S::get_nodal_density_gradient(gp_Pnt p) const{
    const gp_Pnt pnorm = this->normalize(p);

    const double xi = pnorm.X();
    const double eta = pnorm.Y();

    return std::vector<double>{-(b-eta)/(4*a*b), (b-eta)/(4*a*b), (b+eta)/(4*a*b), -(b+eta)/(4*a*b),
                               -(a-xi)/(4*a*b), -(a-xi)/(4*a*b), (a+xi)/(4*a*b), -(b+eta)/(4*a*b)};
}

}


