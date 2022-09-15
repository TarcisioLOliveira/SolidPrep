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
    MeshElementCommon2DQuad<Q4S>(s.nodes){

    constexpr size_t N = Q4S::NODES_PER_ELEM;
    double maxx = this->nodes[0]->point.X();
    double minx = this->nodes[0]->point.X();
    double maxy = this->nodes[0]->point.Y();
    double miny = this->nodes[0]->point.Y();
    for(size_t i = 1; i < N; ++i){
        double x = this->nodes[i]->point.X();
        double y = this->nodes[i]->point.Y();
        if(x > maxx){
            maxx = x;
        } else if(x < minx){
            minx = x;
        }
        if(y > maxy){
            maxy = y;
        } else if(y < miny){
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
    constexpr double eps = 1e-1;
    while(!this->nodes[0]->point.IsEqual({x0, y0, 0.0}, eps)){
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
    logger::quick_log(x[0], y[0], x[1], y[1]);

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


