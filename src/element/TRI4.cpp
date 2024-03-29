/*
 *   Copyright (C) 2023 Tarcísio Ladeia de Oliveira.
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

#include "element/TRI4.hpp"
#include "cblas.h"
#include "logger.hpp"
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS.hxx>
#include <IntTools_EdgeEdge.hxx>
#include <lapacke.h>

namespace element{

TRI4::TRI4(ElementShape s):
    MeshElementCommon2DTri<TRI4>(s.nodes){

    constexpr size_t N = TRI4::NODES_PER_ELEM;
    std::array<double, N> x, y;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
    }

    std::array<double, N*N> M = 
        {1, x[0], y[0], x[0]*y[0],
         1, x[1], y[1], x[1]*y[1],
         1, x[2], y[2], x[2]*y[2],
         1, x[3], y[3], x[3]*y[3]};

    std::array<int, N> ipiv;

    // M*C = I -> C=M^-1
    // C = {a[0], a[1], a[2], a[3],
    //      b[0], b[1], b[2], b[3],
    //      c[0], c[1], c[2], c[3],
    //      d[0], d[1], d[2], d[3]}
    int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating LU in Q4.", info);
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating computing inverse from LU in Q4.", info);

    this->delta = this->get_volume(1.0);

    std::copy(M.begin(), M.begin()+4, a);
    std::copy(M.begin()+4, M.begin()+8, b);
    std::copy(M.begin()+8, M.begin()+12, c);
    std::copy(M.begin()+12, M.begin()+16, d);
}

std::vector<double> TRI4::get_k(const std::vector<double>& D, const double t) const{
    std::vector<double> k{
    3*delta*t*(6*D[0] + 5*D[2] + 5*D[6] + 6*D[8])/4
    ,
    3*delta*t*(5*D[1] + 6*D[2] + 6*D[7] + 5*D[8])/4
    ,
    delta*t*(2*D[0] + 7*D[2] - D[6] + 10*D[8])/4
    ,
    delta*t*(7*D[1] + 2*D[2] + 10*D[7] - D[8])/4
    ,
    delta*t*(10*D[0] - D[2] + 7*D[6] + 2*D[8])/4
    ,
    delta*t*(-D[1] + 10*D[2] + 2*D[7] + 7*D[8])/4
    ,
    3*delta*t*(-10*D[0] - 7*D[2] - 7*D[6] - 10*D[8])/4
    ,
    3*delta*t*(-7*D[1] - 10*D[2] - 10*D[7] - 7*D[8])/4
    ,
    3*delta*t*(5*D[3] + 6*D[5] + 6*D[6] + 5*D[8])/4
    ,
    3*delta*t*(6*D[4] + 5*D[5] + 5*D[7] + 6*D[8])/4
    ,
    delta*t*(-D[3] + 10*D[5] + 2*D[6] + 7*D[8])/4
    ,
    delta*t*(10*D[4] - D[5] + 7*D[7] + 2*D[8])/4
    ,
    delta*t*(7*D[3] + 2*D[5] + 10*D[6] - D[8])/4
    ,
    delta*t*(2*D[4] + 7*D[5] - D[7] + 10*D[8])/4
    ,
    3*delta*t*(-7*D[3] - 10*D[5] - 10*D[6] - 7*D[8])/4
    ,
    3*delta*t*(-10*D[4] - 7*D[5] - 7*D[7] - 10*D[8])/4
    ,
    delta*t*(2*D[0] - D[2] + 7*D[6] + 10*D[8])/4
    ,
    delta*t*(-D[1] + 2*D[2] + 10*D[7] + 7*D[8])/4
    ,
    delta*t*(2*D[0] - D[2] - D[6] + 6*D[8])/4
    ,
    delta*t*(-D[1] + 2*D[2] + 6*D[7] - D[8])/4
    ,
    delta*t*(2*D[0] - D[2] + 3*D[6] + 2*D[8])/4
    ,
    delta*t*(-D[1] + 2*D[2] + 2*D[7] + 3*D[8])/4
    ,
    3*delta*t*(-2*D[0] + D[2] - 3*D[6] - 6*D[8])/4
    ,
    3*delta*t*(D[1] - 2*D[2] - 6*D[7] - 3*D[8])/4
    ,
    delta*t*(7*D[3] + 10*D[5] + 2*D[6] - D[8])/4
    ,
    delta*t*(10*D[4] + 7*D[5] - D[7] + 2*D[8])/4
    ,
    delta*t*(-D[3] + 6*D[5] + 2*D[6] - D[8])/4
    ,
    delta*t*(6*D[4] - D[5] - D[7] + 2*D[8])/4
    ,
    delta*t*(3*D[3] + 2*D[5] + 2*D[6] - D[8])/4
    ,
    delta*t*(2*D[4] + 3*D[5] - D[7] + 2*D[8])/4
    ,
    3*delta*t*(-3*D[3] - 6*D[5] - 2*D[6] + D[8])/4
    ,
    3*delta*t*(-6*D[4] - 3*D[5] + D[7] - 2*D[8])/4
    ,
    delta*t*(10*D[0] + 7*D[2] - D[6] + 2*D[8])/4
    ,
    delta*t*(7*D[1] + 10*D[2] + 2*D[7] - D[8])/4
    ,
    delta*t*(2*D[0] + 3*D[2] - D[6] + 2*D[8])/4
    ,
    delta*t*(3*D[1] + 2*D[2] + 2*D[7] - D[8])/4
    ,
    delta*t*(6*D[0] - D[2] - D[6] + 2*D[8])/4
    ,
    delta*t*(-D[1] + 6*D[2] + 2*D[7] - D[8])/4
    ,
    3*delta*t*(-6*D[0] - 3*D[2] + D[6] - 2*D[8])/4
    ,
    3*delta*t*(-3*D[1] - 6*D[2] - 2*D[7] + D[8])/4
    ,
    delta*t*(-D[3] + 2*D[5] + 10*D[6] + 7*D[8])/4
    ,
    delta*t*(2*D[4] - D[5] + 7*D[7] + 10*D[8])/4
    ,
    delta*t*(-D[3] + 2*D[5] + 2*D[6] + 3*D[8])/4
    ,
    delta*t*(2*D[4] - D[5] + 3*D[7] + 2*D[8])/4
    ,
    delta*t*(-D[3] + 2*D[5] + 6*D[6] - D[8])/4
    ,
    delta*t*(2*D[4] - D[5] - D[7] + 6*D[8])/4
    ,
    3*delta*t*(D[3] - 2*D[5] - 6*D[6] - 3*D[8])/4
    ,
    3*delta*t*(-2*D[4] + D[5] - 3*D[7] - 6*D[8])/4
    ,
    3*delta*t*(-10*D[0] - 7*D[2] - 7*D[6] - 10*D[8])/4
    ,
    3*delta*t*(-7*D[1] - 10*D[2] - 10*D[7] - 7*D[8])/4
    ,
    3*delta*t*(-2*D[0] - 3*D[2] + D[6] - 6*D[8])/4
    ,
    3*delta*t*(-3*D[1] - 2*D[2] - 6*D[7] + D[8])/4
    ,
    3*delta*t*(-6*D[0] + D[2] - 3*D[6] - 2*D[8])/4
    ,
    3*delta*t*(D[1] - 6*D[2] - 2*D[7] - 3*D[8])/4
    ,
    27*delta*t*(2*D[0] + D[2] + D[6] + 2*D[8])/4
    ,
    27*delta*t*(D[1] + 2*D[2] + 2*D[7] + D[8])/4
    ,
    3*delta*t*(-7*D[3] - 10*D[5] - 10*D[6] - 7*D[8])/4
    ,
    3*delta*t*(-10*D[4] - 7*D[5] - 7*D[7] - 10*D[8])/4
    ,
    3*delta*t*(D[3] - 6*D[5] - 2*D[6] - 3*D[8])/4
    ,
    3*delta*t*(-6*D[4] + D[5] - 3*D[7] - 2*D[8])/4
    ,
    3*delta*t*(-3*D[3] - 2*D[5] - 6*D[6] + D[8])/4
    ,
    3*delta*t*(-2*D[4] - 3*D[5] + D[7] - 6*D[8])/4
    ,
    27*delta*t*(D[3] + 2*D[5] + 2*D[6] + D[8])/4
    ,
    27*delta*t*(2*D[4] + D[5] + D[7] + 2*D[8])/4
    };

    return k;
}

std::vector<double> TRI4::get_DB(const std::vector<double>& D, const gp_Pnt& point) const{
    const double x = point.X();
    const double y = point.Y();

    std::vector<double> DB{
    D[0]*b[0] + D[0]*d[0]*y + D[2]*c[0] + D[2]*d[0]*x
    ,
    D[1]*c[0] + D[1]*d[0]*x + D[2]*b[0] + D[2]*d[0]*y
    ,
    D[0]*b[1] + D[0]*d[1]*y + D[2]*c[1] + D[2]*d[1]*x
    ,
    D[1]*c[1] + D[1]*d[1]*x + D[2]*b[1] + D[2]*d[1]*y
    ,
    D[0]*b[2] + D[0]*d[2]*y + D[2]*c[2] + D[2]*d[2]*x
    ,
    D[1]*c[2] + D[1]*d[2]*x + D[2]*b[2] + D[2]*d[2]*y
    ,
    D[0]*b[3] + D[0]*d[3]*y + D[2]*c[3] + D[2]*d[3]*x
    ,
    D[1]*c[3] + D[1]*d[3]*x + D[2]*b[3] + D[2]*d[3]*y
    ,
    D[3]*b[0] + D[3]*d[0]*y + D[5]*c[0] + D[5]*d[0]*x
    ,
    D[4]*c[0] + D[4]*d[0]*x + D[5]*b[0] + D[5]*d[0]*y
    ,
    D[3]*b[1] + D[3]*d[1]*y + D[5]*c[1] + D[5]*d[1]*x
    ,
    D[4]*c[1] + D[4]*d[1]*x + D[5]*b[1] + D[5]*d[1]*y
    ,
    D[3]*b[2] + D[3]*d[2]*y + D[5]*c[2] + D[5]*d[2]*x
    ,
    D[4]*c[2] + D[4]*d[2]*x + D[5]*b[2] + D[5]*d[2]*y
    ,
    D[3]*b[3] + D[3]*d[3]*y + D[5]*c[3] + D[5]*d[3]*x
    ,
    D[4]*c[3] + D[4]*d[3]*x + D[5]*b[3] + D[5]*d[3]*y
    ,
    D[6]*b[0] + D[6]*d[0]*y + D[8]*c[0] + D[8]*d[0]*x
    ,
    D[7]*c[0] + D[7]*d[0]*x + D[8]*b[0] + D[8]*d[0]*y
    ,
    D[6]*b[1] + D[6]*d[1]*y + D[8]*c[1] + D[8]*d[1]*x
    ,
    D[7]*c[1] + D[7]*d[1]*x + D[8]*b[1] + D[8]*d[1]*y
    ,
    D[6]*b[2] + D[6]*d[2]*y + D[8]*c[2] + D[8]*d[2]*x
    ,
    D[7]*c[2] + D[7]*d[2]*x + D[8]*b[2] + D[8]*d[2]*y
    ,
    D[6]*b[3] + D[6]*d[3]*y + D[8]*c[3] + D[8]*d[3]*x
    ,
    D[7]*c[3] + D[7]*d[3]*x + D[8]*b[3] + D[8]*d[3]*y
    };

    return DB;
}

std::vector<double> TRI4::get_Nf(const double t, const std::vector<gp_Pnt>& points) const{
    const double x[]{points[0].X(), points[1].X()};
    const double y[]{points[0].Y(), points[1].Y()};

    std::vector<double> Nf{
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[0] + 3*b[0]*x[0] + 3*b[0]*x[1] + 3*c[0]*y[0] + 3*c[0]*y[1] + 2*d[0]*x[0]*y[0] + d[0]*x[0]*y[1] +d[0]*x[1]*y[0] + 2*d[0]*x[1]*y[1])/6
    ,
    0
    ,
    0
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[0] + 3*b[0]*x[0] + 3*b[0]*x[1] + 3*c[0]*y[0] + 3*c[0]*y[1] + 2*d[0]*x[0]*y[0] + d[0]*x[0]*y[1] +d[0]*x[1]*y[0] + 2*d[0]*x[1]*y[1])/6
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[1] + 3*b[1]*x[0] + 3*b[1]*x[1] + 3*c[1]*y[0] + 3*c[1]*y[1] + 2*d[1]*x[0]*y[0] + d[1]*x[0]*y[1] +d[1]*x[1]*y[0] + 2*d[1]*x[1]*y[1])/6
    ,
    0
    ,
    0
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[1] + 3*b[1]*x[0] + 3*b[1]*x[1] + 3*c[1]*y[0] + 3*c[1]*y[1] + 2*d[1]*x[0]*y[0] + d[1]*x[0]*y[1] +d[1]*x[1]*y[0] + 2*d[1]*x[1]*y[1])/6
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[2] + 3*b[2]*x[0] + 3*b[2]*x[1] + 3*c[2]*y[0] + 3*c[2]*y[1] + 2*d[2]*x[0]*y[0] + d[2]*x[0]*y[1] +d[2]*x[1]*y[0] + 2*d[2]*x[1]*y[1])/6
    ,
    0
    ,
    0
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[2] + 3*b[2]*x[0] + 3*b[2]*x[1] + 3*c[2]*y[0] + 3*c[2]*y[1] + 2*d[2]*x[0]*y[0] + d[2]*x[0]*y[1] +d[2]*x[1]*y[0] + 2*d[2]*x[1]*y[1])/6
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[3] + 3*b[3]*x[0] + 3*b[3]*x[1] + 3*c[3]*y[0] + 3*c[3]*y[1] + 2*d[3]*x[0]*y[0] + d[3]*x[0]*y[1] +d[3]*x[1]*y[0] + 2*d[3]*x[1]*y[1])/6
    ,
    0
    ,
    0
    ,
    t*std::sqrt(x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + y[0]*y[0] - 2*y[0]*y[1] + y[1]*y[1])*(6*a[3] + 3*b[3]*x[0] + 3*b[3]*x[1] + 3*c[3]*y[0] + 3*c[3]*y[1] + 2*d[3]*x[0]*y[0] + d[3]*x[0]*y[1] +d[3]*x[1]*y[0] + 2*d[3]*x[1]*y[1])/6
    };
    return Nf;
}

}
