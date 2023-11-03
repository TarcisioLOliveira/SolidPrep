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

#include "curvature.hpp"
#include "logger.hpp"

Curvature::Curvature(const Material* mat, gp_Dir u, gp_Dir v, gp_Dir w, Eigen::Matrix<double, 2, 2> rot2D, Eigen::Matrix<double, 3, 3> rot3D, Element::Shape elem_shape):
    mat(mat), u(u), v(v), w(w),
    rot2D(rot2D), rot3D(rot3D), 
    elem_shape(elem_shape)
{

}

void Curvature::generate_curvature_3D(const std::vector<std::unique_ptr<MeshElement>>& boundary_mesh, const std::array<double, 3>& M){
    std::array<gp_Pnt, 3> points;

    const auto fn_EA =
        [this](const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EA_base_3D(p, px);
        };

    this->EA = this->integrate_surface_3D(boundary_mesh, fn_EA);
    logger::quick_log("A: ", EA/180000);
    logger::quick_log("EA: ", EA);
}

double Curvature::integrate_surface_3D(const std::vector<std::unique_ptr<MeshElement>>& boundary_mesh, const std::function<double(const gp_Pnt&, const gp_Pnt&)>& fn) const{
    double result = 0;
    if(this->elem_shape == Element::Shape::TRI){
        Eigen::Matrix<double, 3, 3> points{{0,0,0},{0,0,0},{0,0,0}};
        #pragma omp parallel for reduction(+:result)
        for(size_t i = 0; i < boundary_mesh.size(); ++i){
            const auto& e = boundary_mesh[i];
            std::array<gp_Pnt, 3> old_points;
            for(size_t x = 0; x < 3; ++x){
                old_points[x] = e->nodes[x]->point;
                for(size_t y = 0; y < 3; ++y){
                    points(y, x) = old_points[x].Coord(y+1);
                }
            }
            Eigen::Matrix<double, 3, 3> rotd_p = this->rot3D.transpose()*points;
            std::array<gp_Pnt, 3> new_points{
                gp_Pnt(rotd_p(0, 0), rotd_p(1, 0), rotd_p(2, 0)),
                gp_Pnt(rotd_p(0, 1), rotd_p(1, 1), rotd_p(2, 1)),
                gp_Pnt(rotd_p(0, 2), rotd_p(1, 2), rotd_p(2, 2))
            };
            result += this->GS_tri(old_points, new_points, fn);
        }
    } else if(this->elem_shape == Element::Shape::QUAD){

    }
    return result;
}

double Curvature::GS_tri(const std::array<gp_Pnt, 3>& p, const std::array<gp_Pnt, 3>& px, const std::function<double(const gp_Pnt&, const gp_Pnt&)>& fn) const{
    double result = 0;
    gp_Vec v1(p[1], p[0]);
    gp_Vec v2(p[2], p[0]);
    const double drnorm = (v1.Crossed(v2)).Magnitude()/2;
    for(size_t i = 0; i < 12; ++i){
        auto P = this->GS_tri_params.data()+i*4;
        gp_Pnt pi{
            P[0]*p[0].X() + P[1]*p[1].X() + P[2]*p[2].X(),
            P[0]*p[0].Y() + P[1]*p[1].Y() + P[2]*p[2].Y(),
            P[0]*p[0].Z() + P[1]*p[1].Z() + P[2]*p[2].Z()
        };
        gp_Pnt pxi{
            P[0]*px[0].X() + P[1]*px[1].X() + P[2]*px[2].X(),
            P[0]*px[0].Y() + P[1]*px[1].Y() + P[2]*px[2].Y(),
            P[0]*px[0].Z() + P[1]*px[1].Z() + P[2]*px[2].Z()
        };
        result += P[3]*fn(pi, pxi);
    }
    result *= drnorm;

    return result;
}
