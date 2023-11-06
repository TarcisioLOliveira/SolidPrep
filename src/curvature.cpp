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
#include "utils/gauss_legendre.hpp"

Curvature::Curvature(const Material* mat, gp_Dir u, gp_Dir v, gp_Dir w, Eigen::Matrix<double, 2, 2> rot2D, Eigen::Matrix<double, 3, 3> rot3D, const BoundaryMeshElementFactory* elem_info, double V_v, double V_w):
    mat(mat), u(u), v(v), w(w),
    rot2D(rot2D), rot3D(rot3D), 
    elem_info(elem_info),
    V_v(V_v), V_w(V_w)
{

}

void Curvature::generate_curvature_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, size_t phi_size){
    std::array<gp_Pnt, 3> points;

    const auto fn_EA =
        [this](const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EA_base_3D(p, px);
        };
    const auto fn_EA_v =
        [this](const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EA_v_base_3D(p, px);
        };
    const auto fn_EA_w =
        [this](const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EA_w_base_3D(p, px);
        };
    const auto fn_EI_v =
        [this](const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EI_v_base_3D(p, px);
        };
    const auto fn_EI_w =
        [this](const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EI_w_base_3D(p, px);
        };

    this->EA = this->integrate_surface_3D(boundary_mesh, fn_EA);
    const double EA_v = this->integrate_surface_3D(boundary_mesh, fn_EA_v);
    const double EA_w = this->integrate_surface_3D(boundary_mesh, fn_EA_w);
    this->c_v = EA_v/this->EA;
    this->c_w = EA_w/this->EA;
    this->EI_v = this->integrate_surface_3D(boundary_mesh, fn_EI_v);
    this->EI_w = this->integrate_surface_3D(boundary_mesh, fn_EI_w);
    logger::quick_log("A: ", EA/180000);
    logger::quick_log("EA: ", EA);
    logger::quick_log("c_v: ", c_v);
    logger::quick_log("c_w: ", c_w);
    logger::quick_log("EI_v", EI_v/180000);
    logger::quick_log("EI_w", EI_w/180000);

    this->phi_torsion.resize(phi_size,0);
    this->phi_shear.resize(phi_size,0);
}

void Curvature::get_shear_in_3D(const BoundaryMeshElement* e, double& t_xz, double& t_yz) const{
    auto c = e->get_centroid();
    t_xz = 0;
    t_yz = 0;

    Eigen::Vector<double, 3> grad = this->theta*this->rot3D.transpose()*e->grad_1dof(c, this->phi_torsion);
    t_xz += grad[1];
    t_yz -= grad[0];
}

double Curvature::integrate_surface_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const std::function<double(const gp_Pnt&, const gp_Pnt&)>& fn) const{
    double result = 0;
    if(this->elem_info->get_shape_type() == Element::Shape::TRI){
        #pragma omp parallel for reduction(+:result)
        for(size_t i = 0; i < boundary_mesh.size(); ++i){
            Eigen::Matrix<double, 3, 3> points{{0,0,0},{0,0,0},{0,0,0}};
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
    } else if(this->elem_info->get_shape_type() == Element::Shape::QUAD){

    }
    return result;
}

double Curvature::GS_tri(const std::array<gp_Pnt, 3>& p, const std::array<gp_Pnt, 3>& px, const std::function<double(const gp_Pnt&, const gp_Pnt&)>& fn) const{
    double result = 0;
    gp_Vec v1(p[1], p[0]);
    gp_Vec v2(p[2], p[0]);
    const double drnorm = (v1.Crossed(v2)).Magnitude()/2;
    const auto& gsi = utils::GaussLegendreTri<6>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        gp_Pnt pi{
            it->a*p[0].X() + it->b*p[1].X() + it->c*p[2].X(),
            it->a*p[0].Y() + it->b*p[1].Y() + it->c*p[2].Y(),
            it->a*p[0].Z() + it->b*p[1].Z() + it->c*p[2].Z()
        };
        gp_Pnt pxi{
            it->a*px[0].X() + it->b*px[1].X() + it->c*px[2].X(),
            it->a*px[0].Y() + it->b*px[1].Y() + it->c*px[2].Y(),
            it->a*px[0].Z() + it->b*px[1].Z() + it->c*px[2].Z()
        };
        result += it->w*fn(pi, pxi);
    }
    result *= drnorm;

    return result;
}
