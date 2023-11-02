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

#ifndef CURVATURE_HPP
#define CURVATURE_HPP

#include <memory>
#include <Eigen/Core>
#include <gp_Dir.hxx>
#include "material.hpp"
#include "element.hpp"

class Curvature{
    public:
    Curvature(const Material* mat, gp_Dir u, gp_Dir v, gp_Dir w, Eigen::Matrix<double, 2, 2> rot2D, Eigen::Matrix<double, 3, 3> rot3D, Element::Shape elem_shape);

    void generate_curvature_3D(const std::vector<std::unique_ptr<MeshElement>>& boundary_mesh, const std::array<double, 3>& M);

    private:
    const Material* mat;
    const gp_Dir u, v, w;
    const Eigen::Matrix<double, 2, 2> rot2D;
    const Eigen::Matrix<double, 3, 3> rot3D;
    const Element::Shape elem_shape;

    const std::array<double, 4*7> GS_tri_params =
        {1.0/3.0, 1.0/3.0, 1.0/3.0, 0.225,
         0.05961587, 0.47014206, 0.47014206, 0.13239415,
         0.47014206, 0.05961587, 0.47014206, 0.13239415,
         0.47014206, 0.47014206, 0.05961587, 0.13239415,
         0.79742699, 0.10128651, 0.10128651, 0.12593918,
         0.10128651, 0.79742699, 0.10128651, 0.12593918,
         0.10128651, 0.10128651, 0.79742699, 0.12593918};

    double EA;
    gp_Pnt neutral;
    double EI_v;
    double EI_w;
    double curv_v;
    double curv_w;
    std::vector<double> phi;

    double integrate_surface_3D(const std::vector<std::unique_ptr<MeshElement>>& boundary_mesh, const std::function<double(const gp_Pnt&, const gp_Pnt& px)>& fn) const;
    double GS_tri(const std::array<gp_Pnt, 3>& p, const std::array<gp_Pnt, 3>& px, const std::function<double(const gp_Pnt&, const gp_Pnt& px)>& fn) const;
    inline double make_EA_base_3D(const gp_Pnt& p, const gp_Pnt& px) const{
        (void)px;
        return this->mat->beam_E_3D(p, this->u);
    }
};

#endif
