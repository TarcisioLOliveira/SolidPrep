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

#include "logger.hpp"
#include "spring.hpp"
#include "utils.hpp"

Spring::Spring(CrossSection cross_section, gp_Dir normal, gp_Dir v, gp_Dir w, Material* mat, std::array<double, 3> L, std::array<double, 3> F, std::array<double, 3> curv, MeshElementFactory* bound_elem, utils::ProblemType type):
    S(std::move(cross_section)), 
    rot2D{{normal.X(), v.X()},
          {normal.Y(), v.Y()}},
    rot3D{{normal.X(), v.X(), w.X()},
         {normal.Y(), v.Y(), w.Y()},
         {normal.Z(), v.Z(), w.Z()}},
    F(F), curv(curv),
    mat(mat),
    normal(normal), 
    boundary_elem_info(bound_elem),
    v(v), w(w), L(L), 
    type(type){

}

std::vector<double> Spring::get_K(const gp_Pnt& p) const{
    if(type == utils::PROBLEM_TYPE_2D){

        auto EG = mat->beam_EG_2D(p, this->normal);

        Eigen::Matrix<double, 2, 2> Korig{{EG[0]/L[0], 0},
                                          {0, EG[1]/L[1]}};

        Eigen::Matrix<double, 2, 2> Ktmp = rot2D*Korig*rot2D.transpose();

        std::vector<double> K(4);
        std::copy(Ktmp.data(), Ktmp.data()+4, K.begin());

        return K;
    } else if(type == utils::PROBLEM_TYPE_3D){

        auto EG = mat->beam_EG_3D(p, this->normal);

        Eigen::Matrix<double, 3, 3> Korig{{EG[0]/L[0], 0, 0},
                                          {0, EG[1]/L[1], 0},
                                          {0, 0, EG[2]/L[2]}};

        Eigen::Matrix<double, 3, 3> Ktmp = rot3D*Korig*rot3D.transpose();

        std::vector<double> K(9);
        std::copy(Ktmp.data(), Ktmp.data()+9, K.begin());

        return K;
    }

    return std::vector<double>();
}
