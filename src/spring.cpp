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

#include "spring.hpp"
#include "utils.hpp"

Spring::Spring(CrossSection cross_section, gp_Dir normal, Material* mat, std::array<double, 3> L, utils::ProblemType type):
    S(std::move(cross_section)), normal(normal), mat(mat), L(L), type(type){
    
}

std::vector<double> Spring::get_K(const gp_Pnt& p) const{
    if(type == utils::PROBLEM_TYPE_2D){

        auto EG = mat->beam_EG_2D(p, this->normal);

        return std::vector<double>{EG[0]/L[0], 0,
                                   0, EG[1]/L[1]};
    } else if(type == utils::PROBLEM_TYPE_3D){

        auto EG = mat->beam_EG_3D(p, this->normal);

        return std::vector<double>{EG[0]/L[0], 0, 0,
                                   0, EG[1]/L[1], 0,
                                   0, 0, EG[2]/L[2]};
    }

    return std::vector<double>();
}
