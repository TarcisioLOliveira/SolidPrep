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
#ifndef SPRING_HPP
#define SPRING_HPP

#include <array>
#include "material.hpp"
#include "utils.hpp"
#include "cross_section.hpp"

class Spring{
    public:

    Spring(CrossSection cross_section, gp_Dir normal, Material* mat, std::array<double, 3> L, utils::ProblemType type);

    const CrossSection S;
    std::vector<double> get_K(const gp_Pnt& p) const;

    private:
    const gp_Dir normal;
    const Material* mat;
    const std::array<double, 3> L;
    const utils::ProblemType type;
    std::vector<double> generate_K(gp_Dir normal, Material* mat, std::array<double, 3> L, utils::ProblemType type) const;
};

#endif