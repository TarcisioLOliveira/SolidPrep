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
#include "utils.hpp"
#include "cross_section.hpp"

class Spring{
    public:

    Spring(CrossSection cross_section, std::array<double, 3> K, utils::ProblemType type);

    const CrossSection S;
    const std::vector<double> K;

    private:
    std::vector<double> generate_K(std::array<double, 3> K, utils::ProblemType type) const;
};

#endif
