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


#ifndef SUB_PROBLEM_HPP
#define SUB_PROBLEM_HPP

#include "support.hpp"
#include "force.hpp"
#include "spring.hpp"
#include "internal_loads.hpp"

class SubProblem{
    public:
    std::vector<Support*> supports;
    std::vector<Force*> forces;
    std::vector<Spring*> springs;
    std::vector<InternalLoads*> internal_loads;
};

#endif
