/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
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

#ifndef FINITE_ELEMENT_HPP
#define FINITE_ELEMENT_HPP

#include <memory>
#include <vector>
#include "element.hpp"
#include "meshing.hpp"
#include "element_factory.hpp"
#include "geometry.hpp"

class ProjectData;

class FiniteElement{
    public:
    const double K_MIN = 1e-6;

    virtual ~FiniteElement() = default;

    virtual std::vector<double> calculate_displacements(const Meshing* const mesh, std::vector<double> load, const std::vector<double>& density = std::vector<double>(), double pc = 3) = 0;

    inline virtual void set_steps(size_t s){
        this->steps = s;
    }

    virtual std::vector<double> calculate_forces(const Meshing* const mesh, const std::vector<double>& displacements) const;

    protected:
    size_t steps = 1;
    size_t current_step = 0;
};

#endif
