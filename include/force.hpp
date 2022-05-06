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

#ifndef FORCE_HPP
#define FORCE_HPP
#include <vector>
#include <gp_Vec.hxx>
#include "cross_section.hpp"

class Force{
    public:

    /**
     * Creates a Force object.
     *
     * @param cross_section Region of force application.
     * @param force Force vector.
     */
    Force(CrossSection cross_section, gp_Vec force);

    const CrossSection S;
    const gp_Vec vec;
};

#endif
