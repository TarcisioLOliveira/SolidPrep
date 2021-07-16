/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   Foobar is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Foobar is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with SolidPrep.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef BEAM_LINEAR_2D_HPP
#define BEAM_LINEAR_2D_HPP

#include "element.hpp"
#include "utils.hpp"

namespace element{

class BeamLinear2D : public BeamElement{
    public:
    BeamLinear2D(BeamNode* n1, BeamNode* n2, double I, double A, double E);

    virtual std::vector<double> get_k() const override;
    virtual BeamNode* get_internal_loads(size_t node, const std::vector<double>& u) const override;

    const double I, A, E;
};

}

#endif
