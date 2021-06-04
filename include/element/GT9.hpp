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

#ifndef GT9_HPP
#define GT9_HPP

#include "element.hpp"
#include "material.hpp"

namespace element{

class GT9 : public MeshElement{
    public:
    GT9(ElementShape s, Material* mat, std::vector<long> u_pos, double thickness);

    virtual std::vector<double> get_k() const override;
    virtual MeshNode* get_stresses(size_t node, const std::vector<double>& u) const override;
    virtual MeshNode* get_internal_loads(size_t node, const std::vector<double>& u) const override;

    private:
    Material const * const mat;
    double t;
};

}

#endif
