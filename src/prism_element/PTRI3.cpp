/*
 *   Copyright (C) 2026 Tarcísio Ladeia de Oliveira.
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

#include "prism_element/PTRI3.hpp"

namespace prism_element{

PTRI3::PTRI3(const std::vector<MeshNode*> nodes, const MeshElement* const parent):
    PrismMeshElementCommon(nodes, parent)
{

    this->C = math::Matrix(3, 3, 0);
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        this->C(i, 0) = 1;
        this->C(i, 1) = this->nodes[i]->point.X();
        this->C(i, 2) = this->nodes[i]->point.Y();
    }
    this->A = 0.5*C.determinant();
    this->C.invert_LU();
}

math::Matrix PTRI3::get_Ni(const gp_Pnt& p) const{
    const auto Ni = this->get_Ni_1dof(p);

    return math::Matrix({Ni[0], 0, 0, Ni[1], 0, 0, Ni[2], 0, 0,
                         0, Ni[0], 0, 0, Ni[1], 0, 0, Ni[2], 0,
                         0, 0, Ni[0], 0, 0, Ni[1], 0, 0, Ni[2]}, NODE_DOF, K_DIM);
}
math::Vector PTRI3::get_Ni_1dof(const gp_Pnt& p) const{
    return math::Vector({N(p, 0, this->C), N(p, 1, this->C), N(p, 2, this->C)});
}
math::Matrix PTRI3::get_nodal_density_gradient(const gp_Pnt& p) const{
    (void) p;
    const double* const a = C.data();
    const double* const b = a + NODES_PER_ELEM;
    const double* const c = b + NODES_PER_ELEM;
    return math::Matrix({b[0], b[1], b[2],
                         c[0], c[1], c[2]}, 2, 3);
}

}
