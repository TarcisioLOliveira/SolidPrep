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

#ifndef PTRI3_HPP
#define PTRI3_HPP

#include "element.hpp"
#include "element_common.hpp"
#include "math/matrix.hpp"

namespace prism_element{

class PTRI3 : public PrismMeshElementCommon<PTRI3>{
    public:
    static const size_t ORDER          = 1;
    static const size_t NODE_DOF       = 3;
    static const size_t NODES_PER_ELEM = 3;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;
    static const size_t INTEG_ORDER    = 1;
    static const size_t Z_ORDER        = 4;
    static const Element::Shape SHAPE_TYPE = Element::Shape::TRI;
    static const utils::ProblemType PROBLEM_TYPE = utils::PROBLEM_TYPE_3D;

    PTRI3(const std::vector<MeshNode*> nodes, const MeshElement* const parent);
    virtual ~PTRI3() = default;

    virtual math::Matrix get_Ni(const gp_Pnt& p) const override;
    virtual math::Vector get_Ni_1dof(const gp_Pnt& p) const override;
    virtual math::Matrix get_nodal_density_gradient(const gp_Pnt& p) const override;

    private:
    math::Matrix C;
    double A;

    inline double N(const gp_Pnt& p, size_t i, const math::Matrix& M) const{
        const double* const a = M.data();
        const double* const b = a + NODES_PER_ELEM;
        const double* const c = b + NODES_PER_ELEM;
        return a[i] + b[i]*p.X() + c[i]*p.Y();
    }
};

}

#endif
