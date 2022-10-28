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

#ifndef TET4_HPP
#define TET4_HPP

#include "element.hpp"
#include "material.hpp"
#include <vector>
#include "utils.hpp"
#include "element_factory.hpp"
#include "element_common.hpp"
#include <array>

namespace element{

class TET4 : public MeshElementCommon3DTet<TET4>{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 4;
    static const size_t NODE_DOF       = 3;
    static const size_t NODES_PER_ELEM = 4;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;

    static const size_t BOUNDARY_NODES_PER_ELEM = 3;
    static const size_t BOUNDARY_GMSH_TYPE = 2;

    TET4(ElementShape s);

    virtual std::vector<double> get_k(const std::vector<double>& D, const double t) const override;

    virtual inline std::unique_ptr<MeshElementFactory> get_element_info() const override{
        return std::unique_ptr<MeshElementFactory>(new MeshElementFactoryImpl<TET4>());
    }

    private:
    virtual std::vector<double> get_DB(const std::vector<double>& D, const gp_Pnt& point) const override;
    virtual std::vector<double> get_Nf(const double t, const std::vector<gp_Pnt>& points) const override;

    std::array<double, NODES_PER_ELEM*NODES_PER_ELEM> get_coeffs() const;

    const std::array<double, NODES_PER_ELEM*NODES_PER_ELEM> coeffs;
};

}

#endif
