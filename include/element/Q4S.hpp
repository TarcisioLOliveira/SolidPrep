/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
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

#ifndef Q4S_HPP
#define Q4S_HPP

#include "element.hpp"
#include "material.hpp"
#include "utils.hpp"
#include "element_factory.hpp"
#include <memory>
#include "element_common.hpp"
#include <array>

namespace element{

class Q4S : public MeshElementCommon2DQuad<Q4S>{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 3;
    static const size_t NODE_DOF       = 2;
    static const size_t NODES_PER_ELEM = 4;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;

    Q4S(ElementShape s);

    virtual std::vector<double> get_k(const std::vector<double>& D, const double t) const override;

    virtual inline std::unique_ptr<MeshElementFactory> get_element_info() const override{
        return std::unique_ptr<MeshElementFactory>(new MeshElementFactoryImpl<Q4S>());
    }

    private:
    double a, b, x0, y0 = 0;

    virtual std::vector<double> get_DB(const std::vector<double>& D, const gp_Pnt& point) const override;
    virtual std::vector<double> get_Nf(const double t, const std::vector<gp_Pnt>& points) const override;

    inline gp_Pnt normalize(const gp_Pnt& point) const{
        return gp_Pnt(
            point.X() - x0 - a,
            point.Y() - y0 - b,
             0);
    }
};

}

#endif