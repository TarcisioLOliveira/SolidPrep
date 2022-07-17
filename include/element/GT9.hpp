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

#ifndef GT9_HPP
#define GT9_HPP

#include "element.hpp"
#include "material.hpp"
#include "utils.hpp"
#include "element_factory.hpp"
#include <memory>
#include "element_common.hpp"

class ProjectData;

namespace element{

class GT9 : public MeshElementCommon2DTri<GT9>{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 2;
    static const size_t NODE_DOF       = 3;
    static const size_t NODES_PER_ELEM = 3;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;

    GT9(ElementShape s, ProjectData* data);

    virtual std::vector<double> get_k() const override;

    virtual inline std::unique_ptr<MeshElementFactory> get_element_info() const override{
        return std::unique_ptr<MeshElementFactory>(new MeshElementFactoryImpl<GT9>());
    }

    private:
    virtual std::vector<double> get_DB(const gp_Pnt& point) const override;
    virtual std::vector<double> get_Nf(const std::vector<gp_Pnt>& points) const override;
};

}

#endif
