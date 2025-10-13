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
#include "utils.hpp"
#include "element_factory.hpp"
#include "element_common.hpp"

namespace element{

class Q4S : public MeshElementCommon2DQuad<Q4S>{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 3;
    static const size_t NODE_DOF       = 2;
    static const size_t NODES_PER_ELEM = 4;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;
    static const size_t INTEG_ORDER    = 1;

    static const size_t BOUNDARY_NODES_PER_ELEM = 2;
    static const size_t BOUNDARY_GMSH_TYPE = 1;
    static inline std::unique_ptr<BoundaryMeshElementFactory> get_boundary_element_info(){
        logger::log_assert(false, logger::ERROR, "LINE ELEMENT TYPE NOT IMPLEMENTED");
        return std::unique_ptr<BoundaryMeshElementFactory>();
    }
    static std::unique_ptr<ContactMeshElementFactory> get_contact_element_info(){
        logger::log_assert(false, logger::ERROR, "LINE ELEMENT TYPE NOT IMPLEMENTED");
        return std::unique_ptr<ContactMeshElementFactory>();
    }

    static const spview::defs::ElementType SPVIEW_CODE = spview::defs::Q4;

    Q4S(ElementShape s);

    virtual math::Matrix get_k(const math::Matrix& D, const double t) const override;
    virtual math::Matrix get_nodal_density_gradient(gp_Pnt p) const override;
    virtual math::Matrix get_R(const math::Matrix& K, const double t, const std::vector<gp_Pnt>& points) const override;
    virtual math::Matrix get_B(const gp_Pnt& point) const override;

    virtual math::Matrix diffusion_1dof(const double t, const math::Matrix& A) const override;
    virtual math::Matrix advection_1dof(const double t, const math::Vector& v) const override;
    virtual math::Matrix absorption_1dof(const double t) const override;
    virtual math::Vector source_1dof(const double t) const override;
    virtual math::Vector flow_1dof(const double t, const MeshNode** nodes) const override;

    virtual math::Matrix get_Ni(const gp_Pnt& p) const override;

    virtual inline std::unique_ptr<MeshElementFactory> get_element_info() const override{
        return std::unique_ptr<MeshElementFactory>(new MeshElementFactoryImpl<Q4S>());
    }

    private:
    static const bool reg;
    virtual math::Matrix get_Nf(const double t, const std::vector<gp_Pnt>& points) const override;

    double a, b, x0, y0 = 0;

    inline double N(double x, double y, size_t i) const{
        switch(i){
            case 0:
                return (a-x)*(b-y)/(4*a*b);
            case 1:
                return (a+x)*(b-y)/(4*a*b);
            case 2:
                return (a+x)*(b+y)/(4*a*b);
            case 3:
                return (a-x)*(b+y)/(4*a*b);
        }
        return 0;
    }

    inline gp_Pnt normalize(const gp_Pnt& point) const{
        return gp_Pnt(
            point.X() - x0 - a,
            point.Y() - y0 - b,
             0);
    }
};

}

#endif
