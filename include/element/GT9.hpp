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
#include "logger.hpp"
#include "material.hpp"
#include "math/matrix.hpp"
#include "utils.hpp"
#include "element_factory.hpp"
#include <memory>
#include "element_common.hpp"

namespace element{

class GT9 : public MeshElementCommon2DTri<GT9>{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 2;
    static const size_t NODE_DOF       = 3;
    static const size_t NODES_PER_ELEM = 3;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;
    static const size_t INTEG_ORDER    = 2;

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
    static std::unique_ptr<ShapeMeshElementFactory> get_shape_element_info(){
        logger::log_assert(false, logger::ERROR, "LINE ELEMENT TYPE NOT IMPLEMENTED");
        return std::unique_ptr<ShapeMeshElementFactory>();
    }

    static const spview::defs::ElementType SPVIEW_CODE = spview::defs::TRI3;

    GT9(ElementShape s);

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

    virtual math::Matrix get_uu(const MeshElement* const e2, const std::vector<gp_Pnt>& bounds, const gp_Dir n) const override;
    virtual math::Matrix get_MnMn(const MeshElement* const e2, const std::vector<double>& u_ext, const std::vector<gp_Pnt>& bounds, const gp_Dir n) const override;

    virtual inline std::unique_ptr<MeshElementFactory> get_element_info() const override{
        return std::unique_ptr<MeshElementFactory>(new MeshElementFactoryImpl<GT9>());
    }

    private:
    static const bool reg;
    virtual math::Matrix get_Nf(const double t, const std::vector<gp_Pnt>& points) const override;

    double a[3], b[3], c[3], delta;

    inline double N(double x, double y, size_t i) const{
        return (a[i] + b[i]*x + c[i]*y)/(2*delta);
    }

    inline math::Matrix N_mat(double x, double y) const{
        double Ni[NODE_DOF] = {N(x,y,0), N(x,y,1), N(x,y,2)};
        math::Matrix NN(DIM, K_DIM);
        for(size_t i = 0; i < NODE_DOF; ++i){
            const size_t j = (i+1) % NODE_DOF;
            const size_t k = (i+2) % NODE_DOF;

            const double Nut = 0.5*Ni[i]*(b[k]*Ni[j]-b[j]*Ni[k]);
            const double Nvt = 0.5*Ni[i]*(c[k]*Ni[j]-c[j]*Ni[k]);
            NN(0, i*NODE_DOF + 0) = Ni[i];
            NN(0, i*NODE_DOF + 2) = Nut;
            NN(1, i*NODE_DOF + 1) = Ni[i];
            NN(1, i*NODE_DOF + 2) = Nvt;
        }

        return NN;
    }

    inline math::Vector N_mat_1dof(double x, double y) const{
        return math::Vector{N(x, y, 0), N(x, y, 1), N(x, y, 2)};
    }
    inline math::Matrix dN_mat_1dof() const{
        return math::Matrix({b[0], b[1], b[2],
                             c[0], c[1], c[2]}, 2, 3)/(2*delta);
    }
};

}

#endif
