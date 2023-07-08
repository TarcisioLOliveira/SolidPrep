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

#ifndef Q4_HPP
#define Q4_HPP

#include "element.hpp"
#include "material.hpp"
#include "utils.hpp"
#include "element_factory.hpp"
#include <memory>
#include "element_common.hpp"
#include <array>

namespace element{

class Q4 : public MeshElementCommon2DQuad<Q4>{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 3;
    static const size_t NODE_DOF       = 2;
    static const size_t NODES_PER_ELEM = 4;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;

    static const size_t BOUNDARY_NODES_PER_ELEM = 2;
    static const size_t BOUNDARY_GMSH_TYPE = 1;

    static const spview::defs::ElementType SPVIEW_CODE = spview::defs::Q4;

    Q4(ElementShape s);

    virtual std::vector<double> get_k(const std::vector<double>& D, const double t) const override;
    virtual std::vector<double> helmholtz_tensor(const double t, const double r) const override;
    virtual std::vector<double> helmholtz_vector(const double t) const override;
    virtual std::vector<double> get_nodal_density_gradient(gp_Pnt p) const override;

    virtual Eigen::MatrixXd diffusion_1dof(const double t, const std::vector<double>& A) const override{}
    virtual Eigen::MatrixXd advection_1dof(const double t, const std::vector<double>& v) const override{}
    virtual Eigen::MatrixXd absorption_1dof(const double t) const override{}
    virtual Eigen::VectorXd source_1dof(const double t) const override{}

    virtual std::vector<double> get_phi_unidirectional(const double t, const double beta, const double l, const std::vector<double>& v, const double vn) const override{}

    virtual inline std::unique_ptr<MeshElementFactory> get_element_info() const override{
        return std::unique_ptr<MeshElementFactory>(new MeshElementFactoryImpl<Q4>());
    }

    private:
    virtual std::vector<double> get_DB(const std::vector<double>& D, const gp_Pnt& point) const override;
    virtual std::vector<double> get_Nf(const double t, const std::vector<gp_Pnt>& points) const override;
    virtual std::vector<double> get_B(const gp_Pnt& point) const override;

    virtual std::vector<double> get_k_base(const std::vector<double>& D, const double t, const double xi, const double eta, 
                                           const std::array<double, NODES_PER_ELEM>& x, const std::array<double, NODES_PER_ELEM>& y) const;

    virtual std::vector<double> get_h_base(const double r, const double t, const double xi, const double eta, 
                                           const std::array<double, NODES_PER_ELEM>& x, const std::array<double, NODES_PER_ELEM>& y) const;

    virtual std::vector<double> get_J_base(const double xi, const double eta, 
                                           const std::array<double, NODES_PER_ELEM>& x, const std::array<double, NODES_PER_ELEM>& y) const;
    virtual double get_detJ_base(const double xi, const double eta, 
                                 const std::array<double, NODES_PER_ELEM>& x, const std::array<double, NODES_PER_ELEM>& y) const;
    virtual std::vector<double> get_dN_base(const double xi, const double eta) const;

    virtual std::vector<double> get_h2_base(const double t, const double xi, const double eta, 
                                            const std::array<double, NODES_PER_ELEM>& x, const std::array<double, NODES_PER_ELEM>& y) const;

    std::array<double, NODES_PER_ELEM*NODES_PER_ELEM> get_coeffs() const;

    const std::array<double, NODES_PER_ELEM*NODES_PER_ELEM> coeffs;
};

}

#endif
