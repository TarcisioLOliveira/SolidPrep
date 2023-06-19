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

#ifndef TRI3_HPP
#define TRI3_HPP

#include <Eigen/Core>
#include "element.hpp"
#include "material.hpp"
#include <vector>
#include "utils.hpp"
#include "element_factory.hpp"
#include "element_common.hpp"

namespace element{

class TRI3 : public MeshElementCommon2DTri<TRI3>{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 2;
    static const size_t NODE_DOF       = 2;
    static const size_t NODES_PER_ELEM = 3;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;

    static const size_t BOUNDARY_NODES_PER_ELEM = 2;
    static const size_t BOUNDARY_GMSH_TYPE = 1;

    static const spview::defs::ElementType SPVIEW_CODE = spview::defs::TRI3;

    TRI3(ElementShape s);

    virtual std::vector<double> get_k(const std::vector<double>& D, const double t) const override;
    virtual std::vector<double> helmholtz_tensor(const double t, const double r) const override;
    virtual std::vector<double> helmholtz_vector(const double t) const override;
    virtual std::vector<double> get_nodal_density_gradient(gp_Pnt p) const override;
    virtual std::vector<double> get_phi_radial(const double t, const double beta, const double vp, const std::vector<double>& axis, const std::vector<double>& center, const double rho) const override;
    virtual std::vector<double> get_phi_grad(const double t, const double beta) const override;
    virtual std::vector<double> get_phi_unidirectional(const double t, const double beta, const double l, const std::vector<double>& v, const double vn) const override;

    virtual inline std::unique_ptr<MeshElementFactory> get_element_info() const override{
        return std::unique_ptr<MeshElementFactory>(new MeshElementFactoryImpl<TRI3>());
    }

    private:
    virtual std::vector<double> get_DB(const std::vector<double>& D, const gp_Pnt& point) const override;
    virtual std::vector<double> get_Nf(const double t, const std::vector<gp_Pnt>& points) const override;

    Eigen::Matrix<double, 3, 3> get_phi_radial_base(const double x, const double y, const Eigen::Vector<double, 2>& A, const Eigen::Vector<double, 2>& C, const double t, const double beta, const double vp, const double rho) const;

    double a[3], b[3], c[3], delta;

    inline double N(double x, double y, size_t i) const{
        return (a[i] + b[i]*x + c[i]*y)/(2*delta);
    }

    inline Eigen::Vector<double, 3> N_mat_1dof(double x, double y) const{
        return Eigen::Vector<double, 3>(N(x, y, 0), N(x, y, 1), N(x, y, 2));
    }
    inline Eigen::Matrix<double, 2, 3> dN_mat_1dof() const{
        return Eigen::Matrix<double, 2, 3>{{b[0], b[1], b[2]},
                                           {c[0], c[1], c[2]}}/(2*delta);
    }
};

}

#endif
