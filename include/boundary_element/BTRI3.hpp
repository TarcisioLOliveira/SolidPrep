/*
 *   Copyright (C) 2023 Tarcísio Ladeia de Oliveira.
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

#ifndef BTRI3_HPP
#define BTRI3_HPP

#include <Eigen/Core>
#include "element.hpp"
#include "material.hpp"
#include <vector>
#include "utils.hpp"
#include "element_factory.hpp"
#include "element_common.hpp"

namespace boundary_element{

class BTRI3 : public BoundaryMeshElement{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 2;
    static const size_t NODE_DOF       = 2;
    static const size_t NODES_PER_ELEM = 3;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;
    static const size_t S_SIZE         = 3; // Size of the stress and strain vectors
    static const size_t DIM            = 2; // Number of dimensions
    static const Element::Shape SHAPE_TYPE = Element::Shape::TRI;
    static const utils::ProblemType PROBLEM_TYPE = utils::PROBLEM_TYPE_3D;

    BTRI3(ElementShape s);

    virtual Eigen::MatrixXd diffusion_1dof(const Eigen::MatrixXd& A) const override;
    virtual Eigen::MatrixXd advection_1dof(const Eigen::VectorXd& v) const override;
    virtual Eigen::MatrixXd absorption_1dof() const override;
    virtual Eigen::VectorXd source_1dof() const override;

    private:
    double a[3], b[3], c[3], delta;

    inline gp_Pnt GS_point(double c1, double c2, double c3) const{
        return gp_Pnt(
            c1*this->nodes[0]->point.X() + c2*this->nodes[1]->point.X() + c3*this->nodes[2]->point.X(),
            c1*this->nodes[0]->point.Y() + c2*this->nodes[1]->point.Y() + c3*this->nodes[2]->point.Y(),
            c1*this->nodes[0]->point.Z() + c2*this->nodes[1]->point.Z() + c3*this->nodes[2]->point.Z()
        );
    }
    inline double N(const gp_Pnt& p, size_t i) const{
        return (a[i]*p.X() + b[i]*p.Y() + c[i]*p.Z())/(2*delta);
    }
    inline Eigen::Vector<double, 3> N_mat_1dof(const gp_Pnt& p) const{
        return Eigen::Vector<double, 3>(N(p, 0), N(p, 1), N(p, 2));
    }
    inline Eigen::Matrix<double, 3, 3> dN_mat_1dof() const{
        return Eigen::Matrix<double, 3, 3>{{a[0], a[1], a[2]},
                                           {b[0], b[1], b[2]},
                                           {c[0], c[1], c[2]}}/(2*delta);
    }
};

}

#endif
