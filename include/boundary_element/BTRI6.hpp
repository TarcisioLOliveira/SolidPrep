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

#ifndef BTRI6_HPP
#define BTRI6_HPP

#include <Eigen/Core>
#include "element.hpp"
#include "material.hpp"
#include <vector>
#include "utils.hpp"
#include "element_factory.hpp"
#include "element_common.hpp"

namespace boundary_element{

class BTRI6 : public BoundaryMeshElement{
    public:
    static const size_t ORDER          = 2;
    static const size_t GMSH_TYPE      = 9;
    static const size_t NODE_DOF       = 2;
    static const size_t NODES_PER_ELEM = 6;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;
    static const size_t S_SIZE         = 3; // Size of the stress and strain vectors
    static const size_t DIM            = 2; // Number of dimensions
    static const Element::Shape SHAPE_TYPE = Element::Shape::TRI;
    static const utils::ProblemType PROBLEM_TYPE = utils::PROBLEM_TYPE_3D;

    //
    // USES TRANSFORMED COORDINATES
    // (u, v, w) system used for applying Spring/InternalLoads objects
    //

    BTRI6(ElementShape s, const MeshElement* const parent);

    virtual Eigen::MatrixXd diffusion_1dof(const Eigen::MatrixXd& A) const override;
    virtual Eigen::MatrixXd advection_1dof(const Eigen::VectorXd& v) const override;
    virtual Eigen::MatrixXd absorption_1dof() const override;
    virtual Eigen::VectorXd source_1dof() const override;
    virtual Eigen::VectorXd source_1dof(const Eigen::Vector<double, 3>& v) const override;
    virtual Eigen::VectorXd source_grad_1dof(const Eigen::VectorXd& v) const override;
    virtual std::array<Eigen::VectorXd, 2> source_1dof(const double S13, const double Gxy, const double S12, const double Gxz, const gp_Pnt& center) const override;
    virtual Eigen::VectorXd source_1dof(double dcurv_v, double dcurv_w, const gp_Pnt& center) const override;
    virtual Eigen::VectorXd flow_1dof(double dcurv_v, double dcurv_w, const gp_Pnt& center, const gp_Dir& n, const std::vector<gp_Pnt>& edges) const override;

    virtual Eigen::VectorXd grad_1dof_upos(const gp_Pnt& p, const std::vector<double>& phi) const override;
    virtual Eigen::VectorXd grad_1dof_id(const gp_Pnt& p, const std::vector<double>& phi) const override;
    virtual Eigen::MatrixXd int_grad_1dof() const override;

    virtual Eigen::MatrixXd L4(const Eigen::MatrixXd& B) const override;
    virtual Eigen::MatrixXd L3(const Eigen::MatrixXd& B) const override;
    virtual Eigen::MatrixXd L2(const Eigen::MatrixXd& B) const override;

    virtual gp_Pnt get_centroid() const override{
        const size_t N = BTRI6::NODES_PER_ELEM;

        double x = 0;
        double y = 0;
        double z = 0;
        for(size_t i = 0; i < N; ++i){
            const auto& n = this->nodes[i];
            x += n->point.X();
            y += n->point.Y();
            z += n->point.Z();
        }

        return gp_Pnt(x/N, y/N, z/N);
    }
    virtual gp_Dir get_normal() const override{
        gp_Vec v1(this->nodes[0]->point, this->nodes[1]->point);
        gp_Vec v2(this->nodes[0]->point, this->nodes[2]->point);

        return v1.Crossed(v2);
    }

    private:
    double a[6], b[6], c[6], d[6], e[6], f[6], delta;

    inline gp_Pnt GS_point(double c1, double c2, double c3) const{
        return gp_Pnt(
            0,
            c1*this->nodes[0]->point.X() + c2*this->nodes[1]->point.X() + c3*this->nodes[2]->point.X(),
            c1*this->nodes[0]->point.Y() + c2*this->nodes[1]->point.Y() + c3*this->nodes[2]->point.Y()
        );
    }
    inline double N(const gp_Pnt& p, size_t i) const{
        return a[i] + b[i]*p.X() + c[i]*p.Y() + d[i]*p.X()*p.Y() + e[i]*p.X()*p.X() + f[i]*p.Y()*p.Y();
    }
    inline double dNdx(const gp_Pnt& p, size_t i) const{
        return b[i] + d[i]*p.Y() + 2*e[i]*p.X();
    }
    inline double dNdy(const gp_Pnt& p, size_t i) const{
        return c[i] + d[i]*p.X() + 2*f[i]*p.Y();
    }
    inline Eigen::Vector<double, 6> N_mat_1dof(const gp_Pnt& p) const{
        return Eigen::Vector<double, 6>(N(p, 0), N(p, 1), N(p, 2), N(p, 3), N(p, 4), N(p, 5));
    }
    inline Eigen::Matrix<double, 2, 6> dN_mat_1dof(const gp_Pnt& p) const{
        return Eigen::Matrix<double, 2, 6>{{dNdx(p, 0), dNdx(p, 1), dNdx(p, 2), dNdx(p, 3), dNdx(p, 4), dNdx(p, 5)},
                                           {dNdy(p, 0), dNdy(p, 1), dNdy(p, 2), dNdy(p, 3), dNdy(p, 4), dNdy(p, 5)}};
    }

    inline Eigen::Matrix<double, 3, 6> ddN_mat_1dof() const{
        return Eigen::Matrix<double, 3, 6>{{2*e[0], 2*e[1], 2*e[2], 2*e[3], 2*e[4], 2*e[5]},
                                           {2*f[0], 2*f[1], 2*f[2], 2*f[3], 2*f[4], 2*f[5]},
                                           {d[0], d[1], d[2], d[3], d[4], d[5]}};
    }
};

}

#endif
