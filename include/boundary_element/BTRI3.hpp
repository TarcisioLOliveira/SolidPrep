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

    //
    // USES TRANSFORMED COORDINATES
    // (u, v, w) system used for applying Spring/InternalLoads objects
    //

    BTRI3(ElementShape s, const MeshElement* const parent);

    virtual Eigen::MatrixXd diffusion_1dof(const Eigen::MatrixXd& A) const override;
    virtual Eigen::MatrixXd advection_1dof(const Eigen::VectorXd& v) const override;
    virtual Eigen::MatrixXd absorption_1dof() const override;
    virtual Eigen::VectorXd source_1dof() const override;
    virtual Eigen::VectorXd flow_1dof(const std::array<const Node*, 2>& nodes) const override;

    virtual Eigen::VectorXd grad_1dof_upos(const gp_Pnt& p, const std::vector<double>& phi) const override;
    virtual Eigen::VectorXd grad_1dof_id(const gp_Pnt& p, const std::vector<double>& phi) const override;
    virtual Eigen::VectorXd dF_2dof_id(const gp_Pnt& p, const std::vector<double>& phi) const override;
    virtual Eigen::MatrixXd int_grad_phi() const override;
    virtual Eigen::MatrixXd int_grad_phi_x(const gp_Pnt& center) const override;
    virtual Eigen::MatrixXd int_grad_phi_y(const gp_Pnt& center) const override;
    virtual Eigen::MatrixXd int_grad_xi() const override;
    virtual Eigen::MatrixXd int_grad_xi_x(const gp_Pnt& center) const override;
    virtual Eigen::MatrixXd int_grad_xi_y(const gp_Pnt& center) const override;
    virtual Eigen::MatrixXd int_grad_F() const override;
    virtual Eigen::MatrixXd int_grad_F_x(const gp_Pnt& center) const override;
    virtual Eigen::MatrixXd int_grad_F_y(const gp_Pnt& center) const override;
    virtual Eigen::VectorXd int_N_x(const gp_Pnt& center) const override;
    virtual Eigen::VectorXd int_N_y(const gp_Pnt& center) const override;
    virtual Eigen::MatrixXd int_NdN(const std::vector<double>& phi) const override;
    virtual Eigen::MatrixXd int_NdF(const std::vector<double>& phi) const override;

    virtual Eigen::MatrixXd int_grad_F_t2_t1(const Eigen::MatrixXd& B3, const gp_Pnt& center) const override;
    virtual Eigen::MatrixXd int_grad_phi_t2_t1(const Eigen::MatrixXd& B2, const gp_Pnt& center) const override;
    virtual Eigen::MatrixXd int_grad_F_D(const Eigen::MatrixXd& a, const gp_Pnt& center) const override;
    virtual Eigen::MatrixXd int_grad_phi_D(const Eigen::MatrixXd& a, const gp_Pnt& center) const override;

    virtual Eigen::MatrixXd L4(const Eigen::MatrixXd& B) const override;
    virtual Eigen::MatrixXd L3(const Eigen::MatrixXd& B) const override;
    virtual Eigen::MatrixXd L2(const Eigen::MatrixXd& B) const override;

    virtual Eigen::MatrixXd L3xi(const Eigen::MatrixXd& B) const override;
    virtual Eigen::MatrixXd L2xi(const Eigen::MatrixXd& B) const override;

    virtual Eigen::MatrixXd L4chi(const Eigen::MatrixXd& B) const override;
    virtual Eigen::MatrixXd L3Tchi(const Eigen::MatrixXd& B) const override;

    virtual Eigen::MatrixXd L4zeta(const Eigen::MatrixXd& B) const override;
    virtual Eigen::MatrixXd L3Tzeta(const Eigen::MatrixXd& B) const override;

    virtual double get_area() const override{
        return this->delta;
    }

    virtual gp_Pnt get_centroid() const override{
        const size_t N = BTRI3::NODES_PER_ELEM;

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
    double a[3], b[3], c[3], delta;

    inline gp_Pnt GS_point(double c1, double c2, double c3) const{
        return gp_Pnt(
            c1*this->nodes[0]->point.X() + c2*this->nodes[1]->point.X() + c3*this->nodes[2]->point.X(),
            c1*this->nodes[0]->point.Y() + c2*this->nodes[1]->point.Y() + c3*this->nodes[2]->point.Y(),
            0
        );
    }
    inline double N(const gp_Pnt& p, size_t i) const{
        return a[i] + b[i]*p.X() + c[i]*p.Y();
    }
    inline Eigen::Vector<double, 3> N_mat_1dof(const gp_Pnt& p) const{
        return Eigen::Vector<double, 3>(N(p, 0), N(p, 1), N(p, 2));
    }
    inline Eigen::Matrix<double, 2, 3> dN_mat_1dof() const{
        return Eigen::Matrix<double, 2, 3>{{b[0], b[1], b[2]},
                                           {c[0], c[1], c[2]}};
    }
    inline Eigen::Matrix<double, 2, 3> dxi_mat_1dof() const{
        return Eigen::Matrix<double, 2, 3>{{-c[0], -c[1], -c[2]},
                                           {b[0], b[1], b[2]}};
    }
    inline Eigen::Matrix<double, 2, 3> dxi_mat_1dof2() const{
        return Eigen::Matrix<double, 2, 3>{{c[0], c[1], c[2]},
                                           {b[0], b[1], b[2]}};
    }
    inline Eigen::Matrix<double, 3, 6> dF_mat_2dof() const{
        return Eigen::Matrix<double, 3, 6>{{b[0],    0, b[1],    0, b[2],    0},
                                           {   0, c[0],    0, c[1],    0, c[2]},
                                           {0.5*c[0], 0.5*b[0], 0.5*c[1], 0.5*b[1], 0.5*c[2], 0.5*b[2]}};
    }
    inline Eigen::Matrix<double, 3, 3> dchi_mat_1dof() const{
        return Eigen::Matrix<double, 3, 3>{
                                           {-b[0], -b[1], -b[2]},
                                           {b[0], b[1], b[2]},
                                           {-c[0], -c[1], -c[2]}};
    }
    inline Eigen::Matrix<double, 3, 3> dzeta_mat_1dof() const{
        return Eigen::Matrix<double, 3, 3>{
                                           {c[0], c[1], c[2]},
                                           {-c[0], -c[1], -c[2]},
                                           {-b[0], -b[1], -b[2]}};
    }
};

}

#endif
