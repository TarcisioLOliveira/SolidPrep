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

namespace boundary_element{

class BTRI6 : public BoundaryMeshElement{
    public:
    static const size_t ORDER          = 2;
    static const size_t GMSH_TYPE      = 9;
    static const size_t NODE_DOF       = 3;
    static const size_t NODES_PER_ELEM = 6;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;
    static const size_t S_SIZE         = 6; // Size of the stress and strain vectors
    static const size_t DIM            = 2; // Number of dimensions
    static const Element::Shape SHAPE_TYPE = Element::Shape::TRI;
    static const utils::ProblemType PROBLEM_TYPE = utils::PROBLEM_TYPE_3D;

    //
    // USES TRANSFORMED COORDINATES
    // (u, v, w) system used for applying Spring/InternalLoads objects
    //

    BTRI6(ElementShape s, const MeshElement* const parent);

    virtual std::vector<double> get_K_ext(const Eigen::MatrixXd& D, const gp_Pnt& center) const override;
    virtual std::vector<double> get_normal_stresses(const Eigen::MatrixXd& D, const std::vector<double>& u, const gp_Pnt& p, const gp_Pnt& center) const override;
    virtual std::vector<double> get_stress_integrals(const Eigen::MatrixXd& D, const gp_Pnt& center) const override;
    virtual std::vector<double> get_equilibrium_partial(const Eigen::MatrixXd& D, const gp_Pnt& center, const std::vector<size_t>& stresses) const override;
    virtual std::vector<double> get_dz_vector(const Eigen::MatrixXd& S, const Eigen::MatrixXd& D, const double Az, const double Bz, const gp_Pnt& center) const override;
    virtual std::vector<double> get_force_vector(const Eigen::MatrixXd& D, const std::vector<double>& u, const gp_Pnt& center, const Eigen::MatrixXd& rot) const override;

    virtual Eigen::MatrixXd diffusion_1dof(const Eigen::MatrixXd& A) const override;
    virtual Eigen::MatrixXd advection_1dof(const Eigen::VectorXd& v) const override;
    virtual Eigen::MatrixXd absorption_1dof() const override;
    virtual Eigen::VectorXd source_1dof() const override;
    virtual Eigen::VectorXd flow_1dof(const std::array<const Node*, 2>& nodes) const override;

    virtual double get_area() const override{
        return this->delta;
    }

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
            c1*this->nodes[0]->point.X() + c2*this->nodes[1]->point.X() + c3*this->nodes[2]->point.X(),
            c1*this->nodes[0]->point.Y() + c2*this->nodes[1]->point.Y() + c3*this->nodes[2]->point.Y(),
            c1*this->nodes[0]->point.Z() + c2*this->nodes[1]->point.Z() + c3*this->nodes[2]->point.Z()
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

    inline Eigen::Matrix<double, S_SIZE, K_DIM + 6> get_eps(const gp_Pnt& p, const gp_Pnt& center) const{
        const auto dN = this->dN_mat_1dof(p);
        const double x = p.X() - center.X();
        const double y = p.Y() - center.Y();
        Eigen::Matrix<double, S_SIZE, K_DIM + 6> M;
        M.fill(0);

        /*
            {{dN(0,i),       0,       0, 0, 0, 0,  0, 0, 0},
             {      0, dN(1,i),       0, 0, 0, 0,  0, 0, 0},
             {      0,       0,       0, x, y, 1,  0, 0, 0},
             {      0,       0, dN(1,i), 0, 0, 0, -x, 1, 0},
             {      0,       0, dN(0,i), 0, 0, 0,  y, 0, 1},
             {dN(1,i), dN(0,i),       0, 0, 0, 0,  0, 0, 0}};
        */

        M(2, K_DIM + 0) = x;
        M(2, K_DIM + 1) = y;
        M(2, K_DIM + 2) = 1;

        M(3, K_DIM + 3) = -x;
        M(3, K_DIM + 4) = 1;

        M(4, K_DIM + 3) = y;
        M(4, K_DIM + 5) = 1;

        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            M(0, i*NODE_DOF + 0) = dN(0, i);
            M(1, i*NODE_DOF + 1) = dN(1, i);

            M(3, i*NODE_DOF + 2) = dN(1, i);
            M(4, i*NODE_DOF + 2) = dN(0, i);

            M(5, i*NODE_DOF + 0) = dN(1, i);
            M(5, i*NODE_DOF + 1) = dN(0, i);
        }

        return M;

    };
    inline Eigen::Matrix<double, S_SIZE, K_DIM> get_deps(const gp_Pnt& p) const{
        const auto dN = this->dN_mat_1dof(p);
        Eigen::Matrix<double, S_SIZE, K_DIM> M;
        M.fill(0);
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            M(0, i*NODE_DOF + 0) = dN(0, i);
            M(1, i*NODE_DOF + 1) = dN(1, i);

            M(3, i*NODE_DOF + 2) = dN(1, i);
            M(4, i*NODE_DOF + 2) = dN(0, i);

            M(5, i*NODE_DOF + 0) = dN(1, i);
            M(5, i*NODE_DOF + 1) = dN(0, i);
        }

        return M;
    };
};

}

#endif
