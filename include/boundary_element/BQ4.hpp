/*
 *   Copyright (C) 2024 Tarcísio Ladeia de Oliveira.
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

#ifndef BQ4_HPP
#define BQ4_HPP

#include "element.hpp"
#include "material.hpp"
#include <vector>
#include "utils.hpp"

namespace boundary_element{

class BQ4 : public BoundaryMeshElement{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 3;
    static const size_t NODE_DOF       = 3;
    static const size_t NODES_PER_ELEM = 4;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;
    static const size_t S_SIZE         = 6; // Size of the stress and strain vectors
    static const size_t DIM            = 2; // Number of dimensions
    static const Element::Shape SHAPE_TYPE = Element::Shape::QUAD;
    static const utils::ProblemType PROBLEM_TYPE = utils::PROBLEM_TYPE_3D;

    //
    // USES TRANSFORMED COORDINATES
    // (u, v, w) system used for applying Spring/InternalLoads objects
    //

    BQ4(ElementShape s, const MeshElement* const parent);

    virtual math::Matrix get_K_ext(const math::Matrix& D, const gp_Pnt& center) const override;
    virtual math::Vector get_normal_stresses(const math::Matrix& D, const std::vector<double>& u, const gp_Pnt& p, const gp_Pnt& center) const override;
    virtual math::Matrix get_stress_integrals(const math::Matrix& D, const gp_Pnt& center) const override;
    virtual math::Matrix get_equilibrium_partial(const math::Matrix& D, const gp_Pnt& center, const std::vector<size_t>& stresses) const override;
    virtual math::Vector get_dz_vector(const math::Matrix& S, const math::Matrix& D, const double Az, const double Bz, const gp_Pnt& center) const override;
    virtual math::Vector get_force_vector(const math::Matrix& D, const std::vector<double>& u, const gp_Pnt& center, const math::Matrix& rot) const override;

    virtual math::Matrix diffusion_1dof(const math::Matrix& A) const override;
    virtual math::Matrix advection_1dof(const math::Vector& v) const override;
    virtual math::Matrix absorption_1dof() const override;
    virtual math::Vector source_1dof() const override;
    virtual math::Vector flow_1dof(const std::array<const Node*, 2>& nodes) const override;

    virtual double get_area() const override{
        return this->A;
    }

    virtual gp_Pnt get_centroid() const override{
        const size_t N = BQ4::NODES_PER_ELEM;

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
    double a[NODES_PER_ELEM], b[NODES_PER_ELEM], c[NODES_PER_ELEM], d[NODES_PER_ELEM], A;

    inline double N(double x, double y, size_t i) const{
        return a[i] + b[i]*x + c[i]*y + d[i]*x*y;
    }
    inline double N_norm(double x, double y, size_t i) const{
        switch(i){
            case 0:
                return 0.25*(1-x)*(1-y);
            case 1:
                return 0.25*(1+x)*(1-y);
            case 2:
                return 0.25*(1+x)*(1+y);
            case 3:
                return 0.25*(1-x)*(1+y);
        }
        return 0;
    }

    inline math::Vector N_mat_1dof(const gp_Pnt& p) const{
        const double x = p.X();
        const double y = p.Y();
        return math::Vector
                {N(x,y,0), N(x,y,1), N(x,y,2), N(x,y,3)};
    }

    inline double dNdx(double x, double y, size_t i) const{
        (void)x;
        return b[i] + d[i]*y;
    }
    inline double dNdy(double x, double y, size_t i) const{
        (void)y;
        return c[i] + d[i]*x;
    }
    inline math::Matrix dN_mat_1dof(const gp_Pnt& p) const{
        const double x = p.X();
        const double y = p.Y();
        return math::Matrix(
                {dNdx(x,y,0),dNdx(x,y,1),dNdx(x,y,2),dNdx(x,y,3),
                 dNdy(x,y,0),dNdy(x,y,1),dNdy(x,y,2),dNdy(x,y,3)}, DIM, NODES_PER_ELEM);
    }
    inline gp_Pnt norm_to_nat(double xi, double eta) const{
        double X = 0, Y = 0, Z = 0;
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            const double Ni = N_norm(xi, eta, i);
            const gp_Pnt& pi = this->nodes[i]->point;
            X += Ni*pi.X();
            Y += Ni*pi.Y();
            Z += Ni*pi.Z();
        }

        return gp_Pnt(X, Y, Z);
    }
    inline math::Matrix J(const double xi, const double eta) const{
        std::array<double, NODES_PER_ELEM> x{
            this->nodes[0]->point.X(),
            this->nodes[1]->point.X(),
            this->nodes[2]->point.X(),
            this->nodes[3]->point.X()
        };
        std::array<double, NODES_PER_ELEM> y{
            this->nodes[0]->point.Y(),
            this->nodes[1]->point.Y(),
            this->nodes[2]->point.Y(),
            this->nodes[3]->point.Y()
        };
        return math::Matrix(
        {
            eta*x[0]/4 - eta*x[1]/4 + eta*x[2]/4 - eta*x[3]/4 - x[0]/4 + x[1]/4 + x[2]/4 - x[3]/4
            ,
            eta*y[0]/4 - eta*y[1]/4 + eta*y[2]/4 - eta*y[3]/4 - y[0]/4 + y[1]/4 + y[2]/4 - y[3]/4
            ,
            x[0]*xi/4 - x[0]/4 - x[1]*xi/4 - x[1]/4 + x[2]*xi/4 + x[2]/4 - x[3]*xi/4 + x[3]/4
            ,
            xi*y[0]/4 - xi*y[1]/4 + xi*y[2]/4 - xi*y[3]/4 - y[0]/4 - y[1]/4 + y[2]/4 + y[3]/4
        }, DIM, DIM);
    }
    inline math::Matrix get_eps(const gp_Pnt& p, const gp_Pnt& center) const{
        const auto dN = this->dN_mat_1dof(p);
        const double x = p.X() - center.X();
        const double y = p.Y() - center.Y();
        math::Matrix M(S_SIZE, K_DIM + 6);

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
    inline math::Matrix get_deps(const gp_Pnt& p) const{
        const auto dN = this->dN_mat_1dof(p);
        math::Matrix M(S_SIZE, K_DIM);
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
