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

#ifndef CTRI3_HPP
#define CTRI3_HPP

#include "element.hpp"
#include "material.hpp"
#include <vector>
#include "math/matrix.hpp"
#include "utils.hpp"

namespace contact_element{

class CTRI3 : public ContactMeshElement{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 2;
    static const size_t NODE_DOF       = 3;
    static const size_t NODES_PER_ELEM = 3;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;
    static const size_t S_SIZE         = 6; // Size of the stress and strain vectors
    static const size_t DIM            = 3; // Number of dimensions
    static const Element::Shape SHAPE_TYPE = Element::Shape::TRI;
    static const utils::ProblemType PROBLEM_TYPE = utils::PROBLEM_TYPE_3D;

    CTRI3(ElementShape s, const MeshElement* const e1, const MeshElement* const e2);

    virtual math::Matrix get_frictionless_Ge() const override;

    virtual math::Matrix fl_uLne(const math::Matrix& D, const math::Vector& ln_e) const override;
    virtual math::Matrix fl_LnLne(const math::Matrix& D, const math::Vector& ln_e) const override;
    virtual math::Matrix fl_LtLne(const math::Matrix& D, const math::Vector& ln_e, size_t ti) const override;
    virtual math::Matrix fl_uLte(const math::Matrix& D, size_t ti) const override;
    virtual math::Matrix fl_LtLte(const math::Matrix& D, size_t t1, size_t t2) const override;

    virtual math::Matrix fl2_uL(const math::Vector& l_e) const override;
    virtual math::Matrix fl2_LL(const math::Vector& l_e, const math::Vector& u1, const math::Vector& u2) const override;

    virtual void fl2_Ku_lambda(const double EPS, const std::vector<long> u1_pos, const std::vector<long> u2_pos, const std::vector<long>& lu_pos, const std::vector<double>& u, std::vector<double>& Ku) const override;

    virtual double get_area() const override{
        return this->delta;
    }

    virtual gp_Pnt get_centroid() const override{
        const size_t N = CTRI3::NODES_PER_ELEM;

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
        return gp_Dir(R(0,2), R(1,2), R(2,2));
    }

    protected:
    double a[3], b[3], c[3], delta;
    size_t E_KW;
    math::Matrix R;

    virtual gp_Dir get_d1() const{
        return gp_Dir(R(0,1), R(1,1), R(2,1));
    }
    virtual gp_Dir get_d2() const{
        return gp_Dir(R(0,0), R(1,0), R(2,0));
    }

    inline gp_Pnt R_GS_point(double c1, double c2, double c3) const{
        const size_t N = CTRI3::NODES_PER_ELEM;
        std::array<double, N> x, y;//, z;
        std::fill(x.begin(), x.end(), 0);
        std::fill(y.begin(), y.end(), 0);
        //std::fill(z.begin(), z.end(), 0);
        for(size_t i = 0; i < N; ++i){
            for(size_t j = 0; j < NODE_DOF; ++j){
                x[i] += R(0,j)*this->nodes[i]->point.X();
                y[i] += R(1,j)*this->nodes[i]->point.Y();
                //z[i] += R(2,j)*this->nodes[i]->point.Z();
            }
        }
        return gp_Pnt(
            c1*x[0] + c2*x[1] + c3*x[2],
            c1*y[0] + c2*y[1] + c3*y[2],
            0 //c1*z[0] + c2*z[1] + c3*z[2]
        );
    }
    inline gp_Pnt GS_point(double c1, double c2, double c3) const{
        return gp_Pnt(
            c1*this->nodes[0]->point.X() + c2*this->nodes[1]->point.X() + c3*this->nodes[2]->point.X(),
            c1*this->nodes[0]->point.Y() + c2*this->nodes[1]->point.Y() + c3*this->nodes[2]->point.Y(),
            c1*this->nodes[0]->point.Z() + c2*this->nodes[1]->point.Z() + c3*this->nodes[2]->point.Z()
        );
    }
    inline double N(const gp_Pnt& p, size_t i) const{
        return a[i] + b[i]*p.X() + c[i]*p.Y();
    }
    inline double dNdx(size_t i) const{
        return b[i];
    }
    inline double dNdy(size_t i) const{
        return c[i];
    }
    inline math::Vector N_mat_1dof(const gp_Pnt& p) const{
        return math::Vector{N(p, 0), N(p, 1), N(p, 2)};
    }
    inline math::Vector dNdx_mat_1dof() const{
        return math::Vector{dNdx(0), dNdx(1), dNdx(2)};
    }
    inline math::Vector dNdy_mat_1dof() const{
        return math::Vector{dNdy(0), dNdy(1), dNdy(2)};
    }
    inline math::Matrix gradN_1dof() const{
        return math::Matrix(
            {dNdx(0), dNdx(1), dNdx(2),
             dNdy(0), dNdy(1), dNdy(2)}, 2, 3);
    }
    inline math::Vector gradN_1dof_vec(const math::Vector& dXI) const{
        return math::Vector
            {dNdx(0)*dXI[0] + dNdy(0)*dXI[1],
             dNdx(1)*dXI[0] + dNdy(1)*dXI[1],
             dNdx(2)*dXI[0] + dNdy(2)*dXI[1]};
    }
    inline math::Matrix eps_mat_lambda(const gp_Dir& n) const
    {
        math::Matrix eps(S_SIZE, NODES_PER_ELEM);
        std::vector<double> didj(NODE_DOF*NODE_DOF*NODES_PER_ELEM);

        math::Vector dXIdx{R(0,0), R(1,0), R(2,0)};
        math::Vector dXIdy{R(0,1), R(1,1), R(2,1)};
        math::Vector dXIdz{R(0,2), R(1,2), R(2,2)};

        std::vector<math::Vector> gradN
        {
            this->gradN_1dof_vec(dXIdx),
            this->gradN_1dof_vec(dXIdy),
            this->gradN_1dof_vec(dXIdz)
        };
        for(size_t i = 0; i < NODE_DOF; ++i){
            for(size_t j = 0; j < NODE_DOF; ++j){
                const math::Vector m(n.Coord(1+i)*gradN[j]);
                for(size_t k = 0; k < NODES_PER_ELEM; ++k){
                    didj[i*NODE_DOF*NODES_PER_ELEM + j*NODES_PER_ELEM + k] = m[k];
                }
            }
        }

        for(size_t i = 0; i < NODE_DOF; ++i){
            for(size_t j = 0; j < NODES_PER_ELEM; ++j){
                eps(i, j) = didj[i*NODE_DOF*NODES_PER_ELEM + i*NODES_PER_ELEM + j];
            }
        }
        for(size_t j = 0; j < NODES_PER_ELEM; ++j){
            eps(3, j) = 
                didj[0*NODE_DOF*NODES_PER_ELEM + 1*NODES_PER_ELEM + j] +
                didj[1*NODE_DOF*NODES_PER_ELEM + 0*NODES_PER_ELEM + j];
        }
        for(size_t j = 0; j < NODES_PER_ELEM; ++j){
            eps(4, j) = 
                didj[0*NODE_DOF*NODES_PER_ELEM + 2*NODES_PER_ELEM + j] +
                didj[2*NODE_DOF*NODES_PER_ELEM + 0*NODES_PER_ELEM + j];
        }
        for(size_t j = 0; j < NODES_PER_ELEM; ++j){
            eps(5, j) = 
                didj[2*NODE_DOF*NODES_PER_ELEM + 1*NODES_PER_ELEM + j] +
                didj[1*NODE_DOF*NODES_PER_ELEM + 2*NODES_PER_ELEM + j];
        }

        return eps;
    }
};

}

#endif
