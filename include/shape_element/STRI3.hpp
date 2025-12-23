/*
 *   Copyright (C) 2025 Tarcísio Ladeia de Oliveira.
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

#ifndef STRI3_HPP
#define STRI3_HPP

#include "element.hpp"
#include "math/matrix.hpp"
#include "utils.hpp"

namespace shape_element{

class STRI3 : public ShapeMeshElement{
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

    STRI3(ElementShape s);

    virtual math::Matrix diffusion_Ndof(const math::Matrix& A) const override;
    virtual math::Matrix absorption_Ndof() const override;

    virtual double get_area() const override{
        return this->delta;
    }

    virtual gp_Pnt get_centroid() const override{
        const size_t N = STRI3::NODES_PER_ELEM;

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
    virtual const math::Matrix& get_R() const override{
        return this->R;
    }

    protected:
    double a[3], b[3], c[3], delta;
    size_t E_KW;
    math::Matrix R;
    const double FL2_OFF = 1e-5;
    const double FL2_K = 1e5;

    double h(const double x) const{
        const double dx = x - FL2_OFF;
        if(std::abs(FL2_K*dx) < 345){
            return std::log(1.0 + std::exp(FL2_K*dx))/FL2_K;
        } else if(x > 0){
            return x;
        } else {
            return 0;
        }
    }

    double dh(const double x) const{
        const double dx = x - FL2_OFF;
        if(std::abs(FL2_K*dx) < 345){
            return 1.0/(1.0 + std::exp(-FL2_K*dx));
        } else if(x > 0){
            return 1;
        } else {
            return 0;
        }
    }

    double ddh(const double x) const{
        const double dx = x - FL2_OFF;
        if(std::abs(FL2_K*dx) < 345){
            const double ekx = std::exp(-FL2_K*dx);
            const double ekx1 = 1.0 + ekx;
            return FL2_K*ekx/(ekx1*ekx1);
        } else {
            return 0;
        }
    }

    virtual gp_Dir get_d1() const{
        return gp_Dir(R(0,1), R(1,1), R(2,1));
    }
    virtual gp_Dir get_d2() const{
        return gp_Dir(R(0,0), R(1,0), R(2,0));
    }

    inline gp_Pnt R_GS_point(double c1, double c2, double c3) const{
        const size_t N = STRI3::NODES_PER_ELEM;
        std::array<double, N> x, y;
        std::fill(x.begin(), x.end(), 0);
        std::fill(y.begin(), y.end(), 0);
        for(size_t i = 0; i < N; ++i){
            for(size_t j = 0; j < NODE_DOF; ++j){
                x[i] += R(j,0)*this->nodes[i]->point.Coord(1+j);
                y[i] += R(j,1)*this->nodes[i]->point.Coord(1+j);
            }
        }
        return gp_Pnt(
            c1*x[0] + c2*x[1] + c3*x[2],
            c1*y[0] + c2*y[1] + c3*y[2],
            0
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
    inline math::Vector N_mat_1dof(const gp_Pnt& p) const{
        return math::Vector{N(p, 0), N(p, 1), N(p, 2)};
    }
    inline math::Matrix N_mat_3dof(const gp_Pnt& p) const{
        math::Matrix N3(DIM, DIM*NODES_PER_ELEM);
        const auto Ni = this->N_mat_1dof(p);
        for(size_t i = 0; i < DIM; ++i){
            for(size_t j = 0; j < NODES_PER_ELEM; ++j){
                N3(i, j*DIM + i) = Ni[j];
            }
        }

        return N3;
    }
    inline math::Matrix dN_mat_1dof_base() const{
        return math::Matrix({b[0], b[1], b[2],
                             c[0], c[1], c[2]}, 2, 3);
    }
    inline double dNidXj(const size_t i, const size_t j) const{
        const math::Vector dNi({b[i], c[i]});
        const math::Matrix Rn(
                {R(0,0), R(1,0), R(2,0),
                 R(0,1), R(1,1), R(2,1)}, 2, 3);
        math::Vector v(DIM);
        v[j] = 1;

        return dNi.T()*(Rn*v);
    }
    inline math::Matrix dN_mat_1dof() const{
        math::Matrix gN(DIM, NODES_PER_ELEM, 0);
        for(size_t i = 0; i < DIM; ++i){
            for(size_t j = 0; j < NODES_PER_ELEM; ++j){
                gN(i,j) = dNidXj(j,i);
            }
        }

        return gN;
    }
    inline math::Matrix dN_mat_dim_dof() const{
        const auto gN = this->dN_mat_1dof();
        math::Matrix dN(DIM*DIM, DIM*NODES_PER_ELEM, 0);
        for(size_t i = 0; i < DIM; ++i){
            for(size_t j = 0; j < NODES_PER_ELEM; ++j){
                dN(DIM*i + 0, j*DIM + i) = gN(0,j);
                dN(DIM*i + 1, j*DIM + i) = gN(1,j);
                dN(DIM*i + 2, j*DIM + i) = gN(2,j);
            }
        }
        return dN;
    }
};

}

#endif
