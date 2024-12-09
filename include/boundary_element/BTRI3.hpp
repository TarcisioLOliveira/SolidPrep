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

#include "element.hpp"
#include "material.hpp"
#include <vector>
#include "utils.hpp"

namespace boundary_element{

class BTRI3 : public BoundaryMeshElement{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 2;
    static const size_t NODE_DOF       = 3;
    static const size_t NODES_PER_ELEM = 3;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;
    static const size_t S_SIZE         = 6; // Size of the stress and strain vectors
    static const size_t DIM            = 2; // Number of dimensions
    static const Element::Shape SHAPE_TYPE = Element::Shape::TRI;
    static const utils::ProblemType PROBLEM_TYPE = utils::PROBLEM_TYPE_3D;

    //
    // USES TRANSFORMED COORDINATES
    // (u, v, w) system used for applying Spring/InternalLoads objects
    //

    BTRI3(ElementShape s, const MeshElement* const parent);

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
            c1*this->nodes[0]->point.Z() + c2*this->nodes[1]->point.Z() + c3*this->nodes[2]->point.Z()
        );
    }
    inline double N(const gp_Pnt& p, size_t i) const{
        return a[i] + b[i]*p.X() + c[i]*p.Y();
    }
    inline math::Vector N_mat_1dof(const gp_Pnt& p) const{
        return math::Vector({N(p, 0), N(p, 1), N(p, 2)});
    }
    inline math::Matrix dN_mat_1dof() const{
        return math::Matrix({b[0], b[1], b[2],
                             c[0], c[1], c[2]}, 2, 3);
    }
    inline math::Matrix get_eps(const gp_Pnt& p, const gp_Pnt& center) const{
        const double x = p.X() - center.X();
        const double y = p.Y() - center.Y();
        return math::Matrix(
            {b[0], 0, 0, b[1], 0, 0, b[2], 0, 0, 0, 0, 0, 0, 0, 0,
             0, c[0], 0, 0, c[1], 0, 0, c[2], 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, x, y, 1, 0, 0, 0,
             0, 0, c[0], 0, 0, c[1], 0, 0, c[2], 0, 0, 0, -x, 1, 0,
             0, 0, b[0], 0, 0, b[1], 0, 0, b[2], 0, 0, 0, y, 0, 1,
             c[0], b[0], 0, c[1], b[1], 0, c[2], c[2], 0, 0, 0, 0, 0, 0, 0}, S_SIZE, K_DIM + 6);
    };
    inline math::Matrix get_deps() const{
        return math::Matrix(
            {b[0], 0, 0, b[1], 0, 0, b[2], 0, 0,
             0, c[0], 0, 0, c[1], 0, 0, c[2], 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, c[0], 0, 0, c[1], 0, 0, c[2],
             0, 0, b[0], 0, 0, b[1], 0, 0, b[2],
             c[0], b[0], 0, c[1], b[1], 0, c[2], c[2], 0}, S_SIZE, K_DIM);
    };
};

}

#endif
