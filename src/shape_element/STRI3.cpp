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

#include "shape_element/STRI3.hpp"
#include "math/matrix.hpp"
#include "utils/gauss_legendre.hpp"
#include "meshing.hpp"

namespace shape_element{

STRI3::STRI3(ElementShape s)
    :ShapeMeshElement(s.nodes){

    const size_t N = STRI3::NODES_PER_ELEM;

    const gp_Dir n(s.normal);
    const gp_Dir d1 = gp_Vec(this->nodes[0]->point, this->nodes[1]->point);
    const gp_Dir d2 = -d1.Crossed(n);

    this->R = math::Matrix(
        {d2.X(), d1.X(), n.X(),
         d2.Y(), d1.Y(), n.Y(),
         d2.Z(), d1.Z(), n.Z()}, 3, 3);

    std::array<double, N> x, y, z;
    std::fill(x.begin(), x.end(), 0);
    std::fill(y.begin(), y.end(), 0);
    std::fill(z.begin(), z.end(), 0);
    for(size_t i = 0; i < N; ++i){
        for(size_t j = 0; j < NODE_DOF; ++j){
            x[i] += R(j,0)*this->nodes[i]->point.Coord(1+j);
            y[i] += R(j,1)*this->nodes[i]->point.Coord(1+j);
            z[i] += R(j,2)*this->nodes[i]->point.Coord(1+j);
        }
    }

    math::Matrix M(
        {1, x[0], y[0],
         1, x[1], y[1],
         1, x[2], y[2]}, N, N);

    // M*C = I -> C=M^-1
    // C = {a[0], a[1], a[2],
    //      b[0], b[1], b[2],
    //      c[0], c[1], c[2]}
    M.invert_LU();

    for(size_t i = 0; i < N; ++i){
        this->a[i] = M(0,i);
        this->b[i] = M(1,i);
        this->c[i] = M(2,i);
    }

    gp_Vec v1(this->nodes[0]->point, this->nodes[1]->point);
    gp_Vec v2(this->nodes[0]->point, this->nodes[2]->point);

    this->delta = 0.5*(v1.Crossed(v2)).Magnitude();
}

math::Vector STRI3::shape_flow(const BoundaryElement* e, const math::Vector& f) const{
    const auto& gli = utils::GaussLegendreTri<2*ORDER>::get();
    const auto elem_info = e->parent->get_element_info();
    const size_t NUM_NODES = elem_info->get_nodes_per_element();
    math::Vector V(NUM_NODES);

    const gp_Dir n = this->get_normal();
    const math::Vector vn((e->normal.Dot(n) > 0) 
            ? 
                math::Vector{n.X(), n.Y(), n.Z()}
            : 
                math::Vector{-n.X(), -n.Y(), -n.Z()});
    math::Vector fi(3);
    for(auto it = gli.begin(); it != gli.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const gp_Pnt rpi = this->R_GS_point(it->a, it->b, it->c);
        const auto Np = e->parent->get_Ni_1dof(pi);
        const auto Nb = this->N_mat_1dof(rpi);
        fi.fill(0);
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            for(size_t j = 0; j < DIM; ++j){
                fi[j] += Nb[i]*f[i*DIM + j];
            }
        }
        V += (it->w*(fi.T()*vn))*(Np);
    }
    V *= this->delta;

    return V;
}

}
