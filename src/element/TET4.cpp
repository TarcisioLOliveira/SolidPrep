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

#include "element/TET4.hpp"
#include "boundary_element/BTRI3.hpp"
#include "contact_element/CTRI3.hpp"
#include "logger.hpp"
#include "math/matrix.hpp"
#include "utils/gauss_legendre.hpp"
#include "project_specification/registry.hpp"

namespace element{

const bool TET4::reg = projspec::ElementRegistry::add("TET4",
    std::make_unique<MeshElementFactoryImpl<TET4>>());

TET4::TET4(ElementShape s):
    MeshElementCommon3DTet<TET4>(s.nodes), C(NODES_PER_ELEM, NODES_PER_ELEM){

    this->calculate_coefficients();    
}

std::unique_ptr<BoundaryMeshElementFactory> TET4::get_boundary_element_info() {
    return std::unique_ptr<BoundaryMeshElementFactory>(new BoundaryMeshElementFactoryImpl<boundary_element::BTRI3>());
}

std::unique_ptr<ContactMeshElementFactory> TET4::get_contact_element_info() {
    return std::unique_ptr<ContactMeshElementFactory>(new ContactMeshElementFactoryImpl<contact_element::CTRI3>());
}

void TET4::calculate_coefficients(){
    constexpr size_t N = TET4::NODES_PER_ELEM;

    for(size_t i = 0; i < N; ++i){
        this->C(i, 0) = 1;
        this->C(i, 1) = this->nodes[i]->point.X();
        this->C(i, 2) = this->nodes[i]->point.Y();
        this->C(i, 3) = this->nodes[i]->point.Z();
    }

    // M*C = I -> C=M^-1
    // C = {a[0], a[1], a[2], a[3],
    //      b[0], b[1], b[2], b[3],
    //      c[0], c[1], c[2], c[3],
    //      d[0], d[1], d[2], d[3]}
    C.invert_LU();

    this->V = this->get_volume(1.0);
}

math::Matrix TET4::get_k(const math::Matrix& D, const double t) const{
    (void)t;

    const auto B = this->B(this->C);

    return this->V*(B.T()*D*B);
}
math::Matrix TET4::get_B(const gp_Pnt& point) const{
    (void)point;
    return this->B(this->C);
}

math::Matrix TET4::get_Nf(const double t, const std::vector<gp_Pnt>& points) const{
    (void)t;

    const gp_Vec v1(points[0], points[1]);
    const gp_Vec v2(points[0], points[2]);

    const double AA = v1.Crossed(v2).Magnitude()/6;
    double A[4] = {0,0,0,0};
    for(size_t i = 0; i < 4; ++i){
        for(size_t j = 0; j < 3; ++j){
            if(points[j].IsEqual(this->nodes[i]->point, Precision::Confusion())){
                A[i] = AA;
                break;
            }
        }
    }

    math::Matrix Nf({
        A[0],  0,  0,
         0, A[0],  0,
         0,  0, A[0],
        A[1],  0,  0,
         0, A[1],  0,
         0,  0, A[1],
        A[2],  0,  0,
         0, A[2],  0,
         0,  0, A[2],
        A[3],  0,  0,
         0, A[3],  0,
         0,  0, A[3]
    }, K_DIM, DIM);

    return Nf;
}

math::Matrix TET4::get_nodal_density_gradient(gp_Pnt p) const{
    (void)p;

    return this->dN_mat_1dof(this->C);
}

math::Matrix TET4::get_R(const math::Matrix& K, const double t, const std::vector<gp_Pnt>& points) const{
    (void)t;
    math::Matrix R(K_DIM, K_DIM);

    auto& gl = utils::GaussLegendreTri<2>::get();

    const auto& p = points;
    gp_Vec v1(p[1], p[0]);
    gp_Vec v2(p[2], p[0]);
    const double drnorm = (v1.Crossed(v2)).Magnitude()/2;
    for(auto it = gl.begin(); it < gl.end(); ++it){
        gp_Pnt pi = this->GL_point_tri(it->a, it->b, it->c, points);
        const auto NN = N_mat(pi.X(), pi.Y(), pi.Z(), this->C);
        R += it->w*(NN.T()*K*NN);
    }
    R *= drnorm;

    return R;
}

math::Matrix TET4::diffusion_1dof(const double t, const math::Matrix& A) const{
    (void)t;
    const auto B = this->dN_mat_1dof(this->C);

    math::Matrix M = this->V*(B.T()*A*B);

    return M;
}
math::Matrix TET4::advection_1dof(const double t, const math::Vector& v) const{
    (void)t;
    math::Matrix M(NODES_PER_ELEM, NODES_PER_ELEM);

    auto& gl = utils::GaussLegendreTet<1>::get();

    for(auto it = gl.begin(); it < gl.end(); ++it){
        const gp_Pnt p = this->GL_point(it->a, it->b, it->c, it->d);
        const auto B = this->dN_mat_1dof(this->C);
        const auto N = this->N_mat_1dof(p.X(), p.Y(), p.Z(), this->C);
        M += it->w*(B.T()*v*N.T());
    }

    return this->V*M;
}
math::Matrix TET4::absorption_1dof(const double t) const{
    (void)t;
    math::Matrix M(NODES_PER_ELEM, NODES_PER_ELEM);

    auto& gl = utils::GaussLegendreTet<1>::get();

    for(auto it = gl.begin(); it < gl.end(); ++it){
        const gp_Pnt p = this->GL_point(it->a, it->b, it->c, it->d);
        const auto N = this->N_mat_1dof(p.X(), p.Y(), p.Z(), this->C);
        M += it->w*(N*N.T());
    }

    return this->V*M;
}

math::Matrix TET4::robin_1dof(const double t, const std::vector<gp_Pnt>& points) const{

    (void)t;
    math::Matrix M(NODES_PER_ELEM, NODES_PER_ELEM);

    auto& gl = utils::GaussLegendreTri<2>::get();

    const auto& p = points;
    gp_Vec v1(p[1], p[0]);
    gp_Vec v2(p[2], p[0]);
    const double drnorm = (v1.Crossed(v2)).Magnitude()/2;
    for(auto it = gl.begin(); it < gl.end(); ++it){
        gp_Pnt p = this->GL_point_tri(it->a, it->b, it->c, points);
        const auto N = this->N_mat_1dof(p.X(), p.Y(), p.Z(), this->C);
        M += it->w*(N*N.T());
    }
    M *= drnorm;

    return M;
}
math::Vector TET4::source_1dof(const double t) const{
    (void)t;

    math::Vector M{
        V/4
        ,
        V/4
        ,
        V/4
        ,
        V/4
    };
    return M;
}

math::Vector TET4::flow_1dof(const double t, const MeshNode** nodes) const{
    (void)t;

    const gp_Vec v1(nodes[0]->point, nodes[1]->point);
    const gp_Vec v2(nodes[0]->point, nodes[2]->point);

    const double AA = v1.Crossed(v2).Magnitude()/6;
    double A[4] = {0,0,0,0};
    for(size_t i = 0; i < 4; ++i){
        for(size_t j = 0; j < 3; ++j){
            if(nodes[j]->point.IsEqual(this->nodes[i]->point, Precision::Confusion())){
                A[i] = AA;
                break;
            }
        }
    }

    math::Vector M{
        A[0],
        A[1],
        A[2],
        A[3]
    };

    return M;
}

math::Matrix TET4::get_Ni(const gp_Pnt& p) const{
    return this->N_mat(p.X(), p.Y(), p.Z(), this->C);
}

math::Matrix TET4::get_dk_sh(const math::Matrix& D, const double t, const size_t n, const size_t dof) const{
    (void)t;
    const auto C2 = this->get_C_derivative(n, dof);

    const auto B1 = this->B(this->C);
    const auto B2 = this->B(C2);

    math::Matrix K(this->V*(B2.T()*D*B1 + B1.T()*D*B2));

    return K;
}

math::Matrix TET4::get_C_derivative(const size_t n, const size_t dof) const{
    math::Matrix C2(NODES_PER_ELEM, NODES_PER_ELEM);
    C2(n, (dof + 1)) = 1;

    return -1.0*(C*C2*C);
}

math::Matrix TET4::get_dB_sh(const gp_Pnt& p, const size_t n, const size_t dof) const{
    (void) p;
    const auto C2 = this->get_C_derivative(n, dof);

    return this->B(C2);
}
math::Matrix TET4::get_dN_sh(const gp_Pnt& p, const size_t n, const size_t dof) const{
    const auto C2 = this->get_C_derivative(n, dof);

    return this->N_mat(p.X(), p.Y(), p.Z(), C2);
}

math::Vector TET4::get_Ni_1dof(const gp_Pnt& p) const {
    return this->N_mat_1dof(p.X(), p.Y(), p.Z(), C);
}

}
