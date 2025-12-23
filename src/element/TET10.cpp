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

#include "element/TET10.hpp"
#include "boundary_element/BTRI6.hpp"
#include "logger.hpp"
#include "math/matrix.hpp"
#include "utils/gauss_legendre.hpp"
#include "project_specification/registry.hpp"

namespace element{

const bool TET10::reg = projspec::ElementRegistry::add("TET10",
    std::make_unique<MeshElementFactoryImpl<TET10>>());

TET10::TET10(ElementShape s):
    MeshElementCommon3DTet<TET10>(s.nodes){

    this->get_coeffs();    
}

std::unique_ptr<BoundaryMeshElementFactory> TET10::get_boundary_element_info() {
    return std::unique_ptr<BoundaryMeshElementFactory>(new BoundaryMeshElementFactoryImpl<boundary_element::BTRI6>());
}
std::unique_ptr<ContactMeshElementFactory> TET10::get_contact_element_info() {
    logger::log_assert(false, logger::ERROR, "CTRI6 ELEMENT TYPE NOT IMPLEMENTED");
    return std::unique_ptr<ContactMeshElementFactory>();
    //return std::unique_ptr<ContactMeshElementFactory>(new ContactMeshElementFactoryImpl<contact_element::CTRI6>());
}
std::unique_ptr<ShapeMeshElementFactory> TET10::get_shape_element_info() {
    logger::log_assert(false, logger::ERROR, "STRI6 ELEMENT TYPE NOT IMPLEMENTED");
    return std::unique_ptr<ShapeMeshElementFactory>();
    //return std::unique_ptr<ShapeMeshElementFactory>(new ShapeMeshElementFactoryImpl<shape_element::STRI6>());
}

void TET10::get_coeffs(){
    constexpr size_t N = 4;
    std::array<double, N> x, y, z;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
        z[i] = this->nodes[i]->point.Z();
    }

    math::Matrix M(
        {1, x[0], y[0], z[0],
         1, x[1], y[1], z[1],
         1, x[2], y[2], z[2],
         1, x[3], y[3], z[3]}, N, N);

    // M*C = I -> C=M^-1
    // C = {a[0], a[1], a[2], a[3],
    //      b[0], b[1], b[2], b[3],
    //      c[0], c[1], c[2], c[3],
    //      d[0], d[1], d[2], d[3]}
    M.invert_LU();

    this->V = this->get_volume(1.0);
    M *= 6*V;

    std::copy(M.data(), M.data()+4, a);
    std::copy(M.data()+4, M.data()+8, b);
    std::copy(M.data()+8, M.data()+12, c);
    std::copy(M.data()+12, M.data()+16, d);
}

math::Matrix TET10::get_k(const math::Matrix& D, const double t) const{
    (void)t;
    
    math::Matrix K(K_DIM, K_DIM);

    const auto& gli = utils::GaussLegendreTet<2>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GL_point(it->a, it->b, it->c, it->d);
        const auto B(this->B_mat(p));
        K += it->w*B.T()*D*B;
    }
    K *= this->V;

    return K;
}
math::Matrix TET10::get_B(const gp_Pnt& point) const{
    return this->B_mat(point);
}

math::Matrix TET10::get_Nf(const double t, const std::vector<gp_Pnt>& points) const{
    (void)t;

    const gp_Vec v1(points[0], points[1]);
    const gp_Vec v2(points[0], points[2]);

    const double delta = v1.Crossed(v2).Magnitude()/2;
    math::Matrix Nf(NODE_DOF, K_DIM);
    const auto& gli = utils::GaussLegendreTri<3>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GL_point_2D(it->a, it->b, it->c, points);
        const auto NN(this->N_mat(p));
        Nf += it->w*NN;
    }
    Nf *= delta;

    return Nf;
}

math::Matrix TET10::get_nodal_density_gradient(gp_Pnt p) const{
    return this->dN_mat_1dof(p);
}

math::Matrix TET10::get_R(const math::Matrix & K, const double t, const std::vector<gp_Pnt>& points) const{
    (void)t;
    math::Matrix R(K_DIM, K_DIM);

    const auto& p = points;
    gp_Vec v1(p[1], p[0]);
    gp_Vec v2(p[2], p[0]);
    const double delta = v1.Crossed(v2).Magnitude()/2;
    const auto& gli = utils::GaussLegendreTri<4>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GL_point_2D(it->a, it->b, it->c, points);
        const auto NN(this->N_mat(p));
        R += it->w*NN.T()*K*NN;
    }
    R *= delta;

    return R;
}

math::Matrix TET10::diffusion_1dof(const double t, const math::Matrix& A) const{
    (void)t;

    math::Matrix M(NODES_PER_ELEM, NODES_PER_ELEM);

    const auto& gli = utils::GaussLegendreTet<2>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GL_point(it->a, it->b, it->c, it->d);
        const auto dN(this->dN_mat_1dof(p));
        M += it->w*dN.T()*A*dN;
    }
    M *= this->V;

    return M;
}
math::Matrix TET10::advection_1dof(const double t, const math::Vector& v) const{
    (void)t;

    math::Matrix M(NODES_PER_ELEM, NODES_PER_ELEM);

    const auto& gli = utils::GaussLegendreTet<3>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GL_point(it->a, it->b, it->c, it->d);
        const auto dN(this->dN_mat_1dof(p));
        const auto N(this->N_mat_1dof(p));
        M += it->w*N*(v.T()*dN);
    }
    M *= this->V;

    return M;
}
math::Matrix TET10::absorption_1dof(const double t) const{
    (void)t;

    math::Matrix M(NODES_PER_ELEM, NODES_PER_ELEM);

    const auto& gli = utils::GaussLegendreTet<4>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GL_point(it->a, it->b, it->c, it->d);
        const auto N(this->N_mat_1dof(p));
        M += it->w*N*N.T();;
    }
    M *= this->V;

    return M;
}
math::Vector TET10::source_1dof(const double t) const{
    (void)t;

    math::Vector M(NODES_PER_ELEM);

    const auto& gli = utils::GaussLegendreTet<2>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GL_point(it->a, it->b, it->c, it->d);
        const auto N(this->N_mat_1dof(p));
        M += it->w*N;
    }
    M *= this->V;

    return M;
}

math::Vector TET10::flow_1dof(const double t, const MeshNode** nodes) const{
    (void)t;

    const gp_Vec v1(nodes[0]->point, nodes[1]->point);
    const gp_Vec v2(nodes[0]->point, nodes[2]->point);

    const double delta = v1.Crossed(v2).Magnitude()/2;
    math::Vector M(NODES_PER_ELEM);

    const std::vector<gp_Pnt> points{nodes[0]->point, nodes[1]->point, nodes[2]->point};

    const auto& gli = utils::GaussLegendreTri<2>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GL_point_2D(it->a, it->b, it->c, points);
        const auto N(this->N_mat_1dof(p));
        M += it->w*N;
    }
    M *= delta;

    return M;
}

math::Matrix TET10::get_Ni(const gp_Pnt& p) const{
    return this->N_mat(p);
}

}
