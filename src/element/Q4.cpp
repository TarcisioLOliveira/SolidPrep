/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
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

#include "math/matrix.hpp"
#include "utils/gauss_legendre.hpp"
#include "element/Q4.hpp"
#include "project_specification/registry.hpp"

namespace element{

const bool Q4::reg = projspec::ElementRegistry::add("Q4",
    std::make_unique<MeshElementFactoryImpl<Q4>>());

Q4::Q4(ElementShape s):
    MeshElementCommon2DQuad<Q4>(s.nodes){
    
    constexpr size_t N = Q4::NODES_PER_ELEM;
    std::array<double, N> x, y;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
    }
    math::Matrix M( 
        {1, x[0], y[0], x[0]*y[0],
         1, x[1], y[1], x[1]*y[1],
         1, x[2], y[2], x[2]*y[2],
         1, x[3], y[3], x[3]*y[3]}, N, N);

    M.invert_LU();

    for(size_t i = 0; i < N; ++i){
        this->a[i] = M.data()[i];
        this->b[i] = M.data()[i+4];
        this->c[i] = M.data()[i+8];
        this->d[i] = M.data()[i+12];
    }

    this->A = this->get_volume(1.0);
}

math::Matrix Q4::get_k(const math::Matrix& D, const double t) const{
    constexpr size_t N = Q4::NODES_PER_ELEM;
    std::array<double, N> x, y;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
    }

    math::Matrix k(K_DIM, K_DIM);
    constexpr size_t GN = 3;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto X = this->norm_to_nat(xi->x, eta->x, x, y);
            const auto B = this->B_mat_nat(X[0], X[1]);
            const double detJ = this->J(xi->x, eta->x, x, y).determinant();

            k += (xi->w*eta->w*detJ)*B.T()*D*B;
        }
    }
    k *= t;

    return k;
}

math::Matrix Q4::get_Nf(const double t, const std::vector<gp_Pnt>& points) const{
    const double x[]{points[0].X(), points[1].X()};
    const double y[]{points[0].Y(), points[1].Y()};
    const double rnorm = 0.5*points[0].Distance(points[1]);

    math::Matrix Nf(NODE_DOF, K_DIM);
    constexpr size_t GN = 1;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        const double s = xi->x;
        const double X = 0.5*(x[0]*(1-s) + x[1]*(1+s));
        const double Y = 0.5*(y[0]*(1-s) + y[1]*(1+s));
        Nf += xi->w*N_mat(X, Y);
    }
    Nf *= t*rnorm;

    return Nf;
}

math::Matrix Q4::get_R(const math::Matrix& K, const double t, const std::vector<gp_Pnt>& points) const{
    const double x[]{points[0].X(), points[1].X()};
    const double y[]{points[0].Y(), points[1].Y()};
    const double rnorm = 0.5*points[0].Distance(points[1]);

    math::Matrix R(K_DIM, K_DIM);
    constexpr size_t GN = 3;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        const double s = xi->x;
        const double X = 0.5*(x[0]*(1-s) + x[1]*(1+s));
        const double Y = 0.5*(y[0]*(1-s) + y[1]*(1+s));
        const auto NN = N_mat(X, Y);
        R += xi->w*NN.T()*K*NN;
    }
    R *= t*rnorm;

    return R;
}

math::Matrix Q4::diffusion_1dof(const double t, const math::Matrix& A) const{
    constexpr size_t N = Q4::NODES_PER_ELEM;
    std::array<double, N> x, y;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
    }

    const math::Matrix Am({A(0,0), A(0,1),
                           A(1,0), A(1,1)}, NODE_DOF, NODE_DOF);

    math::Matrix M(N, N);
    constexpr size_t GN = 3;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto X = this->norm_to_nat(xi->x, eta->x, x, y);
            const auto B = this->dN_mat_1dof(X[0], X[1]);
            const double detJ = this->J(xi->x, eta->x, x, y).determinant();

            M += (xi->w*eta->w*detJ)*B.T()*Am*B;
        }
    }
    M *= t;

    return M;
}
math::Matrix Q4::advection_1dof(const double t, const math::Vector& v) const{
    constexpr size_t N = Q4::NODES_PER_ELEM;
    std::array<double, N> x, y;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
    }

    const math::Vector vv{v[0], v[1]};

    math::Matrix  M(N, N);
    M.fill(0);
    constexpr size_t GN = 4;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto X = this->norm_to_nat(xi->x, eta->x, x, y);
            const auto B = this->dN_mat_1dof(X[0], X[1]);
            const auto Nv = this->N_mat_1dof(X[0], X[1]);
            const double detJ = this->J(xi->x, eta->x, x, y).determinant();

            M += (xi->w*eta->w*detJ)*B.T()*vv*Nv.T();
        }
    }
    M *= t;

    return M;
}
math::Matrix Q4::absorption_1dof(const double t) const{
    constexpr size_t N = Q4::NODES_PER_ELEM;
    std::array<double, N> x, y;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
    }

    math::Matrix M(N, N);
    constexpr size_t GN = 5;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto X = this->norm_to_nat(xi->x, eta->x, x, y);
            const auto Nv = this->N_mat_1dof(X[0], X[1]);
            const double detJ = this->J(xi->x, eta->x, x, y).determinant();

            M += (xi->w*eta->w*detJ)*Nv*Nv.T();
        }
    }
    M *= t;

    return M;
}
math::Vector Q4::source_1dof(const double t) const{
    constexpr size_t N = Q4::NODES_PER_ELEM;
    std::array<double, N> x, y;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
    }

    math::Vector M(N);
    constexpr size_t GN = 3;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto X = this->norm_to_nat(xi->x, eta->x, x, y);
            const auto Nv = this->N_mat_1dof(X[0], X[1]);
            const double detJ = this->J(xi->x, eta->x, x, y).determinant();

            M += (xi->w*eta->w*detJ)*Nv;
        }
    }
    M *= t;

    return M;
}

math::Vector Q4::flow_1dof(const double t, const MeshNode** nodes) const{
    std::array<double, 2> x{nodes[0]->point.X(), nodes[1]->point.X()};
    std::array<double, 2> y{nodes[0]->point.Y(), nodes[1]->point.Y()};

    const double rnorm = 0.5*nodes[0]->point.Distance(nodes[1]->point);

    math::Vector M(NODES_PER_ELEM);
    constexpr size_t GN = 1;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        const double s = xi->x;
        const double X = 0.5*(x[0]*(1-s) + x[1]*(1+s));
        const double Y = 0.5*(y[0]*(1-s) + y[1]*(1+s));
        M += xi->w*N_mat_1dof(X, Y);
    }
    M *= t*rnorm;

    return M;
};

math::Matrix Q4::get_B(const gp_Pnt& point) const{
    const double x = point.X();
    const double y = point.Y();

    math::Matrix B({
    b[0] + d[0]*y
    ,
    0
    ,
    b[1] + d[1]*y
    ,
    0
    ,
    b[2] + d[2]*y
    ,
    0
    ,
    b[3] + d[3]*y
    ,
    0
    ,
    0
    ,
    c[0] + d[0]*x
    ,
    0
    ,
    c[1] + d[1]*x
    ,
    0
    ,
    c[2] + d[2]*x
    ,
    0
    ,
    c[3] + d[3]*x
    ,
    c[0] + d[0]*x
    ,
    b[0] + d[0]*y
    ,
    c[1] + d[1]*x
    ,
    b[1] + d[1]*y
    ,
    c[2] + d[2]*x
    ,
    b[2] + d[2]*y
    ,
    c[3] + d[3]*x
    ,
    b[3] + d[3]*y
    }, S_SIZE, K_DIM);
    return B;
}

math::Matrix Q4::get_nodal_density_gradient(gp_Pnt p) const{
    (void)p;

    return this->dN_mat_1dof(p.X(), p.Y());
}

math::Matrix  Q4::get_Ni(const gp_Pnt& p) const{
    return this->N_mat(p.X(), p.Y());
}

}


