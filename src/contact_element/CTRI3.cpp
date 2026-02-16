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

#include <lapacke.h>
#include "contact_element/CTRI3.hpp"
#include "math/matrix.hpp"
#include "math/slicer.hpp"
#include "utils/gauss_legendre.hpp"
#include "element_factory.hpp"

namespace contact_element{
CTRI3::CTRI3(ElementShape s, const MeshElement* const e1, const MeshElement* const e2, bool e1_base):
    ContactMeshElement(s.nodes, e1, e2){

    this->E_KW = e1->get_element_info()->get_k_dimension();

    const size_t N = CTRI3::NODES_PER_ELEM;
    
    gp_Vec v1(this->nodes[0]->point, this->nodes[1]->point);
    gp_Vec v2(this->nodes[0]->point, this->nodes[2]->point);

    gp_Vec nv = v1.Crossed(v2);
    const gp_Pnt c1 = this->e1->get_centroid();
    const gp_Pnt c2 = this->e2->get_centroid();
    const gp_Dir n2(gp_Vec(c1, c2));

    // SIMPLE
    if(e1_base){
        if(n2.Dot(nv) < 0){
            nv *= -1;
        }
    } else {
    // CONSTR
        if(n2.Dot(nv) > 0){
            nv *= -1;
        }
    }


    const gp_Dir n(nv);
    const gp_Dir d1(v1);
    const gp_Dir d2(d1.Crossed(n));

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

    this->delta = 0.5*(v1.Crossed(v2)).Magnitude();
}

math::Matrix CTRI3::get_frictionless_Ge() const{
    const auto& gli = utils::GaussLegendreTri<2*ORDER>::get();
    math::Vector NN(2*E_KW);
    math::Matrix MnMn(NODES_PER_ELEM, 2*E_KW);

    const gp_Dir n = this->get_normal();

    std::vector<double> up1(DIM), up2(DIM);
    for(auto it = gli.begin(); it != gli.end(); ++it){
        std::fill(up1.begin(), up1.end(), 0);
        std::fill(up2.begin(), up2.end(), 0);
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const gp_Pnt rpi = this->R_GS_point(it->a, it->b, it->c);
        const auto N1 = this->e1->get_Ni(pi);
        const auto N2 = this->e2->get_Ni(pi);
        const auto Nb = this->N_mat_1dof(rpi);
        NN.fill(0);
        for(size_t i = 0; i < E_KW; ++i){
            for(size_t j = 0; j < DIM; ++j){
                NN[i] += N1(j, i)*n.Coord(1+j);
                NN[i + E_KW] += N2(j, i)*n.Coord(1+j);
            }
        }
        MnMn += it->w*(Nb*NN.T());
    }
    MnMn *= this->delta;

    return MnMn;
}

math::Matrix CTRI3::fl2_uu(const math::Vector& l_e, const math::Vector& u1, const math::Vector& u2) const{
    auto e_info = this->e1->get_element_info();
    size_t U_KW = e_info->get_k_dimension();
    const auto& gli = utils::GaussLegendreTri<4*ORDER + 1>::get();
    math::Vector NN(2*U_KW);
    math::Matrix uu(2*U_KW, 2*U_KW);

    const gp_Dir n = -this->get_normal();

    math::Vector up1(DIM), up2(DIM);

    for(auto it = gli.begin(); it != gli.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const gp_Pnt rpi = this->R_GS_point(it->a, it->b, it->c);
        const auto N1 = e1->get_Ni(pi);
        const auto N2 = e2->get_Ni(pi);
        const auto Nl = this->N_mat_1dof(rpi);
        double gp = 0;
        const double l = Nl.T()*l_e;
        up1 = N1*u1;
        up2 = N2*u2;
        for(size_t j = 0; j < DIM; ++j){
            gp += (up2[j] - up1[j])*n.Coord(1+j);
        }
        NN.fill(0);
        for(size_t i = 0; i < U_KW; ++i){
            for(size_t j = 0; j < DIM; ++j){
                NN[i] -= N1(j,i)*n.Coord(1+j);
                NN[i + U_KW] += N2(j,i)*n.Coord(1+j);
            }
        }
        const double mult_uu = l*ddh(gp);
        uu += (it->w*mult_uu)*(NN*NN.T());
    }
    uu *= this->delta;

    return uu;
}
math::Matrix CTRI3::fl2_uL(const math::Vector& u1, const math::Vector& u2) const{
    auto e_info = this->e1->get_element_info();
    size_t U_KW = e_info->get_k_dimension();
    const auto& gli = utils::GaussLegendreTri<4*ORDER + 1>::get();
    math::Vector NN(2*U_KW);
    math::Matrix uL(2*U_KW, NODES_PER_ELEM);

    const gp_Dir n = -this->get_normal();

    math::Vector up1(DIM), up2(DIM);

    for(auto it = gli.begin(); it != gli.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const gp_Pnt rpi = this->R_GS_point(it->a, it->b, it->c);
        const auto N1 = e1->get_Ni(pi);
        const auto N2 = e2->get_Ni(pi);
        const auto Nl = this->N_mat_1dof(rpi);
        double gp = 0;
        up1 = N1*u1;
        up2 = N2*u2;
        for(size_t j = 0; j < DIM; ++j){
            gp += (up2[j] - up1[j])*n.Coord(1+j);
        }
        NN.fill(0);
        for(size_t i = 0; i < U_KW; ++i){
            for(size_t j = 0; j < DIM; ++j){
                NN[i] -= N1(j,i)*n.Coord(1+j);
                NN[i + U_KW] += N2(j,i)*n.Coord(1+j);
            }
        }
        const double mult_uL = dh(gp);
        uL += (it->w*mult_uL)*(NN*Nl.T());
    }
    uL *= this->delta;

    return uL;
}

void CTRI3::fl2_Ku_lambda(const double EPS, const std::vector<long> u1_pos, const std::vector<long> u2_pos, const std::vector<long>& lu_pos, const std::vector<double>& u, std::vector<double>& Ku) const{
    auto e_info = this->e1->get_element_info();
    const size_t bnum = this->NODES_PER_ELEM;
    size_t U_KW = u1_pos.size();
    const auto& gli1 = utils::GaussLegendreTri<4*ORDER+1>::get();
    math::Vector NN(2*U_KW, 0);
    math::Vector uL(2*U_KW, 0);
    math::Vector u1(U_KW), u2(U_KW);
    math::Vector l_e(bnum);
    for(size_t i = 0; i < bnum; ++i){
        l_e[i] = u[lu_pos[i]];
    }
    for(size_t i = 0; i < U_KW; ++i){
        u1[i] = u[u1_pos[i]];
        u2[i] = u[u2_pos[i]];
    }

    const gp_Dir n = -this->get_normal();
    math::Vector up1(DIM), up2(DIM);
    math::Vector LL(NODES_PER_ELEM);
    for(auto it = gli1.begin(); it != gli1.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const gp_Pnt rpi = this->R_GS_point(it->a, it->b, it->c);
        const auto N1 = e1->get_Ni(pi);
        const auto N2 = e2->get_Ni(pi);
        const auto Nl = this->N_mat_1dof(rpi);
        double gp = 0;
        const double l = Nl.T()*l_e;
        up1 = N1*u1;
        up2 = N2*u2;
        for(size_t j = 0; j < DIM; ++j){
            gp += (up2[j] - up1[j])*n.Coord(1+j);
        }
        const double mult_uL = l*dh(gp);
        const double mult_LL = h(gp);
        NN.fill(0);
        for(size_t i = 0; i < U_KW; ++i){
            for(size_t j = 0; j < DIM; ++j){
                NN[i] -= N1(j,i)*n.Coord(1+j);
                NN[i + U_KW] += N2(j,i)*n.Coord(1+j);
            }
        }
        uL += (it->w*mult_uL)*NN;
        LL += (it->w*mult_LL)*Nl;
    }
    for(size_t i = 0; i < U_KW; ++i){
        Ku[u1_pos[i]] += EPS*this->delta*uL[i];
        Ku[u2_pos[i]] += EPS*this->delta*uL[i + U_KW];
    }
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        Ku[lu_pos[i]] += EPS*this->delta*LL[i];
    }
}

void CTRI3::fl2_dKu_lambda(const double EPS, const std::vector<long> u1_pos, const std::vector<long> u2_pos, const std::vector<long>& lu_pos, const std::vector<double>& u, const std::vector<double>& du, std::vector<double>& Ku) const{
    auto e_info = this->e1->get_element_info();
    const size_t bnum = this->NODES_PER_ELEM;
    size_t U_KW = u1_pos.size();
    const auto& gli1 = utils::GaussLegendreTri<4*ORDER+1>::get();
    math::Vector NN(2*U_KW, 0);
    math::Vector uL(2*U_KW, 0);
    math::Vector u1(U_KW), u2(U_KW);
    math::Vector l_e(bnum);
    math::Vector du1(U_KW), du2(U_KW);
    math::Vector dl_e(bnum);
    for(size_t i = 0; i < bnum; ++i){
        l_e[i] = u[lu_pos[i]];
        dl_e[i] = du[lu_pos[i]];
    }
    for(size_t i = 0; i < U_KW; ++i){
        u1[i] = u[u1_pos[i]];
        u2[i] = u[u2_pos[i]];
        du1[i] = du[u1_pos[i]];
        du2[i] = du[u2_pos[i]];
    }

    const gp_Dir n = -this->get_normal();
    math::Vector up1(DIM), up2(DIM);
    math::Vector dup1(DIM), dup2(DIM);
    math::Vector LL(NODES_PER_ELEM);
    for(auto it = gli1.begin(); it != gli1.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const gp_Pnt rpi = this->R_GS_point(it->a, it->b, it->c);
        const auto N1 = e1->get_Ni(pi);
        const auto N2 = e2->get_Ni(pi);
        const auto Nl = this->N_mat_1dof(rpi);
        double gp = 0;
        double dgp = 0;
        const double l = Nl.T()*l_e;
        const double dl = Nl.T()*dl_e;
        up1 = N1*u1;
        up2 = N2*u2;
        dup1 = N1*du1;
        dup2 = N2*du2;
        for(size_t j = 0; j < DIM; ++j){
            gp += (up2[j] - up1[j])*n.Coord(1+j);
            dgp += (dup2[j] - dup1[j])*n.Coord(1+j);
        }
        const double mult_uL = dl*dh(gp) + l*ddh(gp)*dgp;
        const double mult_LL = dh(gp)*dgp;
        NN.fill(0);
        for(size_t i = 0; i < U_KW; ++i){
            for(size_t j = 0; j < DIM; ++j){
                NN[i] -= N1(j,i)*n.Coord(1+j);
                NN[i + U_KW] += N2(j,i)*n.Coord(1+j);
            }
        }
        uL += (it->w*mult_uL)*NN;
        LL += (it->w*mult_LL)*Nl;
    }
    for(size_t i = 0; i < U_KW; ++i){
        Ku[u1_pos[i]] += EPS*this->delta*uL[i];
        Ku[u2_pos[i]] += EPS*this->delta*uL[i + U_KW];
    }
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        Ku[lu_pos[i]] += EPS*this->delta*LL[i];
    }
}

math::Matrix CTRI3::diffusion_Ndof(const math::Matrix& A) const{
    const auto dN = this->dN_mat_dim_dof();

    return this->delta*dN.T()*A*dN;
}
math::Matrix CTRI3::absorption_Ndof() const{
    const auto& gli = utils::GaussLegendreTri<2*ORDER>::get();
    math::Matrix NN(K_DIM, K_DIM);

    for(auto it = gli.begin(); it != gli.end(); ++it){
        const gp_Pnt rpi = this->R_GS_point(it->a, it->b, it->c);
        const auto Nl = this->N_mat_3dof(rpi);
        NN += it->w*(Nl.T()*Nl);
    }
    NN *= this->delta;

    return NN;
}

math::Matrix CTRI3::diffusion_1dof(const math::Matrix& A) const{
    const auto dN = this->dN_mat_1dof();

    return this->delta*dN.T()*A*dN;
}
math::Matrix CTRI3::absorption_1dof() const{
    const auto& gli = utils::GaussLegendreTri<2*ORDER>::get();
    math::Matrix NN(NODES_PER_ELEM, NODES_PER_ELEM);

    for(auto it = gli.begin(); it != gli.end(); ++it){
        const gp_Pnt rpi = this->R_GS_point(it->a, it->b, it->c);
        const auto Nl = this->N_mat_1dof(rpi);
        NN += it->w*(Nl*Nl.T());
    }
    NN *= this->delta;

    return NN;
}

math::Vector CTRI3::lambda_source_log(const std::vector<double>& u_ext, const gp_Dir n, const double K) const{

    auto e_info = this->e1->get_element_info();
    const auto P_N = e_info->get_nodes_per_element();

    math::VectorSliceGeneralView uv1(u_ext, math::slicer::from_node_upos(this->e1->nodes, P_N, NODE_DOF));
    math::VectorSliceGeneralView uv2(u_ext, math::slicer::from_node_upos(this->e2->nodes, P_N, NODE_DOF));
    
    math::Vector NN(NODES_PER_ELEM);

    //const auto& gli = utils::GaussLegendreTri<6>::get();
    //for(auto it = gli.begin(); it != gli.end(); ++it){
    //    const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
    //    const gp_Pnt rpi = this->R_GS_point(it->a, it->b, it->c);
    //    const auto Nl = this->N_mat_1dof(rpi);
    //    const auto N1 = this->e1->get_Ni(pi);
    //    const auto N2 = this->e2->get_Ni(pi);
    //    double gp = 0;
    //    math::Vector up1(N1*uv1);
    //    math::Vector up2(N2*uv2);
    //    for(size_t j = 0; j < DIM; ++j){
    //        gp += (up2[j] - up1[j])*n.Coord(1+j);
    //    }
    //    const double h1 = this->H(gp, K);
    //    NN += (it->w*h1)*Nl;
    //}
    //NN *= this->delta;

    std::array<double, NODES_PER_ELEM> w;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        std::fill(w.begin(), w.end(), 0.0);
        w[i] = 1;
        const gp_Pnt pi = this->GS_point(w[0], w[1], w[2]);
        const gp_Pnt rpi = this->R_GS_point(w[0], w[1], w[2]);
        const auto Nl = this->N_mat_1dof(rpi);
        const auto N1 = this->e1->get_Ni(pi);
        const auto N2 = this->e2->get_Ni(pi);
        double gp = 0;
        math::Vector up1(N1*uv1);
        math::Vector up2(N2*uv2);
        for(size_t j = 0; j < DIM; ++j){
            gp += (up2[j] - up1[j])*n.Coord(1+j);
        }
        const double h1 = this->H(gp, K);
        NN[i] = h1;
    }

    return NN;
}

}
