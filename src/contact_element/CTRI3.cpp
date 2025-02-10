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
#include "utils/gauss_legendre.hpp"
#include "element_factory.hpp"

namespace contact_element{
CTRI3::CTRI3(ElementShape s, const MeshElement* const e1, const MeshElement* const e2):
    ContactMeshElement(s.nodes, e1, e2){

    this->E_KW = e1->get_element_info()->get_k_dimension();

    const size_t N = CTRI3::NODES_PER_ELEM;
    
    gp_Vec v1(this->nodes[0]->point, this->nodes[1]->point);
    gp_Vec v2(this->nodes[0]->point, this->nodes[2]->point);

    gp_Vec nv = v1.Crossed(v2);
    const gp_Pnt c1 = this->get_centroid();
    const gp_Pnt c2 = this->e1->get_centroid();
    const gp_Dir n2(gp_Vec(c1, c2));

    //logger::quick_log(c1.X(), c1.Y(), c1.Z());
    //logger::quick_log(c2.X(), c2.Y(), c2.Z());
    //logger::quick_log(n2.X(), n2.Y(), n2.Z());
    //logger::quick_log(nv.X(), nv.Y(), nv.Z());
    //logger::quick_log(n2.Dot(nv));
    //exit(0);

    if(n2.Dot(nv) > 0){
        nv *= -1;
    }


    const gp_Dir n(nv);
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
math::Matrix CTRI3::fl_uLne(const math::Matrix & D, const math::Vector& ln_e) const{
    auto e_info = this->e1->get_element_info();

    const auto& gli = utils::GaussLegendreTri<1>::get();

    gp_Pnt p = this->R_GS_point(gli.begin()->a, gli.begin()->b, gli.begin()->c);

    auto B = this->eps_mat_lambda(this->get_normal());
    auto Bu = this->e1->get_B(GS_point(gli.begin()->a, gli.begin()->b, gli.begin()->c));
    auto N = this->N_mat_1dof(p);
    double l_i = N.T()*ln_e;

    math::Matrix M(this->delta*(Bu.T()*D*(l_i*B + (B*ln_e)*N.T())));

    return M;
}
math::Matrix CTRI3::fl_LnLne(const math::Matrix& D, const math::Vector& ln_e) const{
    auto e_info = this->e1->get_element_info();

    const auto& gli = utils::GaussLegendreTri<2>::get();

    auto B = this->eps_mat_lambda(this->get_normal());

    math::Matrix M(NODES_PER_ELEM, NODES_PER_ELEM);
    for(auto gl = gli.begin(); gl < gli.end(); ++gl){
        gp_Pnt p = this->R_GS_point(gli.begin()->a, gli.begin()->b, gli.begin()->c);
        auto N = this->N_mat_1dof(p);
        const double l_i = N.T()*ln_e;
        math::Matrix Bs = l_i*B + (B*ln_e)*N.T();
        M += gl->w*Bs.T()*D*Bs;
    }
    M *= this->delta;

    return M;
}
math::Matrix CTRI3::fl_LtLne(const math::Matrix& D, const math::Vector& ln_e, size_t ti) const{
    auto e_info = this->e1->get_element_info();

    const auto& gli = utils::GaussLegendreTri<1>::get();

    gp_Pnt p = this->R_GS_point(gli.begin()->a, gli.begin()->b, gli.begin()->c);

    const gp_Dir nt = (ti == 0) ? this->get_d2() : this->get_d1();

    auto B = this->eps_mat_lambda(this->get_normal());
    auto Bt = this->eps_mat_lambda(nt);
    auto N = this->N_mat_1dof(p);
    const double l_i = N.T()*ln_e;
    math::Matrix M = this->delta*(Bt.T()*D*(l_i*B + (B*ln_e)*N.T()));

    return M;
}
math::Matrix CTRI3::fl_uLte(const math::Matrix& D, size_t ti) const{
    auto e_info = this->e1->get_element_info();

    const auto& gli = utils::GaussLegendreTri<1>::get();

    const gp_Dir nt = (ti == 0) ? this->get_d2() : this->get_d1();

    auto Bt = this->eps_mat_lambda(nt);
    auto Bu = this->e1->get_B(this->GS_point(gli.begin()->a, gli.begin()->b, gli.begin()->c));
    math::Matrix M = this->delta*(Bu.T()*D*Bt);

    return M;
}
math::Matrix CTRI3::fl_LtLte(const math::Matrix& D, size_t t1, size_t t2) const{
    const gp_Dir nt1 = (t1 == 0) ? this->get_d2() : this->get_d1();
    const gp_Dir nt2 = (t2 == 0) ? this->get_d2() : this->get_d1();

    auto B1 = this->eps_mat_lambda(nt1);
    auto B2 = this->eps_mat_lambda(nt2);

    math::Matrix M = this->delta*(B1.T()*D*B2);

    return M;
}

math::Matrix CTRI3::fl2_uL(const math::Vector& l_e) const{
    auto e_info = this->e1->get_element_info();
    size_t U_KW = e_info->get_k_dimension();
    const auto& gli = utils::GaussLegendreTri<3*ORDER>::get();
    math::Vector NN(2*U_KW);
    math::Matrix uL(2*U_KW, NODES_PER_ELEM);

    const gp_Dir n = this->get_normal();

    for(auto it = gli.begin(); it != gli.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const gp_Pnt rpi = this->R_GS_point(it->a, it->b, it->c);
        const auto N1 = e1->get_Ni(pi);
        const auto N2 = e2->get_Ni(pi);
        const auto Nl = this->N_mat_1dof(rpi);
        const double l = Nl.T()*l_e;
        NN.fill(0);
        for(size_t i = 0; i < U_KW; ++i){
            for(size_t j = 0; j < DIM; ++j){
                NN[i] -= N1(j,i)*n.Coord(1+j);
                NN[i + U_KW] += N2(j,i)*n.Coord(1+j);
            }
        }
        uL += (it->w*(-l))*(NN*Nl.T());
    }
    uL *= this->delta;

    return uL;
}
math::Matrix CTRI3::fl2_LL(const math::Vector& l_e, const math::Vector& u1, const math::Vector& u2) const{
    auto e_info = this->e1->get_element_info();
    const auto& gli = utils::GaussLegendreTri<4*ORDER>::get();
    math::Matrix LL(NODES_PER_ELEM, NODES_PER_ELEM);

    const gp_Dir n = this->get_normal();

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
        const double mult = 3*l*l/2 - gp;
        LL += (it->w*mult)*(Nl*Nl.T());
    }
    LL *= this->delta;

    return LL;
}
void CTRI3::fl2_Ku_lambda(const double EPS, const std::vector<long> u1_pos, const std::vector<long> u2_pos, const std::vector<long>& lu_pos, const std::vector<double>& u, std::vector<double>& Ku) const{
    auto e_info = this->e1->get_element_info();
    const size_t bnum = this->NODES_PER_ELEM;
    size_t U_KW = u1_pos.size();
    const auto& gli1 = utils::GaussLegendreTri<4*ORDER>::get();
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

    const gp_Dir n = this->get_normal();
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
        const double mult = (gp - l*l/2);
        NN.fill(0);
        for(size_t i = 0; i < U_KW; ++i){
            for(size_t j = 0; j < DIM; ++j){
                NN[i] -= N1(j,i)*n.Coord(1+j);
                NN[i + U_KW] += N2(j,i)*n.Coord(1+j);
            }
        }
        uL += (it->w*mult)*NN;
        LL += (it->w*(-l*mult))*Nl;
    }
    for(size_t i = 0; i < U_KW; ++i){
        Ku[u1_pos[i]] += EPS*this->delta*uL[i];
        Ku[u2_pos[i]] += EPS*this->delta*uL[i + U_KW];
    }
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        Ku[lu_pos[i]] += EPS*this->delta*LL[i];
    }
}

}
