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
#include "logger.hpp"
#include "math/matrix.hpp"
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
    const gp_Pnt c1 = this->get_centroid();
    const gp_Pnt c2 = this->e1->get_centroid();
    const gp_Dir n2(gp_Vec(c1, c2));

    // SIMPLE
    if(e1_base){
        if(n2.Dot(nv) > 0){
            nv *= -1;
        }
    } else {
    // CONSTR
        if(n2.Dot(nv) < 0){
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
math::Matrix CTRI3::fl3_uL(const math::Matrix& D, const math::Vector& ln_e) const{
    auto e_info = this->e1->get_element_info();
    size_t U_KW = e_info->get_k_dimension();

    const gp_Dir n1(this->get_normal());
    const gp_Dir n2(this->get_d1());
    const gp_Dir n3(this->get_d2());

    const auto& gli = utils::GaussLegendreTri<2>::get();

    math::Matrix M(U_KW, 3*NODES_PER_ELEM);

    for(auto gl = gli.begin(); gl < gli.end(); ++gl){
        gp_Pnt p = this->R_GS_point(gl->a, gl->b, gl->c);
        gp_Pnt pn = this->GS_point(gl->a, gl->b, gl->c);

        auto B = this->delta_eps_1(p, ln_e, n1, n2, n3);
        auto Bu = this->e1->get_B(pn);

        M += gl->w*(Bu.T()*D*B);
    }

    M *= this->delta;

    return M;
}
math::Matrix CTRI3::fl3_LL(const math::Matrix& D, const math::Vector& ln_e, const math::Vector& lp1_e, const math::Vector& lp2_e, const std::vector<double>& u) const{
    auto e_info = this->e1->get_element_info();

    const gp_Dir n1(this->get_normal());
    const gp_Dir n2(this->get_d1());
    const gp_Dir n3(this->get_d2());

    const auto& gli = utils::GaussLegendreTri<3>::get();

    math::Matrix M(3*NODES_PER_ELEM, 3*NODES_PER_ELEM);
    math::Matrix Mn(3*NODES_PER_ELEM, 3*NODES_PER_ELEM);

    for(auto gl = gli.begin(); gl < gli.end(); ++gl){
        gp_Pnt p = this->R_GS_point(gl->a, gl->b, gl->c);
        gp_Pnt pn = this->GS_point(gl->a, gl->b, gl->c);

        auto eps_l = this->eps_vec(p, ln_e, lp1_e, lp2_e, n1, n2, n3);
        auto strain = this->e1->get_strain_vector(pn, u);

        auto B = this->delta_eps_1(p, ln_e, n1, n2, n3);

        //auto N = this->N_mat_1dof(p);
        //auto Bl = this->eps_mat_lambda(n1);
        
        //math::Matrix Mn_tmp(N*(eps_l + strain).T()*D*Bl);
        //for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        //    for(size_t j = 0; j < NODES_PER_ELEM; ++j){
        //        Mn(i,j) = 2*(Mn_tmp(i,j) + Mn_tmp(j,i));
        //    }
        //}
        Mn = this->delta_eps_2(p, n1, eps_l + strain, D);

        M += gl->w*((B.T()*D*B) + Mn);
    }

    M *= this->delta;

    return M;
}
void CTRI3::fl3_Ku(const math::Matrix& D, const std::vector<long> u_pos, const std::vector<long>& lu_pos, const std::vector<double>& u, std::vector<double>& Ku) const{
    auto e_info = this->e1->get_element_info();
    size_t U_KW = e_info->get_k_dimension();

    math::Vector u_e(U_KW);
    math::Vector ln_e(NODES_PER_ELEM);
    math::Vector lp1_e(NODES_PER_ELEM);
    math::Vector lp2_e(NODES_PER_ELEM);
    for(size_t i = 0; i < U_KW; ++i){
        u_e[i] = u[u_pos[i]];
    }
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        size_t j = i + NODES_PER_ELEM;
        size_t k = j + NODES_PER_ELEM;
        ln_e[i] = u[lu_pos[i]];
        lp1_e[i] = u[lu_pos[j]];
        lp2_e[i] = u[lu_pos[k]];
    }

    const gp_Dir n1(this->get_normal());
    const gp_Dir n2(this->get_d1());
    const gp_Dir n3(this->get_d2());

    const auto& gli = utils::GaussLegendreTri<3>::get();

    math::Vector uL(U_KW);
    math::Vector LL(3*NODES_PER_ELEM);

    for(auto gl = gli.begin(); gl < gli.end(); ++gl){
        gp_Pnt p = this->R_GS_point(gl->a, gl->b, gl->c);
        gp_Pnt pn = this->GS_point(gl->a, gl->b, gl->c);

        const auto eps_l = this->eps_vec(p, ln_e, lp1_e, lp2_e, n1, n2, n3);
        const auto strain = this->e1->get_strain_vector(pn, u);
        const auto full_strain = eps_l + strain;

        auto B = this->delta_eps_1(p, ln_e, n1, n2, n3);
        auto Bu = this->e1->get_B(pn);
        uL += gl->w*(Bu.T()*D*eps_l);
        LL += gl->w*(B.T()*D*full_strain);
    }
    for(size_t i = 0; i < U_KW; ++i){
        Ku[u_pos[i]] += this->delta*uL[i];
    }
    for(size_t i = 0; i < 3*NODES_PER_ELEM; ++i){
        Ku[lu_pos[i]] += this->delta*LL[i];
    }
}
void CTRI3::fl3_dKu(const math::Matrix& D, const double eta, const std::vector<long> u_pos, const std::vector<long>& lu_pos, const std::vector<double>& u, const std::vector<double>& du, std::vector<double>& Ku) const{
    auto e_info = this->e1->get_element_info();
    size_t U_KW = e_info->get_k_dimension();

    math::Vector u_e(U_KW);
    math::Vector du_e(U_KW);
    math::Vector ln_e(NODES_PER_ELEM);
    math::Vector lp1_e(NODES_PER_ELEM);
    math::Vector lp2_e(NODES_PER_ELEM);
    math::Vector dln_e(NODES_PER_ELEM);
    math::Vector dlp1_e(NODES_PER_ELEM);
    math::Vector dlp2_e(NODES_PER_ELEM);
    for(size_t i = 0; i < U_KW; ++i){
        u_e[i]  = u[u_pos[i]];
        du_e[i] = du[u_pos[i]];
    }
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        size_t j = i + NODES_PER_ELEM;
        size_t k = j + NODES_PER_ELEM;
        ln_e[i]  = u[lu_pos[i]];
        lp1_e[i] = u[lu_pos[j]];
        lp2_e[i] = u[lu_pos[k]];
        dln_e[i]  = du[lu_pos[i]];
        dlp1_e[i] = du[lu_pos[j]];
        dlp2_e[i] = du[lu_pos[k]];
    }

    const gp_Dir n1(this->get_normal());
    const gp_Dir n2(this->get_d1());
    const gp_Dir n3(this->get_d2());

    const auto& gli = utils::GaussLegendreTri<3>::get();

    math::Vector uL(U_KW);
    math::Vector LL(3*NODES_PER_ELEM);

    for(auto gl = gli.begin(); gl < gli.end(); ++gl){
        gp_Pnt p = this->R_GS_point(gl->a, gl->b, gl->c);
        gp_Pnt pn = this->GS_point(gl->a, gl->b, gl->c);

        const auto eps_l = this->eps_vec(p, ln_e, lp1_e, lp2_e, n1, n2, n3);
        const auto deps_l = this->eps_vec_eta(p, eta, ln_e, dln_e, dlp1_e, dlp2_e, n1, n2, n3);

        const auto strain = this->e1->get_strain_vector(pn, u);
        const auto dstrain = this->e1->get_strain_vector(pn, du);
        const auto full_strain = eps_l + strain;
        const auto dfull_strain = deps_l + dstrain;

        const auto eps_eta = this->eps_D_eta(p, dln_e, n1, full_strain, D);

        auto B = this->delta_eps_1(p, ln_e, n1, n2, n3);
        auto Bu = this->e1->get_B(pn);

        uL += gl->w*(Bu.T()*D*deps_l);
        LL += gl->w*(B.T()*D*dfull_strain + eps_eta);
    }
    for(size_t i = 0; i < U_KW; ++i){
        Ku[u_pos[i]] += this->delta*uL[i];
    }
    for(size_t i = 0; i < 3*NODES_PER_ELEM; ++i){
        Ku[lu_pos[i]] += this->delta*LL[i];
    }
}

math::Vector CTRI3::fl3_eq(const math::Vector& ln_e, const math::Vector& lp1_e, const math::Vector& lp2_e, const math::Vector& u_e) const{
    auto e_info = this->e1->get_element_info();
    size_t U_KW = e_info->get_k_dimension();
    math::Vector V(U_KW);
    const auto& gli = utils::GaussLegendreTri<4>::get();

    for(auto gl = gli.begin(); gl < gli.end(); ++gl){
        gp_Pnt p = this->R_GS_point(gl->a, gl->b, gl->c);
        gp_Pnt pn = this->GS_point(gl->a, gl->b, gl->c);

        const math::Vector N(this->N_mat_1dof(p));
        const double l_n = N.T()*ln_e;
        const double l_p1 = N.T()*lp1_e;
        const double l_p2 = N.T()*lp2_e;

        const math::Vector lv({l_p2, l_p1, l_n*l_n});
        
        const math::Matrix Nu1(this->e1->get_Ni(pn));
        const math::Matrix Nu2(this->e2->get_Ni(pn));
        const math::Vector u1(Nu1*u_e);

        //V += gl->w*(Nu2.T()*(u1 + this->R*lv));
        V += gl->w*(Nu2.T()*(this->R*lv));
    }
    V *= this->delta;

    return V;
}

math::Vector CTRI3::fl3_eq(const math::Vector& ln_e, const math::Vector& lp1_e, const math::Vector& lp2_e, const math::Vector& u_e, const size_t dof) const{
    (void)u_e;
    auto e_info = this->e1->get_element_info();
    //size_t U_KW = e_info->get_k_dimension();
    size_t U_NE = e_info->get_nodes_per_element();
    math::Vector V(U_NE);
    const auto& gli = utils::GaussLegendreTri<4>::get();

    for(auto gl = gli.begin(); gl < gli.end(); ++gl){
        gp_Pnt p = this->R_GS_point(gl->a, gl->b, gl->c);
        gp_Pnt pn = this->GS_point(gl->a, gl->b, gl->c);

        const math::Vector N(this->N_mat_1dof(p));
        const double l_n = N.T()*ln_e;
        const double l_p1 = N.T()*lp1_e;
        const double l_p2 = N.T()*lp2_e;

        const math::Vector lv({l_p2, l_p1, l_n*l_n});
        
        const math::Vector Nu2(this->e2->get_Ni_1dof(pn));
        //const math::Vector Nu1(this->e1->get_Ni_1dof(pn));
        //const math::Vector u1(Nu1*u_e);
        double l_i = 0;
        for(size_t i = 0; i < 3; ++i){
            l_i += this->R(dof, i)*lv[i];
        }

        //V += gl->w*(Nu2*(u1[dof] + l_i));
        V += gl->w*(Nu2*l_i);
    }
    V *= this->delta;

    return V;
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

math::Vector CTRI3::fl2_LAG(const math::Vector& l_e, const math::Vector& u1, const math::Vector& u2) const{
    auto e_info = this->e1->get_element_info();
    size_t U_KW = u1.get_N();
    const auto& gli1 = utils::GaussLegendreTri<4*ORDER>::get();
    math::Vector NN(2*U_KW, 0);
    math::Vector uL(2*U_KW, 0);

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
        const double mult_uL = (gp-l*l/2);
        const double mult_LL = -l*(gp - l*l/2);
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
    math::Vector result(2*U_KW + NODES_PER_ELEM);
    for(size_t i = 0; i < 2*U_KW; ++i){
        result[i] = this->delta*uL[i];
    }
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        result[i + 2*U_KW] = this->delta*LL[i];
    }

    return result;
}
double CTRI3::fl2_int(const math::Vector& l_e, const math::Vector& u1, const math::Vector& u2) const{
    auto e_info = this->e1->get_element_info();
    size_t U_KW = u1.get_N();
    const auto& gli1 = utils::GaussLegendreTri<4*ORDER>::get();
    math::Vector NN(2*U_KW, 0);
    math::Vector uL(2*U_KW, 0);

    const gp_Dir n = this->get_normal();
    math::Vector up1(DIM), up2(DIM);
    math::Vector LL(NODES_PER_ELEM);
    double g = 0;
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
        const double g_tmp = gp-l*l/2;
        g += it->w*g_tmp*g_tmp/2;
    }

    return this->delta*g;
}
double CTRI3::fl2_int_deriv(const math::Vector& l_e, const math::Vector& u1, const math::Vector& u2, const math::Vector& dl_e, const math::Vector& du1, const math::Vector& du2, const double eta) const{
    (void) eta;
    auto e_info = this->e1->get_element_info();
    size_t U_KW = u1.get_N();
    const auto& gli1 = utils::GaussLegendreTri<4*ORDER>::get();
    math::Vector NN(2*U_KW, 0);
    math::Vector uL(2*U_KW, 0);

    const gp_Dir n = this->get_normal();
    math::Vector up1(DIM), up2(DIM);
    math::Vector dup1(DIM), dup2(DIM);
    math::Vector LL(NODES_PER_ELEM);
    double dg = 0;
    for(auto it = gli1.begin(); it != gli1.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const gp_Pnt rpi = this->R_GS_point(it->a, it->b, it->c);
        const auto N1 = e1->get_Ni(pi);
        const auto N2 = e2->get_Ni(pi);
        const auto Nl = this->N_mat_1dof(rpi);
        double gp = 0, dgp = 0;
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
        const double g_tmp = gp - l*l/2;
        const double dg_tmp = dgp - l*dl;

        dg += it->w*g_tmp*dg_tmp;
    }

    return this->delta*dg;
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
        const double mult_uL = (gp-l*l/2);
        const double mult_LL = -l*(gp - l*l/2);
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


void CTRI3::fl2_dKu_lambda(const double EPS, const std::vector<long> u1_pos, const std::vector<long> u2_pos, const std::vector<long>& lu_pos, const std::vector<double>& u, const std::vector<double>& du, const double eta, std::vector<double>& dKu) const{
    (void) eta;
    auto e_info = this->e1->get_element_info();
    const size_t bnum = this->NODES_PER_ELEM;
    size_t U_KW = u1_pos.size();
    const auto& gli1 = utils::GaussLegendreTri<4*ORDER>::get();
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

    const gp_Dir n = this->get_normal();
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
        const double mult_uL = (dgp - l*dl);
        const double mult_LL = (3*l*l*dl/2 - (dl*gp + l*dgp));
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
        dKu[u1_pos[i]] += EPS*this->delta*uL[i];
        dKu[u2_pos[i]] += EPS*this->delta*uL[i + U_KW];
    }
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        dKu[lu_pos[i]] += EPS*this->delta*LL[i];
    }
}

}
