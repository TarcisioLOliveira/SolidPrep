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

    this->R = Eigen::Matrix<double, 3, 3>
        {{d2.X(), d1.X(), n.X()},
         {d2.Y(), d1.Y(), n.Y()},
         {d2.Z(), d1.Z(), n.Z()}};

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

    std::array<double, N*N> M = 
        {1, x[0], y[0],
         1, x[1], y[1],
         1, x[2], y[2]};

    std::array<int, N> ipiv;

    // M*C = I -> C=M^-1
    // C = {a[0], a[1], a[2],
    //      b[0], b[1], b[2],
    //      c[0], c[1], c[2]}
    int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating LU in CTRI3.", info);
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating computing inverse from LU in CTRI3.", info);

    for(size_t i = 0; i < N; ++i){
        this->a[i] = M[i];
        this->b[i] = M[i+3];
        this->c[i] = M[i+6];
    }

    this->delta = 0.5*(v1.Crossed(v2)).Magnitude();
}

std::vector<double> CTRI3::get_frictionless_Ge() const{
    const auto& gli = utils::GaussLegendreTri<2*ORDER>::get();
    std::vector<double> NN(2*E_KW, 0);
    std::vector<double> MnMn(NODES_PER_ELEM*2*E_KW, 0);

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
        std::fill(NN.begin(), NN.end(), 0);
        for(size_t i = 0; i < E_KW; ++i){
            for(size_t j = 0; j < DIM; ++j){
                NN[i] += N1[i + j*E_KW]*n.Coord(1+j);
                NN[i + E_KW] += N2[i + j*E_KW]*n.Coord(1+j);
            }
        }
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            for(size_t j = 0; j < 2*E_KW; ++j){
                MnMn[2*E_KW*i + j] += it->w*Nb[i]*NN[j];
            }
        }
    }
    cblas_dscal(MnMn.size(), this->delta, MnMn.data(), 1);

    return MnMn;
}
std::vector<double> CTRI3::fl_uLne(const std::vector<double>& D, const std::vector<double>& ln_e) const{
    Eigen::Matrix<double, S_SIZE, S_SIZE> Dm;
    for(size_t i = 0; i < S_SIZE; ++i){
        for(size_t j = 0; j < S_SIZE; ++j){
            Dm(i, j) = D[i*S_SIZE + j];
        }
    }
    auto e_info = this->e1->get_element_info();
    size_t U_KW = e_info->get_k_dimension();

    const auto& gli = utils::GaussLegendreTri<1>::get();

    gp_Pnt p = this->R_GS_point(gli.begin()->a, gli.begin()->b, gli.begin()->c);

    auto B = this->eps_mat_lambda(this->get_normal());
    auto Bu = this->e1->get_B(GS_point(gli.begin()->a, gli.begin()->b, gli.begin()->c));
    auto N = this->N_mat_1dof(p);
    Eigen::Vector<double, NODES_PER_ELEM> lv;
    double l_i = 0;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        lv[i] = ln_e[i];
        l_i += N[i]*ln_e[i];
    }
    Eigen::MatrixXd Bum;
    Bum.resize(S_SIZE, U_KW);
    for(size_t i = 0; i < S_SIZE; ++i){
        for(size_t j = 0; j < U_KW; ++j){
            Bum(i, j) = Bu[i*U_KW + j];
        }
    }
    Eigen::MatrixXd Mm = this->delta*(Bum.transpose()*Dm*(l_i*B + (B*lv)*N.transpose()));
    std::vector<double> M(U_KW*NODES_PER_ELEM);
    for(size_t i = 0; i < U_KW; ++i){
        for(size_t j = 0; j < NODES_PER_ELEM; ++j){
            M[i*NODES_PER_ELEM + j] = Mm(i, j);
        }
    }

    return M;
}
std::vector<double> CTRI3::fl_LnLne(const std::vector<double>& D, const std::vector<double>& ln_e) const{
    Eigen::Matrix<double, S_SIZE, S_SIZE> Dm;
    for(size_t i = 0; i < S_SIZE; ++i){
        for(size_t j = 0; j < S_SIZE; ++j){
            Dm(i, j) = D[i*S_SIZE + j];
        }
    }
    auto e_info = this->e1->get_element_info();
    size_t U_KW = e_info->get_k_dimension();

    const auto& gli = utils::GaussLegendreTri<2>::get();

    auto B = this->eps_mat_lambda(this->get_normal());
    Eigen::Vector<double, NODES_PER_ELEM> lv;
    Eigen::MatrixXd Mm;
    Mm.resize(NODES_PER_ELEM, NODES_PER_ELEM);
    Mm.fill(0);

    for(auto gl = gli.begin(); gl < gli.end(); ++gl){
        gp_Pnt p = this->R_GS_point(gli.begin()->a, gli.begin()->b, gli.begin()->c);
        auto N = this->N_mat_1dof(p);
        double l_i = 0;
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            lv[i] = ln_e[i];
            l_i += N[i]*ln_e[i];
        }
        Eigen::MatrixXd Bs = l_i*B + (B*lv)*N.transpose();
        Mm += gl->w*Bs.transpose()*Dm*Bs;
    }
    std::vector<double> M(U_KW*NODES_PER_ELEM);
    for(size_t i = 0; i < U_KW; ++i){
        for(size_t j = 0; j < NODES_PER_ELEM; ++j){
            M[i*NODES_PER_ELEM + j] = this->delta*Mm(i, j);
        }
    }

    return M;
}
std::vector<double> CTRI3::fl_LtLne(const std::vector<double>& D, const std::vector<double>& ln_e, size_t ti) const{
    Eigen::Matrix<double, S_SIZE, S_SIZE> Dm;
    for(size_t i = 0; i < S_SIZE; ++i){
        for(size_t j = 0; j < S_SIZE; ++j){
            Dm(i, j) = D[i*S_SIZE + j];
        }
    }
    auto e_info = this->e1->get_element_info();
    size_t U_KW = e_info->get_k_dimension();

    const auto& gli = utils::GaussLegendreTri<1>::get();

    gp_Pnt p = this->R_GS_point(gli.begin()->a, gli.begin()->b, gli.begin()->c);

    const gp_Dir nt = (ti == 0) ? this->get_d2() : this->get_d1();

    auto B = this->eps_mat_lambda(this->get_normal());
    auto Bt = this->eps_mat_lambda(nt);
    auto N = this->N_mat_1dof(p);
    Eigen::Vector<double, NODES_PER_ELEM> lv;
    double l_i = 0;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        lv[i] = ln_e[i];
        l_i += N[i]*ln_e[i];
    }
    Eigen::MatrixXd Mm = this->delta*(Bt.transpose()*Dm*(l_i*B + (B*lv)*N.transpose()));
    std::vector<double> M(U_KW*NODES_PER_ELEM);
    for(size_t i = 0; i < U_KW; ++i){
        for(size_t j = 0; j < NODES_PER_ELEM; ++j){
            M[i*NODES_PER_ELEM + j] = Mm(i, j);
        }
    }

    return M;
}
std::vector<double> CTRI3::fl_uLte(const std::vector<double>& D, size_t ti) const{
    Eigen::Matrix<double, S_SIZE, S_SIZE> Dm;
    for(size_t i = 0; i < S_SIZE; ++i){
        for(size_t j = 0; j < S_SIZE; ++j){
            Dm(i, j) = D[i*S_SIZE + j];
        }
    }
    auto e_info = this->e1->get_element_info();
    size_t U_KW = e_info->get_k_dimension();

    const auto& gli = utils::GaussLegendreTri<1>::get();

    const gp_Dir nt = (ti == 0) ? this->get_d2() : this->get_d1();

    auto Bt = this->eps_mat_lambda(nt);
    auto Bu = this->e1->get_B(this->GS_point(gli.begin()->a, gli.begin()->b, gli.begin()->c));
    Eigen::MatrixXd Bum;
    Bum.resize(S_SIZE, U_KW);
    for(size_t i = 0; i < S_SIZE; ++i){
        for(size_t j = 0; j < U_KW; ++j){
            Bum(i, j) = Bu[i*U_KW + j];
        }
    }
    Eigen::MatrixXd Mm = this->delta*(Bum.transpose()*Dm*Bt);
    std::vector<double> M(U_KW*NODES_PER_ELEM);
    for(size_t i = 0; i < U_KW; ++i){
        for(size_t j = 0; j < NODES_PER_ELEM; ++j){
            M[i*NODES_PER_ELEM + j] = Mm(i, j);
        }
    }

    return M;
}
std::vector<double> CTRI3::fl_LtLte(const std::vector<double>& D, size_t t1, size_t t2) const{
    Eigen::Matrix<double, S_SIZE, S_SIZE> Dm;
    for(size_t i = 0; i < S_SIZE; ++i){
        for(size_t j = 0; j < S_SIZE; ++j){
            Dm(i, j) = D[i*S_SIZE + j];
        }
    }

    const gp_Dir nt1 = (t1 == 0) ? this->get_d2() : this->get_d1();
    const gp_Dir nt2 = (t2 == 0) ? this->get_d2() : this->get_d1();

    auto B1 = this->eps_mat_lambda(nt1);
    auto B2 = this->eps_mat_lambda(nt2);

    Eigen::Matrix<double, NODES_PER_ELEM, NODES_PER_ELEM> Mm = this->delta*(B1.transpose()*Dm*B2);

    std::vector<double> M(NODES_PER_ELEM*NODES_PER_ELEM);
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        for(size_t j = 0; j < NODES_PER_ELEM; ++j){
            M[i*NODES_PER_ELEM + j] = Mm(i, j);
        }
    }

    return M;
}

std::vector<double> CTRI3::fl2_uL(const std::vector<double>& l_e) const{
    auto e_info = this->e1->get_element_info();
    size_t U_KW = e_info->get_k_dimension();
    const auto& gli = utils::GaussLegendreTri<3*ORDER>::get();
    std::vector<double> NN(2*U_KW, 0);
    std::vector<double> uL(2*U_KW*NODES_PER_ELEM, 0);

    const gp_Dir n = this->get_normal();

    for(auto it = gli.begin(); it != gli.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const gp_Pnt rpi = this->R_GS_point(it->a, it->b, it->c);
        const auto N1 = e1->get_Ni(pi);
        const auto N2 = e2->get_Ni(pi);
        const auto Nl = this->N_mat_1dof(rpi);
        double l = 0;
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            l += Nl[i]*l_e[i];
        }
        std::fill(NN.begin(), NN.end(), 0);
        for(size_t i = 0; i < U_KW; ++i){
            for(size_t j = 0; j < DIM; ++j){
                NN[i] += N1[i + j*U_KW]*n.Coord(1+j);
                NN[i + U_KW] -= N2[i + j*U_KW]*n.Coord(1+j);
            }
        }
        for(size_t i = 0; i < 2*U_KW; ++i){
            for(size_t j = 0; j < NODES_PER_ELEM; ++j){
                uL[i*NODES_PER_ELEM + j] += it->w*l*NN[i]*Nl[j];
            }
        }
    }
    cblas_dscal(uL.size(), this->delta, uL.data(), 1);

    return uL;
}
std::vector<double> CTRI3::fl2_LL(const std::vector<double>& l_e, const std::vector<double>& u1, const std::vector<double>& u2) const{
    auto e_info = this->e1->get_element_info();
    size_t U_KW = u1.size();
    const auto& gli = utils::GaussLegendreTri<4*ORDER>::get();
    std::vector<double> LL(NODES_PER_ELEM*NODES_PER_ELEM, 0);

    const gp_Dir n = this->get_normal();

    std::vector<double> up1(DIM), up2(DIM);
    for(auto it = gli.begin(); it != gli.end(); ++it){
        std::fill(up1.begin(), up1.end(), 0);
        std::fill(up2.begin(), up2.end(), 0);
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const gp_Pnt rpi = this->R_GS_point(it->a, it->b, it->c);
        const auto N1 = e1->get_Ni(pi);
        const auto N2 = e2->get_Ni(pi);
        const auto Nl = this->N_mat_1dof(rpi);
        double gp = 0;
        double l = 0;
        for(size_t j = 0; j < DIM; ++j){
            for(size_t i = 0; i < U_KW; ++i){
                up1[j] += N1[j*U_KW + i]*u1[i];
                up2[j] += N2[j*U_KW + i]*u2[i];
            }
            gp += (up2[j] - up1[j])*n.Coord(1+j);
        }
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            l += Nl[i]*l_e[i];
        }
        const double mult = 3*l*l/2 - gp;
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            for(size_t j = 0; j < NODES_PER_ELEM; ++j){
                LL[i*NODES_PER_ELEM + j] += it->w*mult*Nl[i]*Nl[j];
            }
        }
    }
    cblas_dscal(LL.size(), this->delta, LL.data(), 1);

    return LL;
}
void CTRI3::fl2_Ku_lambda(const double EPS, const std::vector<long> u1_pos, const std::vector<long> u2_pos, const std::vector<long>& lu_pos, const std::vector<double>& u, std::vector<double>& Ku) const{
    auto e_info = this->e1->get_element_info();
    const size_t bnum = this->NODES_PER_ELEM;
    size_t U_KW = u1_pos.size();
    const auto& gli1 = utils::GaussLegendreTri<3*ORDER>::get();
    std::vector<double> NN(2*U_KW, 0);
    std::vector<double> uL(2*U_KW, 0);
    std::vector<double> u1(U_KW), u2(U_KW);
    std::vector<double> l_e(bnum);
    for(size_t i = 0; i < bnum; ++i){
        l_e[i] = u[lu_pos[i]];
    }
    for(size_t i = 0; i < U_KW; ++i){
        u1[i] = u[u1_pos[i]];
        u2[i] = u[u2_pos[i]];
    }

    const gp_Dir n = this->get_normal();

    for(auto it = gli1.begin(); it != gli1.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const gp_Pnt rpi = this->R_GS_point(it->a, it->b, it->c);
        const auto N1 = e1->get_Ni(pi);
        const auto N2 = e2->get_Ni(pi);
        const auto Nl = this->N_mat_1dof(rpi);
        double l = 0;
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            l += Nl[i]*l_e[i];
        }
        std::fill(NN.begin(), NN.end(), 0);
        for(size_t i = 0; i < U_KW; ++i){
            for(size_t j = 0; j < DIM; ++j){
                NN[i] += N1[i + j*U_KW]*n.Coord(1+j);
                NN[i + U_KW] -= N2[i + j*U_KW]*n.Coord(1+j);
            }
        }
        for(size_t i = 0; i < 2*U_KW; ++i){
            uL[i] += it->w*l*l*NN[i]/2;
        }
    }
    for(size_t i = 0; i < U_KW; ++i){
        Ku[u1_pos[i]] += EPS*this->delta*uL[i];
        Ku[u2_pos[i]] += EPS*this->delta*uL[i + U_KW];
    }

    const auto& gli2 = utils::GaussLegendreTri<4*ORDER>::get();

    std::vector<double> LL(NODES_PER_ELEM, 0);

    std::vector<double> up1(DIM), up2(DIM);
    for(auto it = gli2.begin(); it != gli2.end(); ++it){
        std::fill(up1.begin(), up1.end(), 0);
        std::fill(up2.begin(), up2.end(), 0);
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const gp_Pnt rpi = this->R_GS_point(it->a, it->b, it->c);
        const auto N1 = e1->get_Ni(pi);
        const auto N2 = e2->get_Ni(pi);
        const auto Nl = this->N_mat_1dof(rpi);
        double gp = 0;
        double l = 0;
        for(size_t j = 0; j < DIM; ++j){
            for(size_t i = 0; i < U_KW; ++i){
                up1[j] += N1[j*U_KW + i]*u1[i];
                up2[j] += N2[j*U_KW + i]*u2[i];
            }
            gp += (up2[j] - up1[j])*n.Coord(1+j);
        }
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            l += Nl[i]*l_e[i];
        }
        const double mult = l*(l*l/2 - gp);
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            LL[i] += it->w*mult*Nl[i];
        }
    }

    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        Ku[lu_pos[i]] += EPS*this->delta*LL[i];
    }
}

}
