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

#include <lapacke.h>
#include "boundary_element/BTRI6.hpp"
#include "logger.hpp"
#include "utils/gauss_legendre.hpp"
#include "utils/boundary_nullifier.hpp"

namespace boundary_element{

BTRI6::BTRI6(ElementShape s, const MeshElement* const parent):
    BoundaryMeshElement(s.nodes, parent){
    const size_t N = BTRI6::NODES_PER_ELEM;
    
    std::array<double, N> x, y, z;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
        z[i] = this->nodes[i]->point.Z();
    }
    std::array<double, N*N> M = 
        {1, x[0], y[0], x[0]*y[0], x[0]*x[0], y[0]*y[0],
         1, x[1], y[1], x[1]*y[1], x[1]*x[1], y[1]*y[1],
         1, x[2], y[2], x[2]*y[2], x[2]*x[2], y[2]*y[2],
         1, x[3], y[3], x[3]*y[3], x[3]*x[3], y[3]*y[3],
         1, x[4], y[4], x[4]*y[4], x[4]*x[4], y[4]*y[4],
         1, x[5], y[5], x[5]*y[5], x[5]*x[5], y[5]*y[5]};

    std::array<int, N> ipiv;

    // M*C = I -> C=M^-1
    // C = {a[0], a[1], a[2],
    //      b[0], b[1], b[2],
    //      c[0], c[1], c[2]}
    int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating LU in BTRI6.", info);
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating computing inverse from LU in BTRI6.", info);

    for(size_t i = 0; i < N; ++i){
        this->a[i] = M[i];
        this->b[i] = M[i+N];
        this->c[i] = M[i+2*N];
        this->d[i] = M[i+3*N];
        this->e[i] = M[i+4*N];
        this->f[i] = M[i+5*N];
    }

    gp_Vec v1(this->nodes[0]->point, this->nodes[1]->point);
    gp_Vec v2(this->nodes[0]->point, this->nodes[2]->point);

    this->delta = 0.5*(v1.Crossed(v2)).Magnitude();
}

Eigen::MatrixXd BTRI6::diffusion_1dof(const Eigen::MatrixXd& A) const{
    const auto& gsi = utils::GaussLegendreTri<2>::get();
    Eigen::Matrix<double, 6, 6> result;
    result.fill(0);
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto dNN = this->dN_mat_1dof(p);
        result += it->w*dNN.transpose()*A*dNN;
    }
    return result*delta;
}
Eigen::MatrixXd BTRI6::advection_1dof(const Eigen::VectorXd& v) const{
    const auto& gsi = utils::GaussLegendreTri<3>::get();
    Eigen::Matrix<double, 6, 6> result;
    result.fill(0);
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto dNN = this->dN_mat_1dof(p);
        const auto  NN = this->N_mat_1dof(p);
        result += it->w*NN*(v.transpose()*dNN);
    }
    return result*delta;
}
Eigen::MatrixXd BTRI6::absorption_1dof() const{
    const auto& gsi = utils::GaussLegendreTri<4>::get();
    Eigen::Matrix<double, 6, 6> result;
    result.fill(0);
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto NN = this->N_mat_1dof(p);
        result += it->w*NN*NN.transpose();
    }
    return result*delta;
}
Eigen::VectorXd BTRI6::source_1dof() const{
    const auto& gsi = utils::GaussLegendreTri<2>::get();
    Eigen::Vector<double, 6> result{0,0,0,0,0,0};
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto NN = this->N_mat_1dof(p);
        result += it->w*NN;
    }
    return result*delta;
}
Eigen::VectorXd BTRI6::grad_1dof_upos(const gp_Pnt& p, const std::vector<double>& phi) const{
    (void)p;
    Eigen::Vector<double, 6> phiv{0, 0, 0, 0, 0, 0};
    for(size_t i = 0; i < 6; ++i){
        const auto p = this->nodes[i]->u_pos[0];
        if(p > -1){
            phiv[i] = phi[p];
        }
    }
    return this->dN_mat_1dof(p)*phiv;
};
Eigen::VectorXd BTRI6::grad_1dof_id(const gp_Pnt& p, const std::vector<double>& phi) const{
    (void)p;
    Eigen::Vector<double, 6> phiv{0, 0, 0, 0, 0, 0};
    for(size_t i = 0; i < 6; ++i){
        const auto p = this->nodes[i]->id;
        phiv[i] = phi[p];
    }
    return this->dN_mat_1dof(p)*phiv;
};
Eigen::VectorXd BTRI6::dF_2dof_id(const gp_Pnt& p, const std::vector<double>& phi) const{
    (void)p;
    Eigen::Vector<double, 12> phiv;
    phiv.fill(0);
    for(size_t i = 0; i < 6; ++i){
        const auto p = this->nodes[i]->id;
        phiv[2*i+0] = phi[2*p+0];
        phiv[2*i+1] = phi[2*p+1];
    }
    return this->dF_mat_2dof(p)*phiv;
}

Eigen::MatrixXd BTRI6::int_grad_phi() const{
    Eigen::Matrix<double, 2, 6> result;
    result.fill(0);

    const auto& gli = utils::GaussLegendreTri<2>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto dN = this->dN_mat_1dof(p);
        result += it->w*dN;
    }
    result *= delta;

    return result;
};
Eigen::MatrixXd BTRI6::int_grad_phi_x(const gp_Pnt& center) const{
    Eigen::Matrix<double, 2, 6> result;
    result.fill(0);

    const auto& gli = utils::GaussLegendreTri<2>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dx = p.X() - center.X();
        const auto dN = this->dN_mat_1dof(p);
        result += it->w*dN*dx;
    }
    result *= delta;

    return result;
};
Eigen::MatrixXd BTRI6::int_grad_phi_y(const gp_Pnt& center) const{
    Eigen::Matrix<double, 2, 6> result;
    result.fill(0);

    const auto& gli = utils::GaussLegendreTri<2>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dy = p.Y() - center.Y();
        const auto dN = this->dN_mat_1dof(p);
        result += it->w*dN*dy;
    }
    result *= delta;

    return result;
};
Eigen::MatrixXd BTRI6::int_grad_xi() const{
    Eigen::Matrix<double, 2, 6> result;
    result.fill(0);

    const auto& gli = utils::GaussLegendreTri<2>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto dxi = this->dxi_mat_1dof2(p);
        result += it->w*dxi;
    }
    result *= delta;

    return result;
};
Eigen::MatrixXd BTRI6::int_grad_xi_x(const gp_Pnt& center) const{
    Eigen::Matrix<double, 2, 6> result;
    result.fill(0);

    const auto& gli = utils::GaussLegendreTri<2>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dx = p.X() - center.X();
        const auto dxi = this->dxi_mat_1dof2(p);
        result += it->w*dxi*dx;
    }
    result *= delta;

    return result;
};
Eigen::MatrixXd BTRI6::int_grad_xi_y(const gp_Pnt& center) const{
    Eigen::Matrix<double, 2, 6> result;
    result.fill(0);

    const auto& gli = utils::GaussLegendreTri<2>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dy = p.Y() - center.Y();
        const auto dxi = this->dxi_mat_1dof2(p);
        result += it->w*dxi*dy;
    }
    result *= delta;

    return result;
};
Eigen::MatrixXd BTRI6::int_grad_F() const{
    Eigen::Matrix<double, 3, 12> result;
    result.fill(0);

    const auto& gli = utils::GaussLegendreTri<2>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto dN = this->dF_mat_2dof(p);
        result += it->w*dN;
    }
    result *= delta;

    return result;
};
Eigen::MatrixXd BTRI6::int_grad_F_x(const gp_Pnt& center) const{
    Eigen::Matrix<double, 3, 12> result;
    result.fill(0);

    const auto& gli = utils::GaussLegendreTri<2>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dx = p.X() - center.X();
        const auto dN = this->dF_mat_2dof(p);
        result += it->w*dN*dx;
    }
    result *= delta;

    return result;
};
Eigen::MatrixXd BTRI6::int_grad_F_y(const gp_Pnt& center) const{
    Eigen::Matrix<double, 3, 12> result;
    result.fill(0);

    const auto& gli = utils::GaussLegendreTri<2>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dy = p.Y() - center.Y();
        const auto dN = this->dF_mat_2dof(p);
        result += it->w*dN*dy;
    }
    result *= delta;

    return result;
};
Eigen::VectorXd BTRI6::int_N_AzBz(const gp_Pnt& center, const double Az, const double Bz) const{
    Eigen::Vector<double, 6> result;
    result.fill(0);

    const auto& gli = utils::GaussLegendreTri<3>::get();
    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dx = p.X() - center.X();
        const double dy = p.Y() - center.Y();
        const auto N = this->N_mat_1dof(p);
        result += it->w*N*(Az*dx+Bz*dy);
    }
    result *= delta;

    return result;
}
Eigen::VectorXd BTRI6::int_N_x(const gp_Pnt& center) const{
    const auto& gsi = utils::GaussLegendreTri<3>::get();
    Eigen::Vector<double, 6>  result{0, 0, 0, 0, 0, 0};
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dx = p.X() - center.X();
        const auto NN = this->N_mat_1dof(p);
        result += it->w*NN*dx;
    }

    return delta*result;
}
Eigen::VectorXd BTRI6::int_N_y(const gp_Pnt& center) const{
    const auto& gsi = utils::GaussLegendreTri<3>::get();
    Eigen::Vector<double, 6>  result{0, 0, 0, 0, 0, 0};
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dy = p.Y() - center.Y();
        const auto NN = this->N_mat_1dof(p);
        result += it->w*NN*dy;
    }

    return delta*result;
}
Eigen::MatrixXd BTRI6::int_grad_F_t2_t1(const Eigen::MatrixXd& B3, const gp_Pnt& center) const{
    Eigen::Matrix<double, 12, 4> result;
    result.fill(0);
    const auto& gsi = utils::GaussLegendreTri<3>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dx = p.X() - center.X();
        const double dy = p.Y() - center.Y();
        const double dx2 = dx*dx;
        const double dy2 = dy*dy;
        const auto dF = this->dF_mat_2dof(p);
        Eigen::Matrix<double, 2, 4> T
            {{dy2, dx2,   0,   0},
             {  0,   0, dx2, dy2}};
        result += it->w*dF.transpose()*B3*T;
    }
    return delta*result;
}
Eigen::MatrixXd BTRI6::int_grad_phi_t2_t1(const Eigen::MatrixXd& B2, const gp_Pnt& center) const{
    Eigen::Matrix<double, 6, 4> result;
    result.fill(0);
    const auto& gsi = utils::GaussLegendreTri<3>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dx = p.X() - center.X();
        const double dy = p.Y() - center.Y();
        const double dx2 = dx*dx;
        const double dy2 = dy*dy;
        const auto dphi = this->dN_mat_1dof(p);
        Eigen::Matrix<double, 2, 4> T
            {{dy2, dx2,   0,   0},
             {  0,   0, dx2, dy2}};
        result += it->w*dphi.transpose()*B2*T;
    }
    return delta*result;
}
Eigen::MatrixXd BTRI6::int_grad_F_D(const Eigen::MatrixXd& a, const gp_Pnt& center) const{
    Eigen::Matrix<double, 12, 3> result;
    result.fill(0);
    const auto& gsi = utils::GaussLegendreTri<2>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dx = p.X() - center.X();
        const double dy = p.Y() - center.Y();
        const auto dF = this->dF_mat_2dof(p);
        Eigen::Matrix<double, 3, 3> ABC
            {{dx, dy, 1},
             {dx, dy, 1},
             {dx, dy, 1}};
        result += it->w*dF.transpose()*a*ABC;
    }
    return delta*result;
}
Eigen::MatrixXd BTRI6::int_grad_phi_D(const Eigen::MatrixXd& a, const gp_Pnt& center) const{
    Eigen::Matrix<double, 6, 3> result;
    result.fill(0);
    const auto& gsi = utils::GaussLegendreTri<2>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dx = p.X() - center.X();
        const double dy = p.Y() - center.Y();
        const auto dphi = this->dN_mat_1dof(p);
        Eigen::Matrix<double, 2, 3> ABC
            {{dx, dy, 1},
             {dx, dy, 1}};
        result += it->w*dphi.transpose()*a*ABC;
    }
    return delta*result;
}


Eigen::MatrixXd BTRI6::L4(const Eigen::MatrixXd& B) const{
    const auto& gli = utils::GaussLegendreTri<2>::get();
    Eigen::Matrix<double, 12, 12> result;
    result.fill(0);

    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto dF = this->dF_mat_2dof(p);

        result += it->w*dF.transpose()*B*dF;
    }
    result *= delta;

    return result;
}
Eigen::MatrixXd BTRI6::L3(const Eigen::MatrixXd& B) const{
    const auto& gli = utils::GaussLegendreTri<2>::get();
    Eigen::Matrix<double, 12, 6> result;
    result.fill(0);

    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto dF = this->dF_mat_2dof(p);
        const auto dN = this->dN_mat_1dof(p);

        result += it->w*dF.transpose()*B*dN;
    }
    result *= delta;

    return result;
}
Eigen::MatrixXd BTRI6::L2(const Eigen::MatrixXd& B) const{
    const auto& gli = utils::GaussLegendreTri<2>::get();
    Eigen::Matrix<double, 6, 6> result;
    result.fill(0);

    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto dN = this->dN_mat_1dof(p);

        result += it->w*dN.transpose()*B*dN;
    }
    result *= delta;

    return result;
}

Eigen::MatrixXd BTRI6::L3xi(const Eigen::MatrixXd& B) const{
    const auto& gli = utils::GaussLegendreTri<2>::get();
    Eigen::Matrix<double, 12, 6> result;
    result.fill(0);

    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto dF = this->dF_mat_2dof(p);
        const auto dxi = this->dxi_mat_1dof(p);

        result += it->w*dF.transpose()*B*dxi;
    }
    result *= delta;

    return result;
}
Eigen::MatrixXd BTRI6::L2xi(const Eigen::MatrixXd& B) const{
    const auto& gli = utils::GaussLegendreTri<2>::get();
    Eigen::Matrix<double, 6, 6> result;
    result.fill(0);

    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto dN = this->dN_mat_1dof(p);
        const auto dxi = this->dxi_mat_1dof(p);

        result += it->w*dN.transpose()*B*dxi;
    }
    result *= delta;

    return result;
}
Eigen::MatrixXd BTRI6::L4chi(const Eigen::MatrixXd& B) const{
    const auto& gli = utils::GaussLegendreTri<2>::get();
    Eigen::Matrix<double, 12, 6> result;
    result.fill(0);

    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto dF = this->dF_mat_2dof(p);
        const auto dchi = this->dchi_mat_1dof(p);

        result += it->w*dF.transpose()*B*dchi;
    }
    result *= delta;

    return result;
}
Eigen::MatrixXd BTRI6::L3Tchi(const Eigen::MatrixXd& B) const{
    const auto& gli = utils::GaussLegendreTri<2>::get();
    Eigen::Matrix<double, 6, 6> result;
    result.fill(0);

    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto dN = this->dN_mat_1dof(p);
        const auto dchi = this->dchi_mat_1dof(p);

        result += it->w*dN.transpose()*B.transpose()*dchi;
    }
    result *= delta;

    return result;
}
Eigen::MatrixXd BTRI6::L4zeta(const Eigen::MatrixXd& B) const{
    const auto& gli = utils::GaussLegendreTri<2>::get();
    Eigen::Matrix<double, 12, 6> result;
    result.fill(0);

    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto dF = this->dF_mat_2dof(p);
        const auto dzeta = this->dzeta_mat_1dof(p);

        result += it->w*dF.transpose()*B*dzeta;
    }
    result *= delta;

    return result;
}
Eigen::MatrixXd BTRI6::L3Tzeta(const Eigen::MatrixXd& B) const{
    const auto& gli = utils::GaussLegendreTri<2>::get();
    Eigen::Matrix<double, 6, 6> result;
    result.fill(0);

    for(auto it = gli.begin(); it < gli.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const auto dN = this->dN_mat_1dof(p);
        const auto dzeta = this->dzeta_mat_1dof(p);

        result += it->w*dN.transpose()*B.transpose()*dzeta;
    }
    result *= delta;

    return result;
}
}
