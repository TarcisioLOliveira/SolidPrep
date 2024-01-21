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
#include "boundary_element/BTRI3.hpp"
#include "logger.hpp"
#include "utils/gauss_legendre.hpp"
#include "utils/boundary_nullifier.hpp"

namespace boundary_element{

BTRI3::BTRI3(ElementShape s, const MeshElement* const parent):
    BoundaryMeshElement(s.nodes, parent){
    const size_t N = BTRI3::NODES_PER_ELEM;
    
    std::array<double, N> x, y, z;
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X();
        y[i] = this->nodes[i]->point.Y();
        z[i] = this->nodes[i]->point.Z();
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
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating LU in BTRI3.", info);
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR, N, M.data(), N, ipiv.data());
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating computing inverse from LU in BTRI3.", info);

    for(size_t i = 0; i < N; ++i){
        this->a[i] = M[i];
        this->b[i] = M[i+3];
        this->c[i] = M[i+6];
    }

    gp_Vec v1(this->nodes[0]->point, this->nodes[1]->point);
    gp_Vec v2(this->nodes[0]->point, this->nodes[2]->point);

    this->delta = 0.5*(v1.Crossed(v2)).Magnitude();
}

Eigen::MatrixXd BTRI3::diffusion_1dof(const Eigen::MatrixXd& A) const{
    const auto B = this->dN_mat_1dof();
    return delta*B.transpose()*A*B;
}
Eigen::MatrixXd BTRI3::advection_1dof(const Eigen::VectorXd& v) const{
    const auto NN = this->N_mat_1dof(this->GS_point(1.0/3.0, 1.0/3.0, 1.0/3.0));
    const auto B = this->dN_mat_1dof();
    return this->delta*NN*(v.transpose()*B);
}
Eigen::MatrixXd BTRI3::absorption_1dof() const{
    const auto& gsi = utils::GaussLegendreTri<2>::get();
    Eigen::Matrix<double, 3, 3> result{{0, 0, 0},{0, 0, 0},{0, 0, 0}};
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const auto NN = this->N_mat_1dof(this->GS_point(it->a, it->b, it->c));
        result += it->w*NN*NN.transpose();
    }

    return delta*result;
}
Eigen::VectorXd BTRI3::source_1dof() const{
    return this->delta*this->N_mat_1dof(this->GS_point(1.0/3.0, 1.0/3.0, 1.0/3.0));
}
Eigen::VectorXd BTRI3::flow_1dof(const std::array<const Node*, 2>& nodes) const{
    gp_Pnt p = nodes[0]->point;
    p.BaryCenter(1, nodes[1]->point, 1);
    const double d = nodes[0]->point.Distance(nodes[1]->point);

    return d*this->N_mat_1dof(p);
}
Eigen::VectorXd BTRI3::grad_1dof_upos(const gp_Pnt& p, const std::vector<double>& phi) const{
    (void)p;
    Eigen::Vector<double, 3> phiv{0, 0, 0};
    for(size_t i = 0; i < 3; ++i){
        const auto p = this->nodes[i]->u_pos[0];
        if(p > -1){
            phiv[i] = phi[p];
        }
    }
    return this->dN_mat_1dof()*phiv;
};
Eigen::VectorXd BTRI3::grad_1dof_id(const gp_Pnt& p, const std::vector<double>& phi) const{
    (void)p;
    Eigen::Vector<double, 3> phiv{0, 0, 0};
    for(size_t i = 0; i < 3; ++i){
        const auto p = this->nodes[i]->id;
        phiv[i] = phi[p];
    }
    return this->dN_mat_1dof()*phiv;
};

Eigen::VectorXd BTRI3::dF_2dof_id(const gp_Pnt& p, const std::vector<double>& phi) const{
    (void)p;
    Eigen::Vector<double, 6> phiv{0, 0, 0, 0, 0, 0};
    for(size_t i = 0; i < 3; ++i){
        const auto p = this->nodes[i]->id;
        phiv[2*i+0] = phi[2*p+0];
        phiv[2*i+1] = phi[2*p+1];
    }
    return this->dF_mat_2dof()*phiv;
}

Eigen::MatrixXd BTRI3::int_grad_phi() const{
    return delta*this->dN_mat_1dof();
};
Eigen::MatrixXd BTRI3::int_grad_phi_x(const gp_Pnt& center) const{
    const auto p = this->GS_point(1.0/3.0, 1.0/3.0, 1.0/3.0);
    const double dx = p.X() - center.X();
    auto dN = this->dN_mat_1dof();
    for(size_t i = 0; i < 3; ++i){
        dN(0, i) *= delta*dx;
        dN(1, i) *= delta*dx;
    }

    return dN;
};
Eigen::MatrixXd BTRI3::int_grad_phi_y(const gp_Pnt& center) const{
    const auto p = this->GS_point(1.0/3.0, 1.0/3.0, 1.0/3.0);
    const double dy = p.Y() - center.Y();
    auto dN = this->dN_mat_1dof();
    for(size_t i = 0; i < 3; ++i){
        dN(0, i) *= delta*dy;
        dN(1, i) *= delta*dy;
    }

    return dN;
};
Eigen::MatrixXd BTRI3::int_grad_xi() const{
    return delta*this->dxi_mat_1dof2();
};
Eigen::MatrixXd BTRI3::int_grad_xi_x(const gp_Pnt& center) const{
    const auto p = this->GS_point(1.0/3.0, 1.0/3.0, 1.0/3.0);
    const double dx = p.X() - center.X();
    auto dxi = this->dxi_mat_1dof2();
    for(size_t i = 0; i < 3; ++i){
        dxi(0, i) *= delta*dx;
        dxi(1, i) *= delta*dx;
    }

    return dxi;
};
Eigen::MatrixXd BTRI3::int_grad_xi_y(const gp_Pnt& center) const{
    const auto p = this->GS_point(1.0/3.0, 1.0/3.0, 1.0/3.0);
    const double dy = p.Y() - center.Y();
    auto dxi = this->dxi_mat_1dof2();
    for(size_t i = 0; i < 3; ++i){
        dxi(0, i) *= delta*dy;
        dxi(1, i) *= delta*dy;
    }

    return dxi;
};
Eigen::MatrixXd BTRI3::int_grad_F() const{
    return delta*this->dF_mat_2dof();
};
Eigen::MatrixXd BTRI3::int_grad_F_x(const gp_Pnt& center) const{
    const auto p = this->GS_point(1.0/3.0, 1.0/3.0, 1.0/3.0);
    const double dx = p.X() - center.X();
    auto dF = this->dF_mat_2dof();
    for(size_t i = 0; i < 3; ++i){
        dF(0, i) *= delta*dx;
        dF(1, i) *= delta*dx;
        dF(2, i) *= delta*dx;
    }

    return dF;
};
Eigen::MatrixXd BTRI3::int_grad_F_y(const gp_Pnt& center) const{
    const auto p = this->GS_point(1.0/3.0, 1.0/3.0, 1.0/3.0);
    const double dy = p.Y() - center.Y();
    auto dF = this->dF_mat_2dof();
    for(size_t i = 0; i < 3; ++i){
        dF(0, i) *= delta*dy;
        dF(1, i) *= delta*dy;
        dF(2, i) *= delta*dy;
    }

    return dF;
};
Eigen::VectorXd BTRI3::int_N_x(const gp_Pnt& center) const{
    const auto& gsi = utils::GaussLegendreTri<2>::get();
    Eigen::Vector<double, 3>  result{0, 0, 0};
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dx = p.X() - center.X();
        const auto NN = this->N_mat_1dof(p);
        result += it->w*NN*dx;
    }

    return delta*result;
}
Eigen::VectorXd BTRI3::int_N_y(const gp_Pnt& center) const{
    const auto& gsi = utils::GaussLegendreTri<2>::get();
    Eigen::Vector<double, 3>  result{0, 0, 0};
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dy = p.Y() - center.Y();
        const auto NN = this->N_mat_1dof(p);
        result += it->w*NN*dy;
    }

    return delta*result;
}
Eigen::MatrixXd BTRI3::int_NdN(const std::vector<double>& phi) const{
    const auto p = this->GS_point(1.0/3.0, 1.0/3.0, 1.0/3.0);
    const auto N = this->N_mat_1dof(p);
    const auto dN = this->dN_mat_1dof();

    Eigen::Vector<double, 3> phiv{0, 0, 0};
    for(size_t i = 0; i < 3; ++i){
        const auto p = this->nodes[i]->id;
        phiv[i] = phi[p];
    }

    return this->delta*N*(dN*phiv).transpose();
}
Eigen::MatrixXd BTRI3::int_NdF(const std::vector<double>& phi) const{
    const auto p = this->GS_point(1.0/3.0, 1.0/3.0, 1.0/3.0);
    const auto N = this->N_mat_1dof(p);
    const auto dF = this->dF_mat_2dof();

    Eigen::Vector<double, 6> phiv{0, 0, 0, 0, 0, 0};
    for(size_t i = 0; i < 3; ++i){
        const auto p = this->nodes[i]->id;
        phiv[2*i+0] = phi[2*p+0];
        phiv[2*i+1] = phi[2*p+1];
    }

    return this->delta*N*(dF*phiv).transpose();
}

Eigen::MatrixXd BTRI3::int_grad_F_t2_t1(const Eigen::MatrixXd& B3, const gp_Pnt& center) const{
    Eigen::Matrix<double, 6, 4> result
        {{0, 0, 0, 0},
         {0, 0, 0, 0},
         {0, 0, 0, 0},
         {0, 0, 0, 0},
         {0, 0, 0, 0},
         {0, 0, 0, 0}};
    const auto& gsi = utils::GaussLegendreTri<3>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dx = p.X() - center.X();
        const double dy = p.Y() - center.Y();
        const double dx2 = dx*dx;
        const double dy2 = dy*dy;
        const auto dF = this->dF_mat_2dof();
        Eigen::Matrix<double, 2, 4> T
            {{dy2, dx2,   0,   0},
             {  0,   0, dx2, dy2}};
        result += it->w*dF.transpose()*B3*T;
    }
    return delta*result;
}
Eigen::MatrixXd BTRI3::int_grad_phi_t2_t1(const Eigen::MatrixXd& B2, const gp_Pnt& center) const{
    Eigen::Matrix<double, 3, 4> result
        {{0, 0, 0, 0},
         {0, 0, 0, 0},
         {0, 0, 0, 0}};
    const auto& gsi = utils::GaussLegendreTri<3>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dx = p.X() - center.X();
        const double dy = p.Y() - center.Y();
        const double dx2 = dx*dx;
        const double dy2 = dy*dy;
        const auto dphi = this->dN_mat_1dof();
        Eigen::Matrix<double, 2, 4> T
            {{dy2, dx2,   0,   0},
             {  0,   0, dx2, dy2}};
        result += it->w*dphi.transpose()*B2*T;
    }
    return delta*result;
}
Eigen::MatrixXd BTRI3::int_grad_F_D(const Eigen::MatrixXd& a, const gp_Pnt& center) const{
    Eigen::Matrix<double, 6, 3> result
        {{0, 0, 0},
         {0, 0, 0},
         {0, 0, 0},
         {0, 0, 0},
         {0, 0, 0},
         {0, 0, 0}};
    const auto& gsi = utils::GaussLegendreTri<2>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dx = p.X() - center.X();
        const double dy = p.Y() - center.Y();
        const auto dF = this->dF_mat_2dof();
        Eigen::Matrix<double, 3, 3> ABC
            {{dx, dy, 1},
             {dx, dy, 1},
             {dx, dy, 1}};
        result += it->w*dF.transpose()*a*ABC;
    }
    return delta*result;
}
Eigen::MatrixXd BTRI3::int_grad_phi_D(const Eigen::MatrixXd& a, const gp_Pnt& center) const{
    Eigen::Matrix<double, 3, 3> result
        {{0, 0, 0},
         {0, 0, 0},
         {0, 0, 0}};
    const auto& gsi = utils::GaussLegendreTri<2>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt p = this->GS_point(it->a, it->b, it->c);
        const double dx = p.X() - center.X();
        const double dy = p.Y() - center.Y();
        const auto dphi = this->dN_mat_1dof();
        Eigen::Matrix<double, 2, 3> ABC
            {{dx, dy, 1},
             {dx, dy, 1}};
        result += it->w*dphi.transpose()*a*ABC;
    }
    return delta*result;
}

Eigen::MatrixXd BTRI3::L4(const Eigen::MatrixXd& B) const{
    const auto dF = this->dF_mat_2dof();

    return delta*dF.transpose()*B*dF;
}
Eigen::MatrixXd BTRI3::L3(const Eigen::MatrixXd& B) const{
    const auto dF = this->dF_mat_2dof();
    const auto dN = this->dN_mat_1dof();

    return delta*dF.transpose()*B*dN;
}
Eigen::MatrixXd BTRI3::L2(const Eigen::MatrixXd& B) const{
    const auto dN = this->dN_mat_1dof();

    return delta*dN.transpose()*B*dN;
}

Eigen::MatrixXd BTRI3::L3xi(const Eigen::MatrixXd& B) const{
    const auto dF = this->dF_mat_2dof();
    const auto dxi = this->dxi_mat_1dof();

    return delta*dF.transpose()*B*dxi;
}
Eigen::MatrixXd BTRI3::L2xi(const Eigen::MatrixXd& B) const{
    const auto dN = this->dN_mat_1dof();
    const auto dxi = this->dxi_mat_1dof();

    return delta*dN.transpose()*B*dxi;
}

Eigen::MatrixXd BTRI3::L4chi(const Eigen::MatrixXd& B) const{
    const auto dF = this->dF_mat_2dof();
    const auto dchi = this->dchi_mat_1dof();

    return delta*dF.transpose()*B*dchi;
}

Eigen::MatrixXd BTRI3::L3Tchi(const Eigen::MatrixXd& B) const{
    const auto dphi = this->dN_mat_1dof();
    const auto dchi = this->dchi_mat_1dof();

    return delta*dphi.transpose()*B.transpose()*dchi;
}

Eigen::MatrixXd BTRI3::L4zeta(const Eigen::MatrixXd& B) const{
    const auto dF = this->dF_mat_2dof();
    const auto dzeta = this->dzeta_mat_1dof();

    return delta*dF.transpose()*B*dzeta;
}

Eigen::MatrixXd BTRI3::L3Tzeta(const Eigen::MatrixXd& B) const{
    const auto dphi = this->dN_mat_1dof();
    const auto dzeta = this->dzeta_mat_1dof();

    return delta*dphi.transpose()*B.transpose()*dzeta;

}

}
