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
#include "element_factory.hpp"
#include "logger.hpp"
#include "utils/gauss_legendre.hpp"

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

std::vector<double> BTRI3::get_K_ext(const Eigen::MatrixXd& D, const gp_Pnt& center) const{
    Eigen::Matrix<double, K_DIM, K_DIM + 6> Km;
    Km.fill(0);
    const auto& gsi = utils::GaussLegendreTri<ORDER>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const auto eps = this->get_eps(pi, center);
        const auto deps = this->get_deps();
        Km += it->w*(deps.transpose()*D*eps);
    }
    Km *= this->delta;
    std::vector<double> K(K_DIM*(K_DIM+6));
    for(size_t i = 0; i < K_DIM; ++i){
        for(size_t j = 0; j < K_DIM + 6; ++j){
            K[i*(K_DIM + 6) + j] = Km(i,j);
        }
    }

    return K;
}

std::vector<double> BTRI3::get_normal_stresses(const Eigen::MatrixXd& D, const std::vector<double>& u, const gp_Pnt& p, const gp_Pnt& center) const{
    const auto eps = this->get_eps(p, center);
    std::vector<double> E(S_SIZE, 0);
    std::vector<double> ES(S_SIZE, 0);
    std::vector<double> S(3, 0);
    std::vector<double> vals(K_DIM + 6, 0);
    const size_t offset = u.size() - 6;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        size_t pos = this->nodes[i]->id;//u_pos[0];
        for(size_t j = 0; j < NODE_DOF; ++j){
            vals[i*NODE_DOF + j] = u[NODE_DOF*pos + j];
        }
    }
    for(size_t i = 0; i < 6; ++i){
        vals[K_DIM + i] = u[offset + i];
    }

    for(size_t i = 0; i < S_SIZE; ++i){
        for(size_t j = 0; j < K_DIM + 6; ++j){
            E[i] += eps(i, j)*vals[j];
        }
    }
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < S_SIZE; ++j){
            S[i] += D((i + 2), j)*E[j] + ES[j];
        }
    }

    return S;
}

std::vector<double> BTRI3::get_stress_integrals(const Eigen::MatrixXd& D, const gp_Pnt& center) const{
    // M_y M_x V_z M_z V_y V_x
    std::vector<double> M(6*(K_DIM + 6), 0);
    Eigen::Matrix<double, S_SIZE, K_DIM + 6> E;
    const auto& gsi = utils::GaussLegendreTri<ORDER + 1>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const double dx = pi.X() - center.X();
        const double dy = pi.Y() - center.Y();
        const auto eps = this->get_eps(pi, center);
        E = D*eps;
        for(size_t i = 0; i < K_DIM + 6; ++i){
            M[0*(K_DIM + 6) + i] += it->w*E(2, i)*dx;
            M[1*(K_DIM + 6) + i] += it->w*E(2, i)*dy;
            M[2*(K_DIM + 6) + i] += it->w*E(2, i);
            M[3*(K_DIM + 6) + i] += it->w*(E(3, i)*dx - E(4, i)*dy);
            M[4*(K_DIM + 6) + i] += it->w*E(3, i);
            M[5*(K_DIM + 6) + i] += it->w*E(4, i);
        }
    }
    cblas_dscal(M.size(), this->delta, M.data(), 1);

    return M;
}

std::vector<double> BTRI3::get_equilibrium_partial(const Eigen::MatrixXd& D, const gp_Pnt& center, const std::vector<size_t>& stresses) const{
    const size_t KW = this->K_DIM; // workaround that's necessary for some reason
    Eigen::MatrixXd Km(KW, K_DIM + 6);
    Km.fill(0);
    const auto& gsi = utils::GaussLegendreTri<ORDER>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const auto eps = this->get_eps(pi, center);
        const auto deps = this->get_deps();
        const Eigen::MatrixXd Ktmp = D*eps;
        for(size_t i = 0; i < K_DIM; ++i){
            for(size_t j = 0; j < K_DIM + 6; ++j){
                for(size_t k = 0; k < stresses.size(); ++k){
                    Km(i,j) += it->w*(deps(stresses[k], i)*Ktmp(stresses[k], j));
                }
            }
        }
    }
    Km *= this->delta;
    std::vector<double> K(K_DIM*(K_DIM+6));
    for(size_t i = 0; i < K_DIM; ++i){
        for(size_t j = 0; j < K_DIM + 6; ++j){
            K[i*(K_DIM + 6) + j] = Km(i,j);
        }
    }

    return K;
}

std::vector<double> BTRI3::get_dz_vector(const Eigen::MatrixXd& S, const Eigen::MatrixXd& D, const double Az, const double Bz, const gp_Pnt& center) const{
    const size_t KW = this->K_DIM; // workaround that's necessary for some reason
    std::vector<double> vec(KW, 0);
    const auto& gsi = utils::GaussLegendreTri<ORDER+1>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const double dx = pi.X() - center.X();
        const double dy = pi.Y() - center.Y();
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            const double N = this->N(pi, i);
            vec[i*NODE_DOF + 0] += -it->w*N*D(4,2)*(Az*dx + Bz*dy)/(S(2,2)*D(2,2));
            vec[i*NODE_DOF + 1] += -it->w*N*D(3,2)*(Az*dx + Bz*dy)/(S(2,2)*D(2,2));
            vec[i*NODE_DOF + 2] += -it->w*N*(Az*dx + Bz*dy)/S(2,2);
        }
    }
    cblas_dscal(vec.size(), this->delta, vec.data(), 1);

    return vec;
}


std::vector<double> BTRI3::get_force_vector(const Eigen::MatrixXd& D, const std::vector<double>& u, const gp_Pnt& center, const Eigen::MatrixXd& rot) const{
    const size_t KW_P = this->parent->get_element_info()->get_k_dimension();

    std::vector<double> F(KW_P, 0);

    const auto& gsi = utils::GaussLegendreTri<ORDER+2>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const Eigen::Vector<double, 3> pv(rot*Eigen::Vector<double, 3>(pi.X(), pi.Y(), pi.Z()));
        const gp_Pnt pr(pv[0], pv[1], pv[2]);
        auto S = this->get_normal_stresses(D, u, pi, center);
        const Eigen::Vector<double, 3> Sv{S[2], S[1], S[0]};
        const Eigen::Vector<double, 3> Sr(rot*Sv);
        const auto N = this->parent->get_Ni(pr);
        for(size_t i = 0; i < 3; ++i){
            for(size_t j = 0; j < KW_P; ++j){
                F[j] += it->w*N[i*KW_P + j]*Sr[i];
            }
        }
    }
    cblas_dscal(KW_P, this->delta, F.data(), 1);

    return F;
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

}
