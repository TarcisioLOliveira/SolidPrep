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
#include "math/slicer.hpp"
#include "math/matrix.hpp"
#include "math/vector_view.hpp"
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
    math::Matrix M(
        {1, x[0], y[0],
         1, x[1], y[1],
         1, x[2], y[2]}, 3, 3);

    M.invert_LU();

    for(size_t i = 0; i < N; ++i){
        this->a[i] = M(0,i);
        this->b[i] = M(1,i);
        this->c[i] = M(2,i);
    }

    gp_Vec v1(this->nodes[0]->point, this->nodes[1]->point);
    gp_Vec v2(this->nodes[0]->point, this->nodes[2]->point);

    this->delta = 0.5*(v1.Crossed(v2)).Magnitude();
}

math::Matrix BTRI3::get_K_ext(const math::Matrix & D, const gp_Pnt& center) const{
    math::Matrix K(K_DIM, K_DIM + 6);
    const auto& gsi = utils::GaussLegendreTri<ORDER>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const auto eps = this->get_eps(pi, center);
        const auto deps = this->get_deps();
        K += it->w*(deps.T()*D*eps);
    }
    K *= this->delta;

    return K;
}

math::Matrix BTRI3::get_stress_integrals(const math::Matrix& D, const gp_Pnt& center) const{
    // M_y M_x V_z M_z V_y V_x
    math::Matrix M(6, (K_DIM + 6), 0);
    math::Matrix E(S_SIZE, K_DIM + 6);
    const auto& gsi = utils::GaussLegendreTri<ORDER + 1>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const double dx = pi.X() - center.X();
        const double dy = pi.Y() - center.Y();
        const auto eps = this->get_eps(pi, center);
        E = D*eps;
        for(size_t i = 0; i < K_DIM + 6; ++i){
            M(0, i) += it->w*E(2, i)*dx;
            M(1, i) += it->w*E(2, i)*dy;
            M(2, i) += it->w*E(2, i);
            M(3, i) += it->w*(E(3, i)*dx - E(4, i)*dy);
            M(4, i) += it->w*E(3, i);
            M(5, i) += it->w*E(4, i);
        }
    }
    M *= this->delta;

    return M;
}

math::Matrix BTRI3::get_equilibrium_partial(const math::Matrix& D, const gp_Pnt& center, const std::vector<size_t>& stresses) const{
    const size_t KW = this->K_DIM; // workaround that's necessary for some reason
    math::Matrix K(KW, K_DIM + 6);
    const auto& gsi = utils::GaussLegendreTri<ORDER>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const auto eps = this->get_eps(pi, center);
        const auto deps = this->get_deps();
        const math::Matrix Ktmp(D*eps);
        for(size_t i = 0; i < K_DIM; ++i){
            for(size_t j = 0; j < K_DIM + 6; ++j){
                for(size_t k = 0; k < stresses.size(); ++k){
                    K(i,j) += it->w*(deps(stresses[k], i)*Ktmp(stresses[k], j));
                }
            }
        }
    }
    K *= this->delta;

    return K;
}

math::Matrix BTRI3::get_dz_vector_matrix(const math::Matrix& S, const math::Matrix& D, const gp_Pnt& center) const{
    const size_t KW = this->K_DIM; // workaround that's necessary for some reason
    math::Matrix M(KW, 2);
    const auto& gsi = utils::GaussLegendreTri<ORDER+1>::get();

    const math::Vector mult({D(4,2)/(S(2,2)*D(2,2)),
                             D(3,2)/(S(2,2)*D(2,2)),
                             1.0/S(2,2)});

    const std::vector<size_t> pos_j({0, 1});
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const double dx = pi.X() - center.X();
        const double dy = pi.Y() - center.Y();

        const math::Vector dxdy({dx, dy});
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            const double N = this->N(pi, i);
            const std::vector<size_t> seq({
                i*NODE_DOF + 0,
                i*NODE_DOF + 1,
                i*NODE_DOF + 2
            });
            M.slice(seq, pos_j) += -it->w*(mult*N*dxdy.T());
        }
    }
    M *= this->delta;

    return M;
}

math::Matrix BTRI3::get_dz_vector_matrix_1d(const math::Matrix& S, const math::Matrix& D, const math::Matrix& dS, const math::Matrix& dD, const gp_Pnt& center) const{
    const size_t KW = this->K_DIM; // workaround that's necessary for some reason
    math::Matrix M(KW, 2);
    const auto& gsi = utils::GaussLegendreTri<ORDER+1>::get();

    const double SD = S(2,2)*D(2,2);
    const double dSD = S(2,2)*dD(2,2) + dS(2,2)*D(2,2);
    const math::Vector mult({(dD(4,2)*SD - D(4,2)*dSD)/(SD*SD),
                             (dD(3,2)*SD - D(3,2)*dSD)/(SD*SD),
                             -dS(2,2)/(S(2,2)*S(2,2))});

    const std::vector<size_t> pos_j({0, 1});
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);
        const double dx = pi.X() - center.X();
        const double dy = pi.Y() - center.Y();

        const math::Vector dxdy({dx, dy});
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            const double N = this->N(pi, i);
            const std::vector<size_t> seq({
                i*NODE_DOF + 0,
                i*NODE_DOF + 1,
                i*NODE_DOF + 2
            });
            M.slice(seq, pos_j) += -it->w*(mult*N*dxdy.T());
        }
    }
    M *= this->delta;

    return M;
}

math::Vector BTRI3::get_dz_vector(const math::Matrix& S, const math::Matrix& D, const double Az, const double Bz, const gp_Pnt& center) const{
    const math::Vector ABz{Az, Bz};

    return this->get_dz_vector_matrix(S, D, center)*ABz;
}

math::Matrix BTRI3::get_force_vector_matrix(const math::Matrix& D, const gp_Pnt& center, const math::Matrix& rot) const{
    const size_t KW_P = this->parent->get_element_info()->get_k_dimension();

    const auto EPS_L = K_DIM + 6;

    const auto slice_i(math::slicer::sequence<size_t>(2, 5));
    const auto slice_j(math::slicer::sequence<size_t>(0, EPS_L));
    const math::Matrix P({0, 0, 1,
                          0, 1, 0,
                          1, 0, 0}, 3, 3);

    math::Matrix M(KW_P, EPS_L);

    const auto& gsi = utils::GaussLegendreTri<ORDER+2>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const gp_Pnt pi = this->GS_point(it->a, it->b, it->c);

        const math::Vector pv(rot*math::Vector{pi.X(), pi.Y(), pi.Z()});
        const gp_Pnt pr(pv[0], pv[1], pv[2]);

        const auto eps = this->get_eps(pi, center);

        const auto S(D*eps);

        const auto N = this->parent->get_Ni(pr);
        M += it->w*(N.T()*rot*P*S.slice(slice_i, slice_j));
    }
    M *= this->delta;

    return M;
}

math::Vector BTRI3::get_force_vector(const math::Matrix& D, const std::vector<double>& u, const gp_Pnt& center, const math::Matrix& rot, const bool transpose) const{

    const auto M(this->get_force_vector_matrix(D, center, rot));

    if(transpose){
        const size_t parent_node_num = this->parent->get_element_info()->get_nodes_per_element();
        std::vector<long> pos(parent_node_num*NODE_DOF);
        math::slicer::from_node_upos(this->parent->nodes, parent_node_num, NODE_DOF, pos);
        const math::VectorSliceGeneralView us(u, pos);
        return M.T()*us;
    } else {
        std::vector<size_t> pos(K_DIM + 6);
        math::slicer::from_node_id(this->nodes, NODES_PER_ELEM, NODE_DOF, pos);
        std::iota(pos.begin()+K_DIM, pos.end(), u.size() - 6);
        const math::VectorSliceView us(u, pos);
        return M*us;
    }
}


math::Matrix BTRI3::diffusion_1dof(const math::Matrix& A) const{
    const auto B = this->dN_mat_1dof();
    return delta*B.T()*A*B;
}
math::Matrix BTRI3::advection_1dof(const math::Vector& v) const{
    const auto NN = this->N_mat_1dof(this->GS_point(1.0/3.0, 1.0/3.0, 1.0/3.0));
    const auto B = this->dN_mat_1dof();
    return this->delta*NN*(v.T()*B);
}
math::Matrix BTRI3::absorption_1dof() const{
    const auto& gsi = utils::GaussLegendreTri<2>::get();
    math::Matrix result(3, 3);
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        const auto NN = this->N_mat_1dof(this->GS_point(it->a, it->b, it->c));
        result += it->w*NN*NN.T();
    }
    result *= delta;

    return result;
}
math::Vector BTRI3::source_1dof() const{
    return this->delta*this->N_mat_1dof(this->GS_point(1.0/3.0, 1.0/3.0, 1.0/3.0));
}
math::Vector BTRI3::flow_1dof(const std::array<const Node*, 2>& nodes) const{
    gp_Pnt p = nodes[0]->point;
    p.BaryCenter(1, nodes[1]->point, 1);
    const double d = nodes[0]->point.Distance(nodes[1]->point);

    return d*this->N_mat_1dof(p);
}

}
