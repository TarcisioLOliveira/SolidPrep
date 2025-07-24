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

#include "element/H8.hpp"
#include "boundary_element/BQ4.hpp"
#include "logger.hpp"
#include "math/matrix.hpp"
#include "utils/gauss_legendre.hpp"
#include "project_specification/registry.hpp"

namespace element{

const bool H8::reg = projspec::ElementRegistry::add("H8",
    std::make_unique<MeshElementFactoryImpl<H8>>());

H8::H8(ElementShape s):
    MeshElementCommon3DHex<H8>(s.nodes){
}

std::unique_ptr<BoundaryMeshElementFactory> H8::get_boundary_element_info() {
    return std::unique_ptr<BoundaryMeshElementFactory>(new BoundaryMeshElementFactoryImpl<boundary_element::BQ4>());
}
std::unique_ptr<ContactMeshElementFactory> H8::get_contact_element_info() {
    logger::log_assert(false, logger::ERROR, "CQ4 ELEMENT TYPE NOT IMPLEMENTED");
    return std::unique_ptr<ContactMeshElementFactory>();
    //return std::unique_ptr<ContactMeshElementFactory>(new ContactMeshElementFactoryImpl<contact_element::CQ4>());
}

math::Matrix H8::get_k(const math::Matrix & D, const double t) const{
    (void)t;
    constexpr size_t N = H8::NODES_PER_ELEM;
    std::array<double, N> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X() - c.X();
        y[i] = this->nodes[i]->point.Y() - c.Y();
        z[i] = this->nodes[i]->point.Z() - c.Z();
    }

    math::Matrix k(K_DIM, K_DIM);
    constexpr size_t GN = 3;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            for(auto zeta = GL.begin(); zeta < GL.end(); ++zeta){
                const auto B = this->B_mat_norm(xi->x, eta->x, zeta->x, x, y, z);
                const double detJ = this->J(xi->x, eta->x, zeta->x, x, y, z).determinant();

                k += B.T()*((xi->w*eta->w*zeta->w*detJ)*D)*B;
            }
        }
    }

    return k;
}

math::Matrix H8::get_nodal_density_gradient(gp_Pnt p) const{
    constexpr size_t N = H8::NODES_PER_ELEM;
    std::array<double, N> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X() - c.X();
        y[i] = this->nodes[i]->point.Y() - c.Y();
        z[i] = this->nodes[i]->point.Z() - c.Z();
    }

    return this->dN_mat_1dof(p.X(), p.Y(), p.Z(), x, y, z);
}

math::Matrix H8::get_B(const gp_Pnt& point) const{
    constexpr size_t N = H8::NODES_PER_ELEM;
    std::array<double, N> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X() - c.X();
        y[i] = this->nodes[i]->point.Y() - c.Y();
        z[i] = this->nodes[i]->point.Z() - c.Z();
    }
    const double px = point.X() - c.X();
    const double py = point.Y() - c.Y();
    const double pz = point.Z() - c.Z();

    return this->B_mat_norm(px, py, pz, x, y, z);
}

math::Matrix H8::get_Nf(const double t, const std::vector<gp_Pnt>& points) const{
    (void)t;

    std::array<double, BOUNDARY_NODES_PER_ELEM> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
        x[i] = points[i].X() - c.X();
        y[i] = points[i].Y() - c.Y();
        z[i] = points[i].Z() - c.Z();
    }
    std::array<double, NODES_PER_ELEM> x2, y2, z2;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        x2[i] = this->nodes[i]->point.X() - c.X();
        y2[i] = this->nodes[i]->point.Y() - c.Y();
        z2[i] = this->nodes[i]->point.Z() - c.Z();
    }

    const CubeSide cs = this->get_cube_side(points);

    math::Matrix Nf(NODE_DOF, K_DIM);
    constexpr size_t GN = 2;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto drnorm = this->surface_drnorm(xi->x, eta->x, x, y, z);
            gp_Pnt pi = this->to_surface_point(xi->x, eta->x, cs);
            Nf += (xi->w*eta->w*drnorm)*N_mat_norm(pi.X(), pi.Y(), pi.Z());
        }
    }

    return Nf;
}

math::Matrix H8::get_R(const math::Matrix & K, const double t, const std::vector<gp_Pnt>& points) const{
    (void)t;

    std::array<double, BOUNDARY_NODES_PER_ELEM> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
        x[i] = points[i].X() - c.X();
        y[i] = points[i].Y() - c.Y();
        z[i] = points[i].Z() - c.Z();
    }
    std::array<double, NODES_PER_ELEM> x2, y2, z2;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        x2[i] = this->nodes[i]->point.X() - c.X();
        y2[i] = this->nodes[i]->point.Y() - c.Y();
        z2[i] = this->nodes[i]->point.Z() - c.Z();
    }

    const CubeSide cs = this->get_cube_side(points);

    math::Matrix R(K_DIM, K_DIM);
    constexpr size_t GN = 5;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto drnorm = this->surface_drnorm(xi->x, eta->x, x, y, z);
            gp_Pnt pi = this->to_surface_point(xi->x, eta->x, cs);
            const auto NN = N_mat_norm(pi.X(), pi.Y(), pi.Z());
            R += (xi->w*eta->w*drnorm)*NN.T()*K*NN;
        }
    }

    return R;
}

math::Matrix H8::diffusion_1dof(const double t, const math::Matrix& A) const{
    (void)t;
    constexpr size_t N = H8::NODES_PER_ELEM;
    std::array<double, N> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X() - c.X();
        y[i] = this->nodes[i]->point.Y() - c.Y();
        z[i] = this->nodes[i]->point.Z() - c.Z();
    }

    math::Matrix M(NODES_PER_ELEM, NODES_PER_ELEM);
    constexpr size_t GN = 4;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            for(auto zeta = GL.begin(); zeta < GL.end(); ++zeta){
                const auto B = this->dN_mat_1dof(xi->x, eta->x, zeta->x, x, y, z);
                const double detJ = this->J(xi->x, eta->x, zeta->x, x, y, z).determinant();

                M += B.T()*((xi->w*eta->w*zeta->w*detJ)*A)*B;
            }
        }
    }

    return M;
}
math::Matrix H8::advection_1dof(const double t, const math::Vector& v) const{
    (void)t;
    constexpr size_t N = H8::NODES_PER_ELEM;
    std::array<double, N> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X() - c.X();
        y[i] = this->nodes[i]->point.Y() - c.Y();
        z[i] = this->nodes[i]->point.Z() - c.Z();
    }

    math::Matrix M(NODES_PER_ELEM, NODES_PER_ELEM);
    constexpr size_t GN = 5;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            for(auto zeta = GL.begin(); zeta < GL.end(); ++zeta){
                const auto B = this->dN_mat_1dof(xi->x, eta->x, zeta->x, x, y, z);
                const auto Nv = this->N_mat_1dof(xi->x, eta->x, zeta->x);
                const double detJ = this->J(xi->x, eta->x, zeta->x, x, y, z).determinant();

                M += B.T()*((xi->w*eta->w*zeta->w*detJ)*v)*Nv.T();
            }
        }
    }

    return M;

}
math::Matrix H8::absorption_1dof(const double t) const{
    (void)t;
    constexpr size_t N = H8::NODES_PER_ELEM;
    std::array<double, N> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X() - c.X();
        y[i] = this->nodes[i]->point.Y() - c.Y();
        z[i] = this->nodes[i]->point.Z() - c.Z();
    }

    math::Matrix M(NODES_PER_ELEM, NODES_PER_ELEM);
    constexpr size_t GN = 6;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            for(auto zeta = GL.begin(); zeta < GL.end(); ++zeta){
                const auto Nv = this->N_mat_1dof(xi->x, eta->x, zeta->x);
                const double detJ = this->J(xi->x, eta->x, zeta->x, x, y, z).determinant();

                M += Nv*((xi->w*eta->w*zeta->w*detJ))*Nv.T();
            }
        }
    }

    return M;
}
math::Vector H8::source_1dof(const double t) const{
    (void)t;
    constexpr size_t N = H8::NODES_PER_ELEM;
    std::array<double, N> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < N; ++i){
        x[i] = this->nodes[i]->point.X() - c.X();
        y[i] = this->nodes[i]->point.Y() - c.Y();
        z[i] = this->nodes[i]->point.Z() - c.Z();
    }

    math::Vector M(NODES_PER_ELEM);
    M.fill(0);
    constexpr size_t GN = 2;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            for(auto zeta = GL.begin(); zeta < GL.end(); ++zeta){
                const auto Nv = this->N_mat_1dof(xi->x, eta->x, zeta->x);
                const double detJ = this->J(xi->x, eta->x, zeta->x, x, y, z).determinant();

                M += (xi->w*eta->w*zeta->w*detJ)*Nv;
            }
        }
    }

    return M;
}

math::Vector H8::flow_1dof(const double t, const MeshNode** nodes) const{
    (void)t;

    std::vector<gp_Pnt> points(BOUNDARY_NODES_PER_ELEM);
    std::array<double, BOUNDARY_NODES_PER_ELEM> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
        x[i] = nodes[i]->point.X() - c.X();
        y[i] = nodes[i]->point.Y() - c.Y();
        z[i] = nodes[i]->point.Z() - c.Z();
        points[i] = nodes[i]->point;
    }
    std::array<double, NODES_PER_ELEM> x2, y2, z2;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        x2[i] = this->nodes[i]->point.X() - c.X();
        y2[i] = this->nodes[i]->point.Y() - c.Y();
        z2[i] = this->nodes[i]->point.Z() - c.Z();
    }

    const CubeSide cs = this->get_cube_side(points);

    math::Vector M(NODES_PER_ELEM);
    constexpr size_t GN = 2;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto drnorm = this->surface_drnorm(xi->x, eta->x, x, y, z);
            gp_Pnt pi = this->to_surface_point(xi->x, eta->x, cs);
            M += (xi->w*eta->w*drnorm)*N_mat_1dof(pi.X(), pi.Y(), pi.Z());
        }
    }

    return M;
}

math::Matrix H8::get_uu(const MeshElement* const e2, const std::vector<gp_Pnt>& bounds, const gp_Dir n) const{
    const size_t KW = K_DIM;

    std::array<double, BOUNDARY_NODES_PER_ELEM> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
        x[i] = bounds[i].X() - c.X();
        y[i] = bounds[i].Y() - c.Y();
        z[i] = bounds[i].Z() - c.Z();
    }

    math::Vector NN(2*KW, 0);
    math::Matrix MnMn(2*KW, 2*KW, 0);

    const CubeSide cs = this->get_cube_side(bounds);

    constexpr size_t GN = 2*ORDER;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){
            const auto drnorm = this->surface_drnorm(xi->x, eta->x, x, y, z);
            const gp_Pnt pi = this->to_surface_point(xi->x, eta->x, cs);
            const auto N1 = this->get_Ni(pi);
            const auto N2 = e2->get_Ni(pi);
            NN.fill(0);
            for(size_t i = 0; i < KW; ++i){
                for(size_t j = 0; j < DIM; ++j){
                    NN[i] += N1(j, i)*n.Coord(1+j);
                    NN[i + KW] += N2(j, i)*n.Coord(1+j);
                }
            }
            MnMn += (xi->w*eta->w*drnorm)*(NN*NN.T());
        }
    }
    for(size_t i = 0; i < KW; ++i){
        for(size_t j = KW; j < 2*KW; ++j){
                MnMn(i, j) *= -1.0;
                MnMn(j, i) *= -1.0;
        }
    }

    return MnMn;
}
math::Matrix H8::get_MnMn(const MeshElement* const e2, const std::vector<double>& u_ext, const std::vector<gp_Pnt>& bounds, const gp_Dir n) const{
    const size_t DOF = NODE_DOF;
    const size_t KW = K_DIM;

    std::array<double, BOUNDARY_NODES_PER_ELEM> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
        x[i] = bounds[i].X() - c.X();
        y[i] = bounds[i].Y() - c.Y();
        z[i] = bounds[i].Z() - c.Z();
    }

    math::Vector uv1(KW), uv2(KW);
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        for(size_t j = 0; j < DOF; ++j){
            uv1[i*DOF + j] = u_ext[this->nodes[i]->u_pos[j]];
            uv2[i*DOF + j] = u_ext[e2->nodes[i]->u_pos[j]];
        }
    }
    math::Vector NN(2*KW, 0);
    math::Matrix MnMn(2*KW, 2*KW, 0);

    const CubeSide cs = this->get_cube_side(bounds);

    constexpr size_t GN = 2*ORDER + 2;
    const auto& GL = utils::GaussLegendre<GN>::get();

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){

            const auto drnorm = this->surface_drnorm(xi->x, eta->x, x, y, z);
            const gp_Pnt pi = this->to_surface_point(xi->x, eta->x, cs);
            const auto N1 = this->get_Ni(pi);
            const auto N2 = e2->get_Ni(pi);
            double gp = 0;
            math::Vector up1(N1*uv1);
            math::Vector up2(N2*uv2);
            for(size_t j = 0; j < DIM; ++j){
                gp += (up2[j] - up1[j])*n.Coord(1+j);
            }
            if(gp < 1e-7){
                NN.fill(0);
                for(size_t i = 0; i < KW; ++i){
                    for(size_t j = 0; j < DIM; ++j){
                        NN[i] += N1(j, i)*n.Coord(1+j);
                        NN[i + KW] += N2(j, i)*n.Coord(1+j);
                    }
                }
                MnMn += (xi->w*eta->w*drnorm)*(NN*NN.T());
            }
        }
    }
    for(size_t i = 0; i < KW; ++i){
        for(size_t j = KW; j < 2*KW; ++j){
                MnMn(i, j) *= -1.0;
                MnMn(j, i) *= -1.0;
        }
    }

    return MnMn;
}
math::Matrix H8::get_MnMn_log(const MeshElement* const e2, const std::vector<double>& u_ext, const double eps, const std::vector<gp_Pnt>& bounds, const gp_Dir n) const{
    const size_t DOF = NODE_DOF;
    const size_t KW = K_DIM;

    std::array<double, BOUNDARY_NODES_PER_ELEM> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
        x[i] = bounds[i].X() - c.X();
        y[i] = bounds[i].Y() - c.Y();
        z[i] = bounds[i].Z() - c.Z();
    }

    math::Vector uv1(KW), uv2(KW);
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        for(size_t j = 0; j < DOF; ++j){
            uv1[i*DOF + j] = u_ext[this->nodes[i]->u_pos[j]];
            uv2[i*DOF + j] = u_ext[e2->nodes[i]->u_pos[j]];
        }
    }
    math::Vector NN(2*KW, 0);
    math::Matrix MnMn(2*KW, 2*KW, 0);

    constexpr size_t GN = 2*ORDER + 3;
    const auto& GL = utils::GaussLegendre<GN>::get();

    const CubeSide cs = this->get_cube_side(bounds);

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){

            const auto drnorm = this->surface_drnorm(xi->x, eta->x, x, y, z);
            const gp_Pnt pi = this->to_surface_point(xi->x, eta->x, cs);

            const auto N1 = this->get_Ni(pi);
            const auto N2 = e2->get_Ni(pi);
            double gp = 0;
            math::Vector up1(N1*uv1);
            math::Vector up2(N2*uv2);
            for(size_t j = 0; j < DIM; ++j){
                gp += (up2[j] - up1[j])*n.Coord(1+j);
            }
            const double tmp = (gp + eps);
            const double log = 1.0/(tmp*tmp);
            NN.fill(0);
            for(size_t i = 0; i < KW; ++i){
                for(size_t j = 0; j < DIM; ++j){
                    NN[i] -= N1(j, i)*n.Coord(1+j);
                    NN[i + KW] += N2(j, i)*n.Coord(1+j);
                }
            }
            MnMn += ((xi->w*eta->w*drnorm*log)*NN)*NN.T();
        }
    }
    //for(size_t i = 0; i < KW; ++i){
    //    for(size_t j = KW; j < 2*KW; ++j){
    //        MnMn(i, j) *= -1.0;
    //        MnMn(j, i) *= -1.0;
    //    }
    //}

    return MnMn;
}
void H8::Ku_log(const MeshElement* const e2, const std::vector<double>& u_ext, const double eps, const double C, const std::vector<gp_Pnt>& bounds, const gp_Dir n, std::vector<double>& Ku) const{
    const size_t DOF = NODE_DOF;
    const size_t KW = K_DIM;

    std::array<double, BOUNDARY_NODES_PER_ELEM> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
        x[i] = bounds[i].X() - c.X();
        y[i] = bounds[i].Y() - c.Y();
        z[i] = bounds[i].Z() - c.Z();
    }

    math::Vector uv1(KW), uv2(KW);
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        for(size_t j = 0; j < DOF; ++j){
            uv1[i*DOF + j] = u_ext[this->nodes[i]->u_pos[j]];
            uv2[i*DOF + j] = u_ext[e2->nodes[i]->u_pos[j]];
        }
    }

    math::Vector NN(2*KW, 0);
    math::Vector Mn(2*KW, 0);

    constexpr size_t GN = 2*ORDER + 3;
    const auto& GL = utils::GaussLegendre<GN>::get();

    const CubeSide cs = this->get_cube_side(bounds);

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){

            const auto drnorm = this->surface_drnorm(xi->x, eta->x, x, y, z);
            const gp_Pnt pi = this->to_surface_point(xi->x, eta->x, cs);
            const auto N1 = this->get_Ni(pi);
            const auto N2 = e2->get_Ni(pi);
            double gp = 0;
            math::Vector up1(N1*uv1);
            math::Vector up2(N2*uv2);
            for(size_t j = 0; j < DIM; ++j){
                gp += (up2[j] - up1[j])*n.Coord(1+j);
            }
            const double log = -1.0/(gp + eps);
            NN.fill(0);
            for(size_t i = 0; i < KW; ++i){
                for(size_t j = 0; j < DIM; ++j){
                    NN[i] -= N1(j, i)*n.Coord(1+j);
                    NN[i + KW] += N2(j, i)*n.Coord(1+j);
                }
            }
            Mn += (xi->w*eta->w*drnorm*log)*NN;
        }
    }
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        for(size_t j = 0; j < DOF; ++j){
            Ku[this->nodes[i]->u_pos[j]] = C*Mn[i*DOF + j];
            Ku[e2->nodes[i]->u_pos[j]] = C*Mn[KW + i*DOF + j];;
        }
    }
}
void H8::dKu_log(const MeshElement* const e2, const std::vector<double>& u_ext, const std::vector<double>& du, const double eps, const double C, const std::vector<gp_Pnt>& bounds, const gp_Dir n, std::vector<double>& dKu) const{
    const size_t DOF = NODE_DOF;
    const size_t KW = K_DIM;

    std::array<double, BOUNDARY_NODES_PER_ELEM> x, y, z;
    const gp_Pnt c(this->get_centroid());
    for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
        x[i] = bounds[i].X() - c.X();
        y[i] = bounds[i].Y() - c.Y();
        z[i] = bounds[i].Z() - c.Z();
    }

    math::Vector uv1(KW), uv2(KW);
    math::Vector duv1(KW), duv2(KW);
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        for(size_t j = 0; j < DOF; ++j){
            uv1[i*DOF + j] = u_ext[this->nodes[i]->u_pos[j]];
            uv2[i*DOF + j] = u_ext[e2->nodes[i]->u_pos[j]];
            duv1[i*DOF + j] = du[this->nodes[i]->u_pos[j]];
            duv2[i*DOF + j] = du[e2->nodes[i]->u_pos[j]];
        }
    }
    math::Vector NN(2*KW, 0);
    math::Vector Mn(2*KW, 0);

    constexpr size_t GN = 2*ORDER + 3;
    const auto& GL = utils::GaussLegendre<GN>::get();

    const CubeSide cs = this->get_cube_side(bounds);

    for(auto xi = GL.begin(); xi < GL.end(); ++xi){
        for(auto eta = GL.begin(); eta < GL.end(); ++eta){

            const auto drnorm = this->surface_drnorm(xi->x, eta->x, x, y, z);
            const gp_Pnt pi = this->to_surface_point(xi->x, eta->x, cs);
            const auto N1 = this->get_Ni(pi);
            const auto N2 = e2->get_Ni(pi);
            double gp = 0;
            double dgp = 0;
            math::Vector up1(N1*uv1);
            math::Vector up2(N2*uv2);
            math::Vector dup1(N1*duv1);
            math::Vector dup2(N2*duv2);
            for(size_t j = 0; j < DIM; ++j){
                gp += (up2[j] - up1[j])*n.Coord(1+j);
                dgp += (dup2[j] - dup1[j])*n.Coord(1+j);
            }
            //const double tmp = eps/(gp + eps);
            //const double log = -10*eps*tmp;
            const double s = (gp + eps);
            const double dlog = dgp/(s*s);
            NN.fill(0);
            for(size_t i = 0; i < KW; ++i){
                for(size_t j = 0; j < DIM; ++j){
                    NN[i] -= N1(j, i)*n.Coord(1+j);
                    NN[i + KW] += N2(j, i)*n.Coord(1+j);
                }
            }
            Mn += (xi->w*eta->w*drnorm*dlog)*NN;
        }
    }
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        for(size_t j = 0; j < DOF; ++j){
            dKu[this->nodes[i]->u_pos[j]] = C*Mn[i*DOF + j];
            dKu[e2->nodes[i]->u_pos[j]] = C*Mn[KW + i*DOF + j];;
        }
    }
}

H8::CubeSide H8::get_cube_side(const std::vector<gp_Pnt>& points) const{
    // As I was unable to make this work using natural coordinates, this method
    // has additional code to find which face of the hex the loads are being
    // applied onto.
    //
    // This is probably (hopefully) not the best solution, but it works for now.

    int xp[8] = {-1, 1, 1, -1, -1, 1, 1, -1};
    int yp[8] = {1, 1, -1, -1, 1, 1, -1, -1};
    int zp[8] = {-1, -1, -1, -1, 1, 1, 1, 1};

    int idx[4];
    int idx_i = 0;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        for(size_t j = 0; j < BOUNDARY_NODES_PER_ELEM; ++j){
            if(this->nodes[i]->point.IsEqual(points[j], Precision::Confusion())){
                idx[idx_i] = i;
                ++idx_i;
                break;
            }
        }
    }
    int xt = xp[idx[0]];
    int yt = yp[idx[0]]; 
    int zt = zp[idx[0]]; 
    for(size_t i = 1; i < 4; ++i){
        if(xp[idx[i]] != xt){
            xt = 0;
        }
        if(yp[idx[i]] != yt){
            yt = 0;
        }
        if(zp[idx[i]] != zt){
            zt = 0;
        }
    }
    logger::log_assert(xt*yt*zt == 0 && xt*yt == 0 && xt*zt == 0 && yt*zt == 0,
                       logger::ERROR,
                       "wrong assumptions getting cube side: {} {} {}", xt, yt, zt);

    if(xt == 1){
        return CubeSide::X_MAX;
    } else if(xt == -1){
        return CubeSide::X_MIN;
    }
    if(yt == 1){
        return CubeSide::Y_MAX;
    } else if(yt == -1){
        return CubeSide::Y_MIN;
    }
    if(zt == 1){
        return CubeSide::Z_MAX;
    } else if(zt == -1){
        return CubeSide::Z_MIN;
    }

    return CubeSide::UNKNOWN;
}
gp_Pnt H8::to_surface_point(const double xi, const double eta, const H8::CubeSide side) const{
    switch(side){
        case CubeSide::X_MAX:
            return {1, xi, eta};
        case CubeSide::X_MIN:
            return {-1, xi, eta};
        case CubeSide::Y_MAX:
            return {xi, 1, eta};
        case CubeSide::Y_MIN:
            return {xi, -1, eta};
        case CubeSide::Z_MAX:
            return {xi, eta, 1};
        case CubeSide::Z_MIN:
            return {xi, eta, -1};
        case CubeSide::UNKNOWN:
            logger::log_assert(false, logger::ERROR,  "unknown cube side for H8");
    }
    return {0,0,0};
}

math::Matrix H8::get_Ni(const gp_Pnt& p) const{
    return this->N_mat_norm(p.X(), p.Y(), p.Z());
}

}
