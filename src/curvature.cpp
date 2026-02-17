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

#include "curvature.hpp"
#include <Eigen/Dense>
#include <vector>
#include "general_solver/mumps_general.hpp"
#include "logger.hpp"
#include "math/matrix.hpp"
#include "utils/gauss_legendre.hpp"
#include "utils/basis_tensor.hpp"
#include "math/vector_view.hpp"
#include "math/slicer.hpp"

Curvature::Curvature(const Material* mat, math::Matrix rot2D, math::Matrix rot3D, std::unique_ptr<BoundaryMeshElementFactory> elem_info, double V_u, double V_v, double V_w, double M_u, double M_v, double M_w):
    mat(mat), 
    rot2D(std::move(rot2D)), 
    //Lek_basis(
    //    {0, 0, -1,
    //     0, 1,  0,
    //     1, 0,  0}, 3, 3),
    rot3D(std::move(rot3D)), 
    elem_info(std::move(elem_info)),
    V_u(V_w), V_v(V_v), V_w(-V_u),
    M_u(M_w), 
    M_v(M_v),
    M_w(-M_u)
{

    //this->permute_shear_3D = math::Matrix(
    //    {1, 0, 0, 0, 0, 0,
    //     0, 1, 0, 0, 0, 0,
    //     0, 0, 1, 0, 0, 0,
    //     0, 0, 0, 0, 0, 1,
    //     0, 0, 0, 0, 1, 0,
    //     0, 0, 0, 1, 0, 0}, 6, 6);

    //this->rotate_X_Z_3D = utils::basis_tensor_3D(Lek_basis);
    //this->rotate_X_Z_3D_inv = utils::basis_tensor_3D_inv_T(Lek_basis);

    //this->permute_and_rotate_3D = permute_shear_3D*rotate_X_Z_3D;
    //this->permute_and_rotate_3D_inv = rotate_X_Z_3D_inv*permute_shear_3D;
}

void Curvature::generate_curvature_3D(const std::vector<std::unique_ptr<MeshNode>>& boundary_nodes, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, size_t reduced_vector_size, size_t full_vector_size)
{
    (void) boundary_nodes;
    std::array<gp_Pnt, 3> points;

    const auto fn_A =
        [this](const math::Matrix& S, const gp_Pnt& px)->double{
            return this->make_A_base_3D(S, px);
        };
    const auto fn_EA =
        [this](const math::Matrix& S, const gp_Pnt& px)->double{
            return this->make_EA_base_3D(S, px);
        };
    const auto fn_EA_u =
        [this](const math::Matrix& S, const gp_Pnt& px)->double{
            return this->make_EA_u_base_3D(S, px);
        };
    const auto fn_EA_v =
        [this](const math::Matrix& S, const gp_Pnt& px)->double{
            return this->make_EA_v_base_3D(S, px);
        };
    const auto fn_EA_w =
        [this](const math::Matrix& S, const gp_Pnt& px)->double{
            return this->make_EA_w_base_3D(S, px);
        };
    const auto fn_EI_vv =
        [this](const math::Matrix& S, const gp_Pnt& px)->double{
            return this->make_EI_vv_base_3D(S, px);
        };
    const auto fn_EI_uu =
        [this](const math::Matrix& S, const gp_Pnt& px)->double{
            return this->make_EI_uu_base_3D(S, px);
        };
    const auto fn_EI_uv =
        [this](const math::Matrix& S, const gp_Pnt& px)->double{
            return this->make_EI_uv_base_3D(S, px);
        };

    const auto result1 = this->integrate_surface_3D(boundary_mesh, {fn_A, fn_EA, fn_EA_u, fn_EA_v, fn_EA_w});
    this->Area = result1[0];
    this->EA   = result1[1];
    this->EA_u = result1[2];
    this->EA_v = result1[3];
    const double EA_w = result1[4];
    this->c_u = EA_u/this->EA;
    this->c_v = EA_v/this->EA;
    this->c_w = EA_w/this->EA;
    this->center = gp_Pnt(c_u, c_v, c_w);
    const auto result2 =
        this->integrate_surface_3D(boundary_mesh, 
                {fn_EI_uu,
                fn_EI_vv,
                fn_EI_uv});
    this->EI_uu = result2[0];
    this->EI_vv = result2[1];
    this->EI_uv = result2[2];
    logger::quick_log("A: ", Area);
    logger::quick_log("EA: ", EA);
    logger::quick_log("c_v: ", c_v);
    logger::quick_log("c_u: ", c_u);
    logger::quick_log("EA_u", EA_u);
    logger::quick_log("EA_v", EA_v);
    logger::quick_log("EI_uu", EI_uu);
    logger::quick_log("EI_vv", EI_vv);
    logger::quick_log("EI_uv", EI_uv);

    this->reduced_vec_len = reduced_vector_size;
    this->full_vector_size = full_vector_size;

    this->u.resize(3*full_vector_size + 6, 0); 

    this->calculate_stress_field_3D(boundary_mesh);
}

math::Vector Curvature::get_force_vector_3D(const BoundaryMeshElement* const e) const{
    gp_Pnt c = e->get_centroid();
    gp_Pnt c2 = utils::change_point(c, this->rot3D);

    const auto D = this->get_D_3D(e->parent, c2);
    const auto S = this->get_S_3D(e->parent, c2);
    const auto DT = this->get_DT_3D(S, D);

    return e->get_force_vector(DT, this->u, this->center, this->rot3D);
}

double Curvature::get_adjoint_gradient(const std::vector<double>& fea_u, const math::Matrix& S, const math::Matrix& D, const math::Matrix& dS, const math::Matrix& dD, const BoundaryMeshElement* const e) const{
    double result = 0;

    const auto parent_elem_info = e->parent->get_element_info();

    const size_t K_DIM = this->elem_info->get_k_dimension();
    const size_t NODE_DOF = this->elem_info->get_dof_per_node();
    const size_t NODES_PER_ELEM = this->elem_info->get_nodes_per_element();

    const size_t PARENT_K_DIM = parent_elem_info->get_k_dimension();
    const size_t PARENT_NODE_DOF = parent_elem_info->get_dof_per_node();
    const size_t PARENT_NODES_PER_ELEM = parent_elem_info->get_nodes_per_element();

    const auto dDT = this->get_DT_3D_1d(S, D, dS, dD);
    const auto dTDT = this->get_TDT_3D_1d(S, D, dS, dD);

    const auto df = e->get_force_vector(dDT, this->u, this->center, this->rot3D);

    std::vector<size_t> posup(PARENT_K_DIM);
    math::slicer::from_node_id(e->parent->nodes, PARENT_NODES_PER_ELEM, PARENT_NODE_DOF, posup);
    math::VectorSliceView fea_s(fea_u, posup);
    result += df.dot(fea_s);

    const auto St = e->get_stress_integrals(dDT, this->center);
    const auto tau_xz_dz = e->get_equilibrium_partial(dTDT, this->center, {0, 5});
    const auto tau_yz_dz = e->get_equilibrium_partial(dTDT, this->center, {1, 5});
    const auto sigma_z_dz = e->get_equilibrium_partial(dTDT, this->center, {3, 4});

    const size_t u_offset = this->u.size() - 6;
    std::vector<size_t> posi(K_DIM);
    std::vector<size_t> posj(K_DIM + 6);
    std::vector<size_t> pos_last(6);

    math::slicer::from_node_id(e->nodes, NODES_PER_ELEM, NODE_DOF, posi);
    std::copy(posi.begin(), posi.end(), posj.begin());
    for(size_t i = 0; i < 6; ++i){
        posj[K_DIM + i] = u_offset + i;
        pos_last[i] = u_offset + i;
    }

    math::VectorSliceView u_s(this->u, posj);
    math::VectorSliceView u_bar_s(u_bar, posi);
    math::VectorSliceView u_bar_end(u_bar, pos_last);

    result -= ((tau_xz_dz + tau_yz_dz + sigma_z_dz)*u_s).dot(u_bar_s);
    result -= ((St)*u_s).dot(u_bar_end);

    const auto Mz = e->get_dz_vector_matrix_1d(S, D, dS, dD, this->center);

    result += (Mz*this->ABz).dot(u_bar_s);

    const auto fn_dEI_vv =
        [this](const math::Matrix& S, const math::Matrix& dS, const gp_Pnt& px)->double{
            return this->make_EI_vv_base_3D_1d(S, dS, px);
        };
    const auto fn_dEI_uu =
        [this](const math::Matrix& S, const math::Matrix& dS, const gp_Pnt& px)->double{
            return this->make_EI_uu_base_3D_1d(S, dS, px);
        };
    const auto fn_dEI_uv =
        [this](const math::Matrix& S, const math::Matrix& dS, const gp_Pnt& px)->double{
            return this->make_EI_uv_base_3D_1d(S, dS, px);
        };

    const std::vector<std::function<double(const math::Matrix& S, const math::Matrix& dS, const gp_Pnt& px)>> fn
        ({fn_dEI_vv, fn_dEI_uu, fn_dEI_uv});

    std::vector<double> dEI(fn.size(), 0);
    if(this->elem_info->get_shape_type() == Element::Shape::TRI){
        constexpr size_t N = 3;
        math::Matrix points(3, N);
        std::array<gp_Pnt, N> rel_points;
        for(size_t x = 0; x < N; ++x){
            rel_points[x] = e->nodes[x]->point;
            for(size_t y = 0; y < 3; ++y){
                points(y, x) = rel_points[x].Coord(y+1);
            }
        }
        math::Matrix rotd_p(this->rot3D*points);
        std::array<gp_Pnt, N> abs_points{
            gp_Pnt(rotd_p(0, 0), rotd_p(1, 0), rotd_p(2, 0)),
            gp_Pnt(rotd_p(0, 1), rotd_p(1, 1), rotd_p(2, 1)),
            gp_Pnt(rotd_p(0, 2), rotd_p(1, 2), rotd_p(2, 2))
        };

        const auto& p = abs_points;
        const auto& px = rel_points;

        gp_Vec v1(p[1], p[0]);
        gp_Vec v2(p[2], p[0]);
        const gp_Pnt centroid(
                (p[0].X() + p[1].X() + p[2].X())/3.0,
                (p[0].Y() + p[1].Y() + p[2].Y())/3.0,
                (p[0].Z() + p[1].Z() + p[2].Z())/3.0);
        const double drnorm = (v1.Crossed(v2)).Magnitude()/2;

        for(size_t ri = 0; ri < fn.size(); ++ri){
            const auto& gsi = utils::GaussLegendreTri<2>::get();
            for(auto it = gsi.begin(); it < gsi.end(); ++it){
                for(size_t i = 0; i < fn.size(); ++i){
                    gp_Pnt pxi(this->GLT_point(*it, px));
                    dEI[ri] += it->w*fn[i](S, dS, pxi);
                }
            }
            dEI[ri] *= drnorm;
        }

    } else if(this->elem_info->get_shape_type() == Element::Shape::QUAD){
        constexpr size_t N = 4;
        math::Matrix points(3, N);
        std::array<gp_Pnt, N> rel_points;
        for(size_t x = 0; x < N; ++x){
            rel_points[x] = e->nodes[x]->point;
            for(size_t y = 0; y < 3; ++y){
                points(y, x) = rel_points[x].Coord(y+1);
            }
        }
        math::Matrix rotd_p(this->rot3D*points);
        std::array<gp_Pnt, N> abs_points{
            gp_Pnt(rotd_p(0, 0), rotd_p(1, 0), rotd_p(2, 0)),
            gp_Pnt(rotd_p(0, 1), rotd_p(1, 1), rotd_p(2, 1)),
            gp_Pnt(rotd_p(0, 2), rotd_p(1, 2), rotd_p(2, 2)),
            gp_Pnt(rotd_p(0, 3), rotd_p(1, 3), rotd_p(2, 3))
        };

        const auto& p = abs_points;
        const auto& px = rel_points;

        const gp_Pnt centroid(
                (p[0].X() + p[1].X() + p[2].X() + p[3].X())/4.0,
                (p[0].Y() + p[1].Y() + p[2].Y() + p[3].Y())/4.0,
                (p[0].Z() + p[1].Z() + p[2].Z() + p[3].Z())/4.0);

        for(size_t ri = 0; ri < fn.size(); ++ri){
            const auto& gsi = utils::GaussLegendre<2>::get();
            for(auto xi = gsi.begin(); xi < gsi.end(); ++xi){
                for(auto eta = gsi.begin(); eta < gsi.end(); ++eta){
                    for(size_t i = 0; i < fn.size(); ++i){
                        gp_Pnt pxi(this->GL_point(*xi, *eta, px));
                        const double detJ = std::abs(this->J(xi->x, eta->x, px).determinant());
                        dEI[ri] += xi->w*eta->w*detJ*fn[i](S, dS, pxi);
                    }
                }
            }
        }
    }


    math::Matrix dEIm({dEI[0], dEI[2],
                       dEI[2], dEI[1]}, 2, 2);

    result -= this->ABz_bar.T()*dEIm*this->ABz;

    return result;
}

void Curvature::calculate_stress_field_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh){
    constexpr size_t MAX_B = 6;

    const size_t u_offset = this->u.size() - MAX_B;

    solver.initialize_matrix(false, this->u.size());

    const size_t K_DIM = this->elem_info->get_k_dimension();
    const size_t NODE_DOF = this->elem_info->get_dof_per_node();
    const size_t NODES_PER_ELEM = this->elem_info->get_nodes_per_element();

    std::vector<long> posi(K_DIM);
    std::vector<long> posj(K_DIM + 6);
    std::vector<long> pos_last(6);

    for(size_t i = 0; i < 6; ++i){
        posj[K_DIM + i] = u_offset + i;
        pos_last[i] = u_offset + i;
    }

    for(const auto& e:boundary_mesh){
        const gp_Pnt centroid = e->get_centroid();
        const gp_Pnt c = utils::change_point(centroid, this->rot3D);
        const auto D = this->get_D_3D(e->parent, c);
        const auto S = this->get_S_3D(e->parent, c);
        const auto TDT = this->get_TDT_3D(S, D);
        const auto DT = this->get_DT_3D(S, D);
        math::slicer::from_node_id(e->nodes, NODES_PER_ELEM, NODE_DOF, posi);
        std::copy(posi.begin(), posi.end(), posj.begin());

        //const auto K = e->get_K_ext(TDT, this->center);
        const auto St = e->get_stress_integrals(DT, this->center);
        const auto tau_xz_dz = e->get_equilibrium_partial(TDT, this->center, {0, 3});
        const auto tau_yz_dz = e->get_equilibrium_partial(TDT, this->center, {1, 3});
        const auto sigma_z_dz = e->get_equilibrium_partial(TDT, this->center, {4, 5});
        //solver.add_element(K, posi, posj);
        solver.add_element(St, pos_last, posj);
        solver.add_element(tau_xz_dz + tau_yz_dz + sigma_z_dz, posi, posj);
    }

    solver.compute();
    if(V_v != 0 || V_u != 0){

        math::Vector vv{V_u, -V_v};
        math::Matrix EI({EI_vv, EI_uv,
                         EI_uv, EI_uu}, 2, 2);

        math::LU lu(EI);
        lu.solve(vv);

        const double Az = vv[0];
        const double Bz = vv[1];
        if(this->calculate_adjoint){
            this->ABz = vv;
            this->ABz_bar = math::Vector(vv.get_N(), 0);
        }

        std::vector<size_t> posu(NODES_PER_ELEM*NODE_DOF);
        for(const auto& e:boundary_mesh){
            const gp_Pnt centroid = e->get_centroid();
            const gp_Pnt c = utils::change_point(centroid, this->rot3D);
            const auto S = this->get_S_3D(e->parent, c);
            const auto D = this->get_D_3D(e->parent, c);

            const auto dz = e->get_dz_vector(S, D, Az, Bz, this->center);

            math::slicer::from_node_id(e->nodes, NODES_PER_ELEM, NODE_DOF, posu);
            math::VectorSlice us(this->u, posu);
            us += dz;
        }
    }

    // A B C theta Kyz Kxz
    this->u[u_offset + 0] = -M_v;
    this->u[u_offset + 1] = -M_u;

    // C
    this->u[u_offset + 2] = -V_w;

    // theta
    this->u[u_offset + 3] = -M_w;

    // t_yz
    this->u[u_offset + 4] = V_v;

    // t_xz
    this->u[u_offset + 5] = -V_u;

    solver.solve(this->u);

    this->A     = this->u[u_offset + 0];
    this->B     = this->u[u_offset + 1];
    this->C     = this->u[u_offset + 2];
    this->theta = this->u[u_offset + 3];
    this->Kyz   = this->u[u_offset + 4];
    this->Kxz   = this->u[u_offset + 5];

    logger::quick_log("A", A);
    logger::quick_log("B", B);
    logger::quick_log("C", C);
    logger::quick_log("theta", theta);
    logger::quick_log("Kyz", Kyz);
    logger::quick_log("Kxz", Kxz);
    logger::quick_log("");

    if(!this->calculate_adjoint){
        this->solver.clear_matrix();
    } else {
        this->u_bar.resize(this->u.size(), 0);
    }
}

void Curvature::solve_adjoint_problem(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const std::vector<double>& fea_u){
    logger::log_assert(this->calculate_adjoint, logger::ERROR, "set_calculate_adjoint() must be called before generate_curvature_3D() to enable adjoint problem solution");

    const size_t K_DIM = this->elem_info->get_k_dimension();
    const size_t NODE_DOF = this->elem_info->get_dof_per_node();
    const size_t NODES_PER_ELEM = this->elem_info->get_nodes_per_element();
    const size_t u_offset = this->u.size() - 6;

    std::vector<size_t> posi(K_DIM + 6);
    for(size_t i = 0; i < 6; ++i){
        posi[K_DIM + i] = u_offset + i;
    }
    for(const auto& e:boundary_mesh){
        gp_Pnt c = e->get_centroid();
        gp_Pnt c2 = utils::change_point(c, this->rot3D);

        const auto D = this->get_D_3D(e->parent, c2);
        const auto S = this->get_S_3D(e->parent, c2);
        const auto DT = this->get_DT_3D(S, D);

        math::slicer::from_node_id(e->nodes, NODES_PER_ELEM, NODE_DOF, posi);
        math::VectorSlice bar_us(this->u_bar, posi);
        const auto v(e->get_force_vector(DT, fea_u, this->center, this->rot3D, true));
        bar_us += v;
    }

    this->solver.solve_adjoint(this->u_bar);
    this->solver.clear_matrix();

    if(V_v != 0 || V_u != 0){

        math::Matrix EI({EI_vv, EI_uv,
                         EI_uv, EI_uu}, 2, 2);

        math::LU lu(EI);

        std::vector<size_t> posu(NODES_PER_ELEM*NODE_DOF);
        for(const auto& e:boundary_mesh){
            const gp_Pnt centroid = e->get_centroid();
            const gp_Pnt c = utils::change_point(centroid, this->rot3D);
            const auto S = this->get_S_3D(e->parent, c);
            const auto D = this->get_D_3D(e->parent, c);

            const auto dzm = e->get_dz_vector_matrix(S, D, this->center);

            math::slicer::from_node_id(e->nodes, NODES_PER_ELEM, NODE_DOF, posu);
            math::VectorSliceView us(this->u_bar, posu);
            this->ABz_bar += dzm.T()*us;
        }

        lu.solve(ABz_bar);
    }

}

class VectorWrapper{
    public:
    math::Vector vec;
    VectorWrapper() = default;
    VectorWrapper(size_t s):vec(s){}

    VectorWrapper& operator+=(const VectorWrapper& v){
        if(this->vec.get_N() == 0){
            this->vec = math::Vector(v.vec.get_N());
        }
        this->vec += v.vec;

        return *this;
    }
};

math::Vector Curvature::integrate_surface_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const std::vector<std::function<double(const math::Matrix& S, const gp_Pnt& px)>>& fn) const{
    VectorWrapper result(fn.size());
    #pragma omp declare reduction(vecsum : VectorWrapper : omp_out += omp_in) 
    if(this->elem_info->get_shape_type() == Element::Shape::TRI){
        constexpr size_t N = 3;
        #pragma omp parallel
        {
            VectorWrapper result_tmp(fn.size());
            #pragma omp for reduction(vecsum:result)
            for(size_t i = 0; i < boundary_mesh.size(); ++i){
                math::Matrix points(3, N);
                const auto& e = boundary_mesh[i];
                std::array<gp_Pnt, N> rel_points;
                for(size_t x = 0; x < N; ++x){
                    rel_points[x] = e->nodes[x]->point;
                    for(size_t y = 0; y < 3; ++y){
                        points(y, x) = rel_points[x].Coord(y+1);
                    }
                }
                math::Matrix rotd_p(this->rot3D*points);
                std::array<gp_Pnt, N> abs_points{
                    gp_Pnt(rotd_p(0, 0), rotd_p(1, 0), rotd_p(2, 0)),
                    gp_Pnt(rotd_p(0, 1), rotd_p(1, 1), rotd_p(2, 1)),
                    gp_Pnt(rotd_p(0, 2), rotd_p(1, 2), rotd_p(2, 2))
                };
                this->GS_tri(e->parent, abs_points, rel_points, fn, result_tmp.vec);
                result += result_tmp;
            }
        }
    } else if(this->elem_info->get_shape_type() == Element::Shape::QUAD){
        constexpr size_t N = 4;
        #pragma omp parallel
        {
            VectorWrapper result_tmp(fn.size());
            #pragma omp for reduction(vecsum:result)
            for(size_t i = 0; i < boundary_mesh.size(); ++i){
                math::Matrix points(3, N);
                const auto& e = boundary_mesh[i];
                std::array<gp_Pnt, N> rel_points;
                for(size_t x = 0; x < N; ++x){
                    rel_points[x] = e->nodes[x]->point;
                    for(size_t y = 0; y < 3; ++y){
                        points(y, x) = rel_points[x].Coord(y+1);
                    }
                }
                math::Matrix rotd_p(this->rot3D*points);
                std::array<gp_Pnt, N> abs_points{
                    gp_Pnt(rotd_p(0, 0), rotd_p(1, 0), rotd_p(2, 0)),
                    gp_Pnt(rotd_p(0, 1), rotd_p(1, 1), rotd_p(2, 1)),
                    gp_Pnt(rotd_p(0, 2), rotd_p(1, 2), rotd_p(2, 2)),
                    gp_Pnt(rotd_p(0, 3), rotd_p(1, 3), rotd_p(2, 3))
                };
                this->GS_quad(e->parent, abs_points, rel_points, fn, result_tmp.vec);
                result += result_tmp;
            }
        }
    }
    return result.vec;
}

math::Matrix Curvature::get_S_3D(const MeshElement* const e, const gp_Pnt& p) const{
    auto S = this->mat->stiffness_inverse_3D(e, p);
    //this->mat->rotate_S(S, this->permute_and_rotate_3D_inv);

    return S;
}
math::Matrix Curvature::get_D_3D(const MeshElement* const e, const gp_Pnt& p) const{
    auto D = this->mat->stiffness_3D(e, p);
    //this->mat->rotate_D(D, this->permute_and_rotate_3D);

    return D;
}
    

math::Matrix Curvature::get_T_D_3D(const math::Matrix& S, const math::Matrix& D) const{
    return math::Matrix(
        {
            1, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0,
            0, 0, 1.0/(S(2,2)*D(2,2)), 0, 0, 0,
            0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 1
        }, 6, 6);
}

math::Matrix Curvature::get_T_D_3D_1d(const math::Matrix& S, const math::Matrix& D, const math::Matrix& dS, const math::Matrix& dD) const{
    const double SD = S(2,2)*D(2,2);
    const double dSD = (S(2,2)*dD(2,2) + dS(2,2)*D(2,2));
    return math::Matrix(
        {
            1, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0,
            0, 0, -dSD/(SD*SD), 0, 0, 0,
            0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 1
        }, 6, 6);
}

math::Matrix Curvature::get_DT_3D(const math::Matrix& S, const math::Matrix& D) const{
    const auto T = this->get_T_D_3D(S, D);

    return D*T;
}
math::Matrix Curvature::get_DT_3D_1d(const math::Matrix& S, const math::Matrix& D, const math::Matrix& dS, const math::Matrix& dD) const{
    const auto T = this->get_T_D_3D(S, D);
    const auto dT = this->get_T_D_3D_1d(S, D, dS, dD);

    return D*dT + dD*T;
}
math::Matrix Curvature::get_TDT_3D(const math::Matrix& S, const math::Matrix& D) const{
    const auto T = this->get_T_D_3D(S, D);

    return T.T()*D*T;
}
math::Matrix Curvature::get_TDT_3D_1d(const math::Matrix& S, const math::Matrix& D, const math::Matrix& dS, const math::Matrix& dD) const{
    const auto T = this->get_T_D_3D(S, D);
    const auto dT = this->get_T_D_3D_1d(S, D, dS, dD);

    return T.T()*(dD*T + D*dT) + dT.T()*D*T;
}

void Curvature::GS_tri(const MeshElement* const e, const std::array<gp_Pnt, 3>& p, const std::array<gp_Pnt, 3>& px, const std::vector<std::function<double(const math::Matrix& S, const gp_Pnt& px)>>& fn, math::Vector& result) const{
    gp_Vec v1(p[1], p[0]);
    gp_Vec v2(p[2], p[0]);
    result.fill(0);
    const gp_Pnt centroid(
            (p[0].X() + p[1].X() + p[2].X())/3.0,
            (p[0].Y() + p[1].Y() + p[2].Y())/3.0,
            (p[0].Z() + p[1].Z() + p[2].Z())/3.0);
    const auto S = this->get_S_3D(e, centroid);
    const double drnorm = (v1.Crossed(v2)).Magnitude()/2;
    const auto& gsi = utils::GaussLegendreTri<2>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        for(size_t i = 0; i < fn.size(); ++i){
            gp_Pnt pxi(this->GLT_point(*it, px));
            result[i] += it->w*fn[i](S, pxi);
        }
    }
    result *= drnorm;
}

void Curvature::GS_quad(const MeshElement* const e, const std::array<gp_Pnt, 4>& p, const std::array<gp_Pnt, 4>& px, const std::vector<std::function<double(const math::Matrix& S, const gp_Pnt& px)>>& fn, math::Vector& result) const{
    result.fill(0);
    const gp_Pnt centroid(
            (p[0].X() + p[1].X() + p[2].X() + p[3].X())/4.0,
            (p[0].Y() + p[1].Y() + p[2].Y() + p[3].Y())/4.0,
            (p[0].Z() + p[1].Z() + p[2].Z() + p[3].Z())/4.0);
    const auto S = this->get_S_3D(e, centroid);
    const auto& gsi = utils::GaussLegendre<2>::get();
    for(auto xi = gsi.begin(); xi < gsi.end(); ++xi){
        for(auto eta = gsi.begin(); eta < gsi.end(); ++eta){
            for(size_t i = 0; i < fn.size(); ++i){
                gp_Pnt pxi(this->GL_point(*xi, *eta, px));
                const double detJ = std::abs(this->J(xi->x, eta->x, px).determinant());
                result[i] += xi->w*eta->w*detJ*fn[i](S, pxi);
            }
        }
    }
}
