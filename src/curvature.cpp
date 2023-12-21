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
#include <Eigen/src/Core/Matrix.h>
#include <type_traits>
#include <vector>
#include "general_solver/mumps_general.hpp"
#include "logger.hpp"
#include "utils/gauss_legendre.hpp"
#include "utils/boundary_nullifier.hpp"
#include "utils/basis_tensor.hpp"

Curvature::Curvature(const Material* mat, gp_Dir u, gp_Dir v, gp_Dir w, Eigen::Matrix<double, 2, 2> rot2D, Eigen::Matrix<double, 3, 3> rot3D, const BoundaryMeshElementFactory* elem_info, double V_u, double V_v, double V_w, double M_u, double M_v, double M_w):
    mat(mat), 
    rot2D(rot2D), 
    Lek_basis
        {{0, 0, -1},
         {0, 1, 0},
         {1, 0, 0}},
    rot3D(rot3D), 
    elem_info(elem_info),
    u(rot3D(0,0), rot3D(1,0), rot3D(2,0)), 
    v(rot3D(0,1), rot3D(1,1), rot3D(2,1)), 
    w(rot3D(0,2), rot3D(1,2), rot3D(2,2)), 
    V_u(-V_w), V_v(V_v), V_w(V_u),
    M_u(-M_w), 
    M_v(M_v),
    M_w(M_u)
{

    this->permute_shear_3D = std::vector<double>
        {1, 0, 0, 0, 0, 0,
         0, 1, 0, 0, 0, 0,
         0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 0, 1,
         0, 0, 0, 0, 1, 0,
         0, 0, 0, 1, 0, 0};

    this->rotate_X_Z_3D = utils::basis_tensor_3D(Lek_basis);
    this->rotate_X_Z_3D_inv = utils::basis_tensor_3D(Lek_basis.transpose());

    this->permute_and_rotate_3D = std::vector<double>(36, 0);
    this->permute_and_rotate_3D_inv = std::vector<double>(36, 0);
    const size_t N = 6;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1, permute_shear_3D.data(), N, rotate_X_Z_3D.data(), N, 0, permute_and_rotate_3D.data(), N);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1, permute_shear_3D.data(), N, rotate_X_Z_3D_inv.data(), N, 0, permute_and_rotate_3D_inv.data(), N);
}

void Curvature::generate_curvature_3D(const std::vector<std::unique_ptr<MeshNode>>& boundary_nodes, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, size_t reduced_vector_size, size_t full_vector_size){
    std::array<gp_Pnt, 3> points;

    const auto fn_A =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_A_base_3D(S, px);
        };
    const auto fn_EA =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_EA_base_3D(S, px);
        };
    const auto fn_EA_u =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_EA_u_base_3D(S, px);
        };
    const auto fn_EA_v =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_EA_v_base_3D(S, px);
        };
    const auto fn_EA_w =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_EA_w_base_3D(S, px);
        };
    const auto fn_EI_vv =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_EI_vv_base_3D(S, px);
        };
    const auto fn_EI_uu =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_EI_uu_base_3D(S, px);
        };
    const auto fn_EI_uv =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_EI_uv_base_3D(S, px);
        };
    const auto fn_EI_vvv =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_EI_vvv_base_3D(S, px);
        };
    const auto fn_EI_uuu =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_EI_uuu_base_3D(S, px);
        };
    const auto fn_EI_uuv =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_EI_uuv_base_3D(S, px);
        };
    const auto fn_EI_uvv =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_EI_uvv_base_3D(S, px);
        };
    const auto fn_S34_vv =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_S34_vv_base_3D(S, px);
        };
    const auto fn_S34_uu =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_S34_uu_base_3D(S, px);
        };
    const auto fn_S34_vvv =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_S34_vvv_base_3D(S, px);
        };
    const auto fn_S34_uuu =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_S34_uuu_base_3D(S, px);
        };
    const auto fn_S34_uuv =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_S34_uuv_base_3D(S, px);
        };
    const auto fn_S34_uvv =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_S34_uvv_base_3D(S, px);
        };
    const auto fn_S35_vv =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_S35_vv_base_3D(S, px);
        };
    const auto fn_S35_uu =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_S35_uu_base_3D(S, px);
        };
    const auto fn_S35_vvv =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_S35_vvv_base_3D(S, px);
        };
    const auto fn_S35_uuu =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_S35_uuu_base_3D(S, px);
        };
    const auto fn_S35_uuv =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_S35_uuv_base_3D(S, px);
        };
    const auto fn_S35_uvv =
        [this](const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)->double{
            return this->make_S35_uvv_base_3D(S, px);
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
    const auto result2 =
        this->integrate_surface_3D(boundary_mesh, 
                {fn_EI_uu,
                fn_EI_vv,
                fn_EI_uv,
                fn_EI_uuu,
                fn_EI_vvv,
                fn_EI_uuv,
                fn_EI_uvv,
                fn_S34_uu,
                fn_S34_vv,
                fn_S34_vvv,
                fn_S34_uuu,
                fn_S34_vvv,
                fn_S34_uuv,
                fn_S34_uvv,
                fn_S35_uu,
                fn_S35_vv,
                fn_S35_uuu,
                fn_S35_vvv,
                fn_S35_uuv,
                fn_S35_uvv});
    this->EI_uu = result2[0];
    this->EI_vv = result2[1];
    this->EI_uv = result2[2];
    this->EI_uuu = result2[3];
    this->EI_vvv = result2[4];
    this->EI_uuv = result2[5];
    this->EI_uvv = result2[6];
    this->S34_uu = result2[7];
    this->S34_vv = result2[8];
    this->S34_vvv = result2[9];
    this->S34_uuu = result2[10];
    this->S34_vvv = result2[11];
    this->S34_uuv = result2[12];
    this->S34_uvv = result2[13];
    this->S35_uu = result2[14];
    this->S35_vv = result2[15];
    this->S35_uuu = result2[16];
    this->S35_vvv = result2[17];
    this->S35_uuv = result2[18];
    this->S35_uvv = result2[19];
    logger::quick_log("A: ", Area);
    logger::quick_log("EA: ", EA);
    logger::quick_log("c_v: ", c_v);
    logger::quick_log("c_u: ", c_u);
    logger::quick_log("EA_u", EA_u);
    logger::quick_log("EA_v", EA_v);
    logger::quick_log("EI_uu", EI_uu);
    logger::quick_log("EI_vv", EI_vv);
    logger::quick_log("EI_uv", EI_uv);
    logger::quick_log("EI_uuu", EI_uuu);
    logger::quick_log("EI_vvv", EI_vvv);
    logger::quick_log("EI_uuv", EI_uuv);
    logger::quick_log("EI_uvv", EI_uvv);
    logger::quick_log("S34_uu", S34_uu);
    logger::quick_log("S34_vv", S34_vv);
    logger::quick_log("S34_uuu", S34_uuu);
    logger::quick_log("S34_vvv", S34_vvv);
    logger::quick_log("S34_uuv", S34_uuv);
    logger::quick_log("S34_uvv", S34_uvv);
    logger::quick_log("S35_uu", S35_uu);
    logger::quick_log("S35_vv", S35_vv);
    logger::quick_log("S35_uuu", S35_uuu);
    logger::quick_log("S35_vvv", S35_vvv);
    logger::quick_log("S35_uuv", S35_uuv);
    logger::quick_log("S35_uvv", S35_uvv);

    this->reduced_vec_len = reduced_vector_size;
    this->F.resize(2*full_vector_size,0);
    this->phi.resize(full_vector_size,0);

    this->calculate_stress_field_3D(boundary_nodes, boundary_mesh);
}

void Curvature::get_stress_3D(const BoundaryMeshElement* e, double& s_w, double& t_uw, double& t_vw) const{
    s_w = 0;
    gp_Pnt c = e->get_centroid();
    gp_Pnt c2 = utils::change_point(c, this->rot3D);
    const double x = c.X() - this->c_u;
    const double y = c.Y() - this->c_v;

    t_uw = 0;
    t_vw = 0;
    Eigen::Vector<double, 2> grad_phi = e->grad_1dof_id(c, this->phi);
    Eigen::Vector<double, 3> grad_F = e->dF_2dof_id(c, this->F);
    t_uw =  grad_phi[1];
    t_vw = -grad_phi[0];

    const auto S = this->get_S_3D(e->parent, c2);
    double s_u = grad_F[1];
    double s_v = grad_F[0];
    double t_uv = -grad_F[2];
    t_uw += (Kxz_1*x*x + Kxz_2*y*y)/S(2,2);
    t_vw += (Kyz_1*y*y + Kyz_2*x*x)/S(2,2);

    s_w = (A*x + B*y + C)/S(2,2) - (S(0,2)*s_u + S(1,2)*s_v + S(2,3)*t_vw + S(2,4)*t_uw + S(2,5)*t_uv)/S(2,2);
}

void Curvature::calculate_stress_field_3D(const std::vector<std::unique_ptr<MeshNode>>& boundary_nodes, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh){
    const size_t num_nodes = this->elem_info->get_nodes_per_element();

    const size_t F_offset = 2*this->reduced_vec_len;
    size_t phi_size = F_offset + this->reduced_vec_len;
    const gp_Pnt center(c_u, c_v, 0);

    // Eigen solvers DO NOT WORK CORRECTLY HERE
    general_solver::MUMPSGeneral solver;
    solver.initialize_matrix(false, phi_size);

    /* 
     * First part: obtain phi_1 and F_1.
     *
     * Both come from: 
     * phi = theta*phi_1 + phi_0
     * F = theta*F_1 + F_0
     *
     * And need to be obtained first so that theta, A and B may be calculated.
     * Based on Lekhnitskii formalism.
     *
     */
    std::vector<double> b(phi_size, 0);

    this->base_matrix_upos(solver, boundary_mesh, num_nodes, F_offset);

    for(const auto& e:boundary_mesh){
        const Eigen::VectorXd N = 2*e->source_1dof();

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->u_pos[0];
            if(id1 < 0){
                continue;
            }
            b[F_offset + id1] += N[i];
        }
    }


    solver.compute();

    solver.solve(b);
    auto phi_tmp = b;

    /*
     * Save phi_1 and F_1 in the global vector.
     *
     */
    for(const auto& n:boundary_nodes){
        const long id1 = n->u_pos[0];
        if(id1 < 0){
            continue;
        }
        this->F[2*n->id] += phi_tmp[2*id1];
        this->F[2*n->id+1] += phi_tmp[2*id1+1];
    }
    for(const auto& n:boundary_nodes){
        const long id1 = n->u_pos[0];
        if(id1 < 0){
            continue;
        }
        this->phi[n->id] += phi_tmp[F_offset + id1];
    }

    /* 
     * Second part: obtain phi_0 and F_0.
     *
     * Matrix is slightly larger to acomodate theta, A and B.
     * Unfortunately, such addition makes it necessary to redo M and refactorize
     * it.
     *
     * It only isn't necessary when S is at most orthotropic, but it is not
     * really necessary to optimize for it as the solver would already ignore
     * the zero-filled rectangles.
     */
    const size_t F_phi_offset = F_offset + this->phi.size();;
    // A B C theta Kyz_2 Kxz_2
    phi_size = F_phi_offset + 4 + 2;
    solver.initialize_matrix(false, phi_size);

    this->base_matrix_id(solver, boundary_mesh, num_nodes, F_offset);

    b.resize(phi_size);
    std::fill(b.begin(), b.end(), 0);

    // Az Bz
    Eigen::Matrix<double, 2, 2> EIM{{EI_vv, EI_uv},
                                    {EI_uv, EI_uu}};
    Eigen::Vector<double, 2> VV{-V_u, V_v};
    Eigen::Vector<double, 2> Vz = EIM.fullPivLu().solve(VV);

    const double Az = Vz[0];
    const double Bz = Vz[1];

    logger::quick_log("Az", Az);
    logger::quick_log("Bz", Bz);
    this->Kxz_1 = -Az/2.0;
    this->Kyz_1 = -Bz/2.0;

    /*
     * Generate the constants which multiply theta and fill the additional
     * rectangles in the matrix.
     */
    double C1 = 0, Cx1 = 0, Cy1 = 0, Ct1 = 0;
    double xC1 = 0, yC1 = 0;
    Eigen::Matrix<double, 2, 2> I{{1,0},
                                  {0,1}};
    for(const auto& e:boundary_mesh){
        const auto N = e->source_1dof();
        const auto Nx = e->int_N_x(center);
        const auto Ny = e->int_N_y(center);
        const auto dphi = e->int_grad_phi();
        const auto dphix = e->int_grad_phi_x(center);
        const auto dphiy = e->int_grad_phi_y(center);
        const auto dF = e->int_grad_F();
        const auto dFx = e->int_grad_F_x(center);
        const auto dFy = e->int_grad_F_y(center);
        gp_Pnt c = utils::change_point(e->get_centroid(), this->rot3D);
        const auto S = this->get_S_3D(e->parent, c);
        const auto B = this->get_B_3D(e->parent, c);

        for(size_t i = 0; i < num_nodes; ++i){
            // phi vector
            const long id2 = e->nodes[i]->id;

            Cx1 += (-S(2,3)*dphix(0,i)+S(2,4)*dphix(1,i))*phi[id2]/S(2,2);
            Cy1 += (-S(2,3)*dphiy(0,i)+S(2,4)*dphiy(1,i))*phi[id2]/S(2,2);
            C1  += (-S(2,3)*dphi (0,i)+S(2,4)*dphi (1,i))*phi[id2]/S(2,2);
            // In this case, phi = 0 at the boundary
            Ct1 += 2*N[i]*phi[id2];

            xC1 += (-dphi(0,i))*phi[id2];
            yC1 += ( dphi(1,i))*phi[id2];

            // Cx0
            solver.add_value(F_phi_offset + 0, F_offset + id2, -(-S(2,3)*dphix(0,i)+S(2,4)*dphix(1,i))/S(2,2));
            // Cy0
            solver.add_value(F_phi_offset + 1, F_offset + id2, -(-S(2,3)*dphiy(0,i)+S(2,4)*dphiy(1,i))/S(2,2));
            // C0                                                    
            solver.add_value(F_phi_offset + 2, F_offset + id2, -(-S(2,3)*dphi (0,i)+S(2,4)*dphi (1,i))/S(2,2));
            // Ct0
            // In this case, phi != 0 at the boundary
            solver.add_value(F_phi_offset + 3, F_offset + id2, -dphix(0,i)-dphiy(1,i));

            // Kyz_2, Kxz_2
            // xC0
            solver.add_value(F_phi_offset + 4, F_offset + id2, (-dphi(0,i)));
            // yC0
            solver.add_value(F_phi_offset + 5, F_offset + id2, ( dphi(1,i)));

            // phi source
            solver.add_value(F_offset + id2, F_phi_offset + 4 + 0,  2*B(3,3)*Nx[i]/S(2,2));
            solver.add_value(F_offset + id2, F_phi_offset + 4 + 1, -2*B(4,4)*Ny[i]/S(2,2));
            b[F_offset + id2] += 2*B(3,4)*Ny[i]*Kyz_1/S(2,2) - 2*B(3,4)*Nx[i]*Kxz_1/S(2,2);
            // A coeffs
            solver.add_value(F_offset + id2, F_phi_offset + 0, -(-S(2,3)*N[i]/S(2,2)));
            // B coeffs
            solver.add_value(F_offset + id2, F_phi_offset + 1, -(S(2,4)*N[i]/S(2,2)));

            // F vector
            const long id1 = e->nodes[i]->u_pos[0];
            if(id1 < 0){
                continue;
            }
            Cx1 +=  (
                     (S(1,2)*dFx(0,2*i+0)-0.5*S(2,5)*dFx(2,2*i+0))*phi_tmp[2*id1+0]+
                     (S(0,2)*dFx(1,2*i+1)-0.5*S(2,5)*dFx(2,2*i+1))*phi_tmp[2*id1+1])
                     /S(2,2);
            Cy1 +=  (
                     (S(1,2)*dFy(0,2*i+0)-0.5*S(2,5)*dFy(2,2*i+0))*phi_tmp[2*id1+0]+
                     (S(0,2)*dFy(1,2*i+1)-0.5*S(2,5)*dFy(2,2*i+1))*phi_tmp[2*id1+1])
                     /S(2,2);
            C1  +=  (
                    (S(1,2)*dF (0,2*i+0)-0.5*S(2,5)*dF (2,2*i+0))*phi_tmp[2*id1+0]+
                    (S(0,2)*dF (1,2*i+1)-0.5*S(2,5)*dF (2,2*i+1))*phi_tmp[2*id1+1])
                    /S(2,2);

            // Cx0
            solver.add_value(F_phi_offset + 0, 2*id1 + 0, -(S(1,2)*dFx(0,2*i+0)-0.5*S(2,5)*dFx(2,2*i+0))/S(2,2));
            solver.add_value(F_phi_offset + 0, 2*id1 + 1, -(S(0,2)*dFx(1,2*i+1)-0.5*S(2,5)*dFx(2,2*i+1))/S(2,2));
            // Cy0
            solver.add_value(F_phi_offset + 1, 2*id1 + 0, -(S(1,2)*dFy(0,2*i+0)-0.5*S(2,5)*dFy(2,2*i+0))/S(2,2));
            solver.add_value(F_phi_offset + 1, 2*id1 + 1, -(S(0,2)*dFy(1,2*i+1)-0.5*S(2,5)*dFy(2,2*i+1))/S(2,2));
            // C0                                                    
            solver.add_value(F_phi_offset + 2, 2*id1 + 0, -(S(1,2)*dF (0,2*i+0)-0.5*S(2,5)*dF (2,2*i+0))/S(2,2));
            solver.add_value(F_phi_offset + 2, 2*id1 + 1, -(S(0,2)*dF (1,2*i+1)-0.5*S(2,5)*dF (2,2*i+1))/S(2,2));

            // F source
            solver.add_value(2*id1 + 0, F_phi_offset + 4 + 0, -2*B(0,4)*N[i]/S(2,2));
            solver.add_value(2*id1 + 0, F_phi_offset + 4 + 1, -2*B(1,3)*N[i]/S(2,2));
            solver.add_value(2*id1 + 1, F_phi_offset + 4 + 0, -2*B(0,4)*N[i]/S(2,2));
            solver.add_value(2*id1 + 1, F_phi_offset + 4 + 1, -2*B(1,3)*N[i]/S(2,2));
            b[2*id1 + 0] += 2*B(0,3)*N[i]*Kyz_1/S(2,2) + 2*B(1,4)*N[i]*Kxz_1/S(2,2);
            b[2*id1 + 1] += 2*B(0,3)*N[i]*Kyz_1/S(2,2) + 2*B(1,4)*N[i]*Kxz_1/S(2,2);
        }
    }

    // A B square
    solver.add_value(F_phi_offset + 0, F_phi_offset + 0, EI_vv);
    solver.add_value(F_phi_offset + 0, F_phi_offset + 1, EI_uv);
    solver.add_value(F_phi_offset + 1, F_phi_offset + 0, EI_uv);
    solver.add_value(F_phi_offset + 1, F_phi_offset + 1, EI_uu);

    solver.add_value(F_phi_offset + 0, F_phi_offset + 4 + 0, -S34_uvv);
    solver.add_value(F_phi_offset + 0, F_phi_offset + 4 + 1, -S35_uuu);
    b[F_phi_offset + 0] = -M_v + Kyz_1*S34_uuu + Kxz_1*S35_uvv;
    solver.add_value(F_phi_offset + 1, F_phi_offset + 4 + 0, -S34_vvv);
    solver.add_value(F_phi_offset + 1, F_phi_offset + 4 + 1, -S35_uuv);
    b[F_phi_offset + 1] = -M_u + Kyz_1*S34_uuv + Kxz_1*S35_vvv;

    // C
    solver.add_value(F_phi_offset + 2, F_phi_offset + 2, EA);
    solver.add_value(F_phi_offset + 2, F_phi_offset + 4 + 0, -S34_vv);
    solver.add_value(F_phi_offset + 2, F_phi_offset + 4 + 1, -S35_uu);
    b[F_phi_offset + 2] = -V_w + Kyz_1*S34_uu + Kxz_1*S35_vv;

    // theta
    solver.add_value(F_phi_offset + 0, F_phi_offset + 3, -Cx1);
    solver.add_value(F_phi_offset + 1, F_phi_offset + 3, -Cy1);
    solver.add_value(F_phi_offset + 2, F_phi_offset + 3, -C1 );
    solver.add_value(F_phi_offset + 3, F_phi_offset + 3,  Ct1);
    solver.add_value(F_phi_offset + 3, F_phi_offset + 4 + 0,  EI_vvv);
    solver.add_value(F_phi_offset + 3, F_phi_offset + 4 + 1, -EI_uuu);
    b[F_phi_offset + 3] = -M_w - Kyz_1*EI_uuv + Kxz_1*EI_uvv;

    // Kyz_2, Kxz_2
    solver.add_value(F_phi_offset + 4, F_phi_offset + 3,  xC1);
    solver.add_value(F_phi_offset + 5, F_phi_offset + 3,  yC1);

    // t_yz
    solver.add_value(F_phi_offset + 4, F_phi_offset + 4 + 0,  EI_vv);
    b[F_phi_offset + 4] = V_v - Kyz_1*EI_uu;
    // t_xz
    solver.add_value(F_phi_offset + 5, F_phi_offset + 4 + 1,  EI_uu);
    b[F_phi_offset + 5] = -V_u - Kxz_1*EI_vv;

    solver.compute();

    solver.solve(b);
    phi_tmp = b;

    double C0 = 0, Cx0 = 0, Cy0 = 0, Ct0 = 0;
    double xC0 = 0, yC0 = 0;
    for(const auto& e:boundary_mesh){
        const auto N = e->source_1dof();
        const auto Nx = e->int_N_x(center);
        const auto Ny = e->int_N_y(center);
        const auto M_e = e->diffusion_1dof(I);
        const auto dphi = e->int_grad_phi();
        const auto dphix = e->int_grad_phi_x(center);
        const auto dphiy = e->int_grad_phi_y(center);
        const auto dF = e->int_grad_F();
        const auto dFx = e->int_grad_F_x(center);
        const auto dFy = e->int_grad_F_y(center);
        gp_Pnt c = utils::change_point(e->get_centroid(), this->rot3D);
        const auto S = this->get_S_3D(e->parent, c);

        for(size_t i = 0; i < num_nodes; ++i){
            const long id2 = e->nodes[i]->id;

            Cx0 += (-S(2,3)*dphix(0,i)+S(2,4)*dphix(1,i))*phi_tmp[F_offset + id2]/S(2,2);
            Cy0 += (-S(2,3)*dphiy(0,i)+S(2,4)*dphiy(1,i))*phi_tmp[F_offset + id2]/S(2,2);
            C0  += (-S(2,3)*dphi (0,i)+S(2,4)*dphi (1,i))*phi_tmp[F_offset + id2]/S(2,2);
            // phi != 0 at the boundary
            Ct0 += (-dphix(0,i)-dphiy(1,i))*phi_tmp[F_offset + id2];

            xC0 += (-dphi(0,i))*phi_tmp[F_offset + id2];
            yC0 += ( dphi(1,i))*phi_tmp[F_offset + id2];
        }
    }

    this->A     = phi_tmp[F_phi_offset + 0];
    this->B     = phi_tmp[F_phi_offset + 1];
    this->C     = phi_tmp[F_phi_offset + 2];
    this->theta = phi_tmp[F_phi_offset + 3];
    this->Kyz_2 = phi_tmp[F_phi_offset + 4];
    this->Kxz_2 = phi_tmp[F_phi_offset + 5];

    /*
     * Multiply phi_1 and F_1 by theta
     */
    for(auto& f : this->F){
        f *= theta;
    }
    for(auto& f : this->phi){
        f *= theta;
    }

    /*
     * Save phi_0 and F_0 in the global vector.
     *
     */
    for(const auto& n:boundary_nodes){
        const long id1 = n->u_pos[0];
        if(id1 < 0){
            continue;
        }
        this->F[2*n->id] += phi_tmp[2*id1];
        this->F[2*n->id+1] += phi_tmp[2*id1+1];
    }
    for(const auto& n:boundary_nodes){
        const long id1 = n->id;
        this->phi[n->id] += phi_tmp[F_offset + id1];
    }

    logger::quick_log("A", A);
    logger::quick_log("B", B);
    logger::quick_log("C", C);
    logger::quick_log("theta", theta);
    logger::quick_log("Kyz_1", Kyz_1);
    logger::quick_log("Kyz_2", Kyz_2);
    logger::quick_log("Kxz_1", Kxz_1);
    logger::quick_log("Kxz_2", Kxz_2);
    logger::quick_log("Cx0", Cx0);
    logger::quick_log("Cy0", Cy0);
    logger::quick_log("C0 ", C0 );
    logger::quick_log("Ct0", Ct0);
    logger::quick_log("xC0", xC0);
    logger::quick_log("yC0", yC0);
    logger::quick_log("Cx1", Cx1);
    logger::quick_log("Cy1", Cy1);
    logger::quick_log("C1 ", C1 );
    logger::quick_log("Ct1", Ct1);
    logger::quick_log("xC1", xC1);
    logger::quick_log("yC1", yC1);
    logger::quick_log("");
}
    
void Curvature::base_matrix_id(general_solver::MUMPSGeneral& M, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const size_t num_nodes, const size_t F_offset) const{
    Eigen::Matrix<double, 3, 3> B4;
    Eigen::Matrix<double, 3, 2> B3;
    Eigen::Matrix<double, 2, 2> B2;

    std::vector<double> L4v(2*num_nodes*2*num_nodes,0);
    std::vector<double> L3v(2*num_nodes*num_nodes,0);
    std::vector<double> L3vT(2*num_nodes*num_nodes,0);
    std::vector<double> L2v(num_nodes*num_nodes,0);

    std::vector<long> pos_L4(num_nodes*2,0);
    std::vector<long> pos_L2(num_nodes,0);

    for(const auto& e:boundary_mesh){
        gp_Pnt c = utils::change_point(e->get_centroid(), this->rot3D);
        this->get_B_tensors_3D(e->parent, c, B4, B3, B2);
        const Eigen::MatrixXd L4 = e->L4(B4);
        const Eigen::MatrixXd L3 = e->L3(B3);
        const Eigen::MatrixXd L2 = e->L2(B2);

        for(size_t i = 0; i < 2*num_nodes; ++i){
            for(size_t j = 0; j < 2*num_nodes; ++j){
                L4v[i*2*num_nodes + j] = L4(i,j);
            }
        }
        for(size_t i = 0; i < 2*num_nodes; ++i){
            for(size_t j = 0; j < num_nodes; ++j){
                L3v[i*num_nodes + j] = L3(i,j);
                L3vT[j*2*num_nodes + i] = L3(i,j);
            }
        }
        for(size_t i = 0; i < num_nodes; ++i){
            for(size_t j = 0; j < num_nodes; ++j){
                L2v[i*num_nodes + j] = L2(i,j);
            }
        }

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->u_pos[0];
            if(id1 < 0){
                pos_L4[2*i] = -1;
                pos_L4[2*i+1] = -1;
            } else {
                pos_L4[2*i] = 2*id1;
                pos_L4[2*i+1] = 2*id1+1;
            }
        }
        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->id;
            pos_L2[i] = F_offset + id1;
        }

        M.add_element(L4v, pos_L4);
        M.add_element(L3v, pos_L4, pos_L2);
        M.add_element(L3vT, pos_L2, pos_L4);
        M.add_element(L2v, pos_L2);
    }
}

void Curvature::base_matrix_upos(general_solver::MUMPSGeneral& M, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const size_t num_nodes, const size_t F_offset) const{
    Eigen::Matrix<double, 3, 3> B4;
    Eigen::Matrix<double, 3, 2> B3;
    Eigen::Matrix<double, 2, 2> B2;

    std::vector<double> L4v(2*num_nodes*2*num_nodes,0);
    std::vector<double> L3v(2*num_nodes*num_nodes,0);
    std::vector<double> L3vT(2*num_nodes*num_nodes,0);
    std::vector<double> L2v(num_nodes*num_nodes,0);

    std::vector<long> pos_L4(num_nodes*2,0);
    std::vector<long> pos_L2(num_nodes,0);

    for(const auto& e:boundary_mesh){
        gp_Pnt c = utils::change_point(e->get_centroid(), this->rot3D);
        this->get_B_tensors_3D(e->parent, c, B4, B3, B2);
        const Eigen::MatrixXd L4 = e->L4(B4);
        const Eigen::MatrixXd L3 = e->L3(B3);
        const Eigen::MatrixXd L2 = e->L2(B2);

        for(size_t i = 0; i < 2*num_nodes; ++i){
            for(size_t j = 0; j < 2*num_nodes; ++j){
                L4v[i*2*num_nodes + j] = L4(i,j);
            }
        }
        for(size_t i = 0; i < 2*num_nodes; ++i){
            for(size_t j = 0; j < num_nodes; ++j){
                L3v[i*num_nodes + j] = L3(i,j);
                L3vT[j*2*num_nodes + i] = L3(i,j);
            }
        }
        for(size_t i = 0; i < num_nodes; ++i){
            for(size_t j = 0; j < num_nodes; ++j){
                L2v[i*num_nodes + j] = L2(i,j);
            }
        }

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->u_pos[0];
            if(id1 < 0){
                pos_L4[2*i] = -1;
                pos_L4[2*i+1] = -1;
            } else {
                pos_L4[2*i] = 2*id1;
                pos_L4[2*i+1] = 2*id1+1;
            }
        }
        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->u_pos[0];
            if(id1 < 0) {
                pos_L2[i] = -1;
            } else {
                pos_L2[i] = F_offset + id1;
            }
        }

        M.add_element(L4v, pos_L4);
        M.add_element(L3v, pos_L4, pos_L2);
        M.add_element(L3vT, pos_L2, pos_L4);
        M.add_element(L2v, pos_L2);
    }
}

class VectorWrapper{
    public:
    Eigen::VectorXd vec;
    VectorWrapper() = default;
    VectorWrapper(size_t s):vec(s){
        vec.fill(0);
    }

    VectorWrapper& operator+=(const VectorWrapper& v){
        if(this->vec.size() == 0){
            this->vec.resize(v.vec.size());
            this->vec.fill(0);
        }
        this->vec += v.vec;

        return *this;
    }
};

Eigen::VectorXd Curvature::integrate_surface_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const std::vector<std::function<double(const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)>>& fn) const{
    VectorWrapper result(fn.size());
    #pragma omp declare reduction(vecsum : VectorWrapper : omp_out += omp_in) 
    if(this->elem_info->get_shape_type() == Element::Shape::TRI){
        #pragma omp parallel
        {
            VectorWrapper result_tmp(fn.size());
            #pragma omp for reduction(vecsum:result)
            for(size_t i = 0; i < boundary_mesh.size(); ++i){
                Eigen::Matrix<double, 3, 3> points{{0,0,0},{0,0,0},{0,0,0}};
                const auto& e = boundary_mesh[i];
                std::array<gp_Pnt, 3> rel_points;
                for(size_t x = 0; x < 3; ++x){
                    rel_points[x] = e->nodes[x]->point;
                    for(size_t y = 0; y < 3; ++y){
                        points(y, x) = rel_points[x].Coord(y+1);
                    }
                }
                Eigen::Matrix<double, 3, 3> rotd_p = this->rot3D*points;
                std::array<gp_Pnt, 3> abs_points{
                    gp_Pnt(rotd_p(0, 0), rotd_p(1, 0), rotd_p(2, 0)),
                    gp_Pnt(rotd_p(0, 1), rotd_p(1, 1), rotd_p(2, 1)),
                    gp_Pnt(rotd_p(0, 2), rotd_p(1, 2), rotd_p(2, 2))
                };
                this->GS_tri(e->parent, abs_points, rel_points, fn, result_tmp.vec);
                result += result_tmp;
            }
        }
    } else if(this->elem_info->get_shape_type() == Element::Shape::QUAD){

    }
    return result.vec;
}

Eigen::Matrix<double, 6, 6> Curvature::get_S_3D(const MeshElement* const e, const gp_Pnt& p) const{
    auto S = this->mat->stiffness_inverse_3D(e, p);
    this->mat->rotate_S_3D(S, this->permute_shear_3D);

    return Eigen::Map<Eigen::Matrix<double, 6, 6>>(S.data(), 6, 6);
}
Eigen::Matrix<double, 6, 6> Curvature::get_B_3D(const MeshElement* const e, const gp_Pnt& p) const{
    const auto S = this->get_S_3D(e, p);
    Eigen::Matrix<double, 6, 6> B;

    for(size_t i = 0; i < 6; ++i){
        for(size_t j = 0; j < 6; ++j){
            B(i,j) = S(i,j) - S(i,2)*S(j,2)/S(2,2);
        }
    }

    return B;
}
    
void Curvature::get_B_tensors_3D(const MeshElement* const e, const gp_Pnt& p, Eigen::Matrix<double, 3, 3>& B4, Eigen::Matrix<double, 3, 2>& B3, Eigen::Matrix<double, 2, 2>& B2) const{
    const auto B = this->get_B_3D(e, p);

    B4 = Eigen::Matrix<double, 3, 3>
         {{ B(1,1),  B(0,1), -B(1,5)/2},
          { B(0,1),  B(0,0), -B(0,5)/2},
          {-B(1,5), -B(0,5),  B(5,5)/2}};

    B3 = Eigen::Matrix<double, 3, 2>
         {{-B(1,3)  ,  B(1,4)  },
          {-B(0,3)  ,  B(0,4)  },
          { B(3,5)/2, -B(4,5)/2}};

    B2 = Eigen::Matrix<double, 2, 2>
         {{ B(3,3), -B(3,4)},
          {-B(3,4),  B(4,4)}};
}

void Curvature::GS_tri(const MeshElement* const e, const std::array<gp_Pnt, 3>& p, const std::array<gp_Pnt, 3>& px, const std::vector<std::function<double(const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)>>& fn, Eigen::VectorXd& result) const{
    gp_Vec v1(p[1], p[0]);
    gp_Vec v2(p[2], p[0]);
    result.fill(0);
    const double drnorm = (v1.Crossed(v2)).Magnitude()/2;
    const auto& gsi = utils::GaussLegendreTri<6>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        for(size_t i = 0; i < fn.size(); ++i){
            gp_Pnt pi{
                it->a*p[0].X() + it->b*p[1].X() + it->c*p[2].X(),
                it->a*p[0].Y() + it->b*p[1].Y() + it->c*p[2].Y(),
                it->a*p[0].Z() + it->b*p[1].Z() + it->c*p[2].Z()
            };
            gp_Pnt pxi{
                it->a*px[0].X() + it->b*px[1].X() + it->c*px[2].X(),
                it->a*px[0].Y() + it->b*px[1].Y() + it->c*px[2].Y(),
                it->a*px[0].Z() + it->b*px[1].Z() + it->c*px[2].Z()
            };
            const auto S = this->get_S_3D(e, pi);
            result[i] += it->w*fn[i](S, pxi);
        }
    }
    result *= drnorm;
}

void Curvature::GS_quad(const MeshElement* const e, const std::array<gp_Pnt, 3>& p, const std::array<gp_Pnt, 3>& px, const std::vector<std::function<double(const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)>>& fn, Eigen::VectorXd& result) const{
    (void)e;
    (void)p;
    (void)px;
    (void)fn;
    (void)result;
}
