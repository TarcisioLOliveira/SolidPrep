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
    this->rotate_X_Z_3D_inv = utils::basis_tensor_3D_inv_T(Lek_basis);

    this->permute_and_rotate_3D = std::vector<double>(36, 0);
    this->permute_and_rotate_3D_inv = std::vector<double>(36, 0);
    const size_t N = 6;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, permute_shear_3D.data(), N, rotate_X_Z_3D.data(), N, 0, permute_and_rotate_3D.data(), N);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, rotate_X_Z_3D_inv.data(), N, permute_shear_3D.data(), N, 0, permute_and_rotate_3D_inv.data(), N);
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
    this->F.resize(2*full_vector_size,0);
    this->phi.resize(full_vector_size,0);
    this->xi.resize(full_vector_size,0);

    this->calculate_stress_field_3D(boundary_nodes, boundary_mesh);
}

void Curvature::get_stress_3D(const BoundaryMeshElement* e, double& t_uw, double& t_vw, double& s_w) const{
    s_w = 0;
    gp_Pnt c = e->get_centroid();
    gp_Pnt c2 = utils::change_point(c, this->rot3D);
    const double x = c.X() - this->c_u;
    const double y = c.Y() - this->c_v;

    t_vw = 0;
    t_uw = 0;
    Eigen::Vector<double, 2> grad_phi = e->grad_1dof_id(c, this->phi);
    Eigen::Vector<double, 2> grad_xi = e->grad_1dof_id(c, this->xi);
    Eigen::Vector<double, 3> grad_F = e->dF_2dof_id(c, this->F);
    t_vw = -grad_phi[0] + grad_xi[1] + Kyz;
    t_uw =  grad_phi[1] + grad_xi[0] + Kxz;

    const auto S = this->get_S_3D(e->parent, c2);
    double s_u = grad_F[1];
    double s_v = grad_F[0];
    double t_uv = -grad_F[2];

    s_w = (A*x + B*y + C)/S(2,2) - (S(0,2)*s_u + S(1,2)*s_v + S(2,3)*t_vw + S(2,4)*t_uw + S(2,5)*t_uv)/S(2,2);
}

void Curvature::calculate_stress_field_3D(const std::vector<std::unique_ptr<MeshNode>>& boundary_nodes, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh){
    const size_t num_nodes = this->elem_info->get_nodes_per_element();

    const gp_Pnt center(c_u, c_v, 0);

    /* 
     * Obtain F, phi, xi, etc.
     *
     * Matrix is slightly larger to acomodate theta, A and B.
     * Unfortunately, such addition makes it necessary to redo M and refactorize
     * it.
     *
     * It only isn't necessary when S is at most orthotropic, but it is not
     * really necessary to optimize for it as the solver would already ignore
     * the zero-filled rectangles.
     */
    const size_t F_offset = 2*this->reduced_vec_len;
    const size_t F_phi_offset = F_offset + this->phi.size();
    const size_t F_phi_xi_offset = F_phi_offset + this->xi.size();
    const size_t phi_size = F_phi_xi_offset;
    // Eigen solvers DO NOT WORK CORRECTLY HERE
    general_solver::MUMPSGeneral solver;
    solver.initialize_matrix(false, phi_size);

    this->base_matrix_id(solver, boundary_mesh, num_nodes, F_offset, F_phi_offset);

    //solver.print_matrix();
    //logger::quick_log("");

    solver.compute();

    // 0 A B C theta Kyz Kxz
    std::vector<std::vector<double>> b(7);
    for(auto& v:b){
        v.resize(phi_size,0);
    }

    // Az Bz
    Eigen::Matrix<double, 2, 2> EIM{{EI_vv, EI_uv},
                                    {EI_uv, EI_uu}};
    Eigen::Vector<double, 2> VV{-V_u, V_v};
    Eigen::Vector<double, 2> Vz = EIM.fullPivLu().solve(VV);

    const double Az = Vz[0];
    const double Bz = Vz[1];

    logger::quick_log("Az", Az);
    logger::quick_log("Bz", Bz);

    /*
     * Generate the constants which multiply theta and fill the additional
     * rectangles in the matrix.
     */
    Eigen::Matrix<double, 2, 2> I{{1,0},
                                  {0,1}};

    // 0
    size_t B_NUM = 0;
    if(Az != 0 || Bz != 0){
        for(const auto& e:boundary_mesh){

            const gp_Pnt centroid = e->get_centroid();
            const gp_Pnt c = utils::change_point(centroid, this->rot3D);
            const auto S = this->get_S_3D(e->parent, c);

            const auto NAzBz = e->int_N_AzBz(center, Az/S(2,2), Bz/S(2,2));

            for(size_t i = 0; i < num_nodes; ++i){
                const long id2 = e->nodes[i]->id;

                b[B_NUM][F_phi_offset + id2] += NAzBz[i];
            }
        }

        solver.solve(b[B_NUM]);
    }

    // A
    B_NUM = 1;
    for(const auto& e:boundary_mesh){
        Eigen::Matrix<double, 3, 3> B4;
        Eigen::Matrix<double, 3, 2> B3;
        Eigen::Matrix<double, 2, 2> B2;

        const gp_Pnt centroid = e->get_centroid();
        const gp_Pnt c = utils::change_point(centroid, this->rot3D);
        const auto S = this->get_S_3D(e->parent, c);
        this->get_B_tensors_3D(e->parent, c, B4, B3, B2);
        const Eigen::Matrix<double, 3, 3> aF
            {{S(1,2)/S(2,2), 0, 0},
             {0, S(0,2)/S(2,2), 0},
             {0, 0, S(2,5)/S(2,2)}};
        const Eigen::Matrix<double, 2, 2> aphi
            {{S(2,3)/S(2,2), 0},
             {0, S(2,4)/S(2,2)}};
        const Eigen::MatrixXd dFD = e->int_grad_F_D(aF, center);
        const Eigen::MatrixXd dphiD = e->int_grad_phi_D(aphi, center);

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->u_pos[0];
            const long id2 = e->nodes[i]->id;

            if(id1 >= 0){
                b[B_NUM][2*id1 + 0] -= dFD(2*i+0,0);
                b[B_NUM][2*id1 + 1] -= dFD(2*i+1,0);
            }
            b[B_NUM][F_offset + id2] -= dphiD(i,0);
        }
    }
    solver.solve(b[B_NUM]);

    // B
    B_NUM = 2;
    for(const auto& e:boundary_mesh){
        Eigen::Matrix<double, 3, 3> B4;
        Eigen::Matrix<double, 3, 2> B3;
        Eigen::Matrix<double, 2, 2> B2;

        const gp_Pnt centroid = e->get_centroid();
        const gp_Pnt c = utils::change_point(centroid, this->rot3D);
        const auto S = this->get_S_3D(e->parent, c);
        this->get_B_tensors_3D(e->parent, c, B4, B3, B2);
        const Eigen::Matrix<double, 3, 3> aF
            {{S(1,2)/S(2,2), 0, 0},
             {0, S(0,2)/S(2,2), 0},
             {0, 0, S(2,5)/S(2,2)}};
        const Eigen::Matrix<double, 2, 2> aphi
            {{S(2,3)/S(2,2), 0},
             {0, S(2,4)/S(2,2)}};
        const Eigen::MatrixXd dFD = e->int_grad_F_D(aF, center);
        const Eigen::MatrixXd dphiD = e->int_grad_phi_D(aphi, center);

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->u_pos[0];
            const long id2 = e->nodes[i]->id;

            if(id1 >= 0){
                b[B_NUM][2*id1 + 0] -= dFD(2*i+0,1);
                b[B_NUM][2*id1 + 1] -= dFD(2*i+1,1);
            }
            b[B_NUM][F_offset + id2] -= dphiD(i,1);
        }
    }
    solver.solve(b[B_NUM]);

    // C
    B_NUM = 3;
    for(const auto& e:boundary_mesh){
        Eigen::Matrix<double, 3, 3> B4;
        Eigen::Matrix<double, 3, 2> B3;
        Eigen::Matrix<double, 2, 2> B2;

        const gp_Pnt centroid = e->get_centroid();
        const gp_Pnt c = utils::change_point(centroid, this->rot3D);
        const auto S = this->get_S_3D(e->parent, c);
        this->get_B_tensors_3D(e->parent, c, B4, B3, B2);
        const Eigen::Matrix<double, 3, 3> aF
            {{S(1,2)/S(2,2), 0, 0},
             {0, S(0,2)/S(2,2), 0},
             {0, 0, S(2,5)/S(2,2)}};
        const Eigen::Matrix<double, 2, 2> aphi
            {{S(2,3)/S(2,2), 0},
             {0, S(2,4)/S(2,2)}};
        const Eigen::MatrixXd dFD = e->int_grad_F_D(aF, center);
        const Eigen::MatrixXd dphiD = e->int_grad_phi_D(aphi, center);

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->u_pos[0];
            const long id2 = e->nodes[i]->id;

            if(id1 >= 0){
                b[B_NUM][2*id1 + 0] -= dFD(2*i+0,2);
                b[B_NUM][2*id1 + 1] -= dFD(2*i+1,2);
            }
            b[B_NUM][F_offset + id2] -= dphiD(i,2);
        }
    }
    solver.solve(b[B_NUM]);

    // theta
    B_NUM = 4;
    for(const auto& e:boundary_mesh){
        const auto N = e->source_1dof();

        for(size_t i = 0; i < num_nodes; ++i){
            const long id2 = e->nodes[i]->id;

            b[B_NUM][F_offset + id2] += 2*N[i];
        }
    }
    solver.solve(b[B_NUM]);

    // Kyz
    B_NUM = 5;
    for(const auto& e:boundary_mesh){
        Eigen::Matrix<double, 3, 3> B4;
        Eigen::Matrix<double, 3, 2> B3;
        Eigen::Matrix<double, 2, 2> B2;

        const gp_Pnt centroid = e->get_centroid();
        const gp_Pnt c = utils::change_point(centroid, this->rot3D);
        const auto dphi = e->int_grad_phi();
        const auto dF = e->int_grad_F();
        this->get_B_tensors_3D(e->parent, c, B4, B3, B2);
        const Eigen::MatrixXd BF = dF.transpose()*B3;
        const Eigen::MatrixXd Bphi = dphi.transpose()*B2;

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->u_pos[0];
            const long id2 = e->nodes[i]->id;

            if(id1 >= 0){
                b[B_NUM][2*id1 + 0] -= BF(2*i+0,0);
                b[B_NUM][2*id1 + 1] -= BF(2*i+1,0);
            }
            b[B_NUM][F_offset + id2] -= Bphi(i,0);
        }
    }
    solver.solve(b[B_NUM]);

    // Kxz
    B_NUM = 6;
    for(const auto& e:boundary_mesh){
        Eigen::Matrix<double, 3, 3> B4;
        Eigen::Matrix<double, 3, 2> B3;
        Eigen::Matrix<double, 2, 2> B2;

        const gp_Pnt centroid = e->get_centroid();
        const gp_Pnt c = utils::change_point(centroid, this->rot3D);
        const auto dphi = e->int_grad_phi();
        const auto dF = e->int_grad_F();
        this->get_B_tensors_3D(e->parent, c, B4, B3, B2);
        const Eigen::MatrixXd BF = dF.transpose()*B3;
        const Eigen::MatrixXd Bphi = dphi.transpose()*B2;

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->u_pos[0];
            const long id2 = e->nodes[i]->id;

            if(id1 >= 0){
                b[B_NUM][2*id1 + 0] -= BF(2*i+0,1);
                b[B_NUM][2*id1 + 1] -= BF(2*i+1,1);
            }
            b[B_NUM][F_offset + id2] -= Bphi(i,1);
        }
    }
    solver.solve(b[B_NUM]);

    Eigen::Matrix<double, 6, 6> Mat;
    Mat.fill(0);
    Eigen::Vector<double, 6> Vec;
    Vec.fill(0);

    // A B square
    Mat(0, 0) = EI_vv;
    Mat(0, 1) = EI_uv;
    Mat(1, 0) = EI_uv;
    Mat(1, 1) = EI_uu;

    Vec[0] = -M_v;
    Vec[1] = -M_u;

    // C
    Mat(2,2) = EA;
    Vec[2] = -V_w;

    // theta
    Vec[3] = -M_w;

    // t_yz
    Mat(4,4) = Area;
    Vec[4] = V_v;
    // t_xz
    Mat(5,5) = Area;
    Vec[5] = -V_u;

    for(const auto& e:boundary_mesh){
        const auto N = e->source_1dof();
        const auto Nx = e->int_N_x(center);
        const auto Ny = e->int_N_y(center);
        const auto M_e = e->diffusion_1dof(I);
        const auto d1dof = e->int_grad_phi();
        const auto d1dofx = e->int_grad_phi_x(center);
        const auto d1dofy = e->int_grad_phi_y(center);
        //const auto dxi = e->int_grad_xi();
        //const auto dxix = e->int_grad_xi_x(center);
        //const auto dxiy = e->int_grad_xi_y(center);
        const auto dF = e->int_grad_F();
        const auto dFx = e->int_grad_F_x(center);
        const auto dFy = e->int_grad_F_y(center);
        const gp_Pnt centroid = e->get_centroid();
        const gp_Pnt c = utils::change_point(centroid, this->rot3D);
        const auto S = this->get_S_3D(e->parent, c);

        const double dx = centroid.X() - center.X();
        const double dy = centroid.Y() - center.Y();
        const double a = e->get_area();

        Mat(0,4) += -(S(2,3)*dx*a)/S(2,2);
        Mat(1,4) += -(S(2,3)*dy*a)/S(2,2);
        Mat(2,4) += -(S(2,3)*a)/S(2,2);

        Mat(0,5) += -(S(2,4)*dx*a)/S(2,2);
        Mat(1,5) += -(S(2,4)*dy*a)/S(2,2);
        Mat(2,5) += -(S(2,4)*a)/S(2,2);

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->u_pos[0];
            if(id1 >= 0){
                Vec[0] -= -(S(1,2)*dFx(0,2*i+0)-S(2,5)*dFx(2,2*i+0))/S(2,2)*b[0][2*id1 + 0];
                Vec[0] -= -(S(0,2)*dFx(1,2*i+1)-S(2,5)*dFx(2,2*i+1))/S(2,2)*b[0][2*id1 + 1];
                Vec[1] -= -(S(1,2)*dFy(0,2*i+0)-S(2,5)*dFy(2,2*i+0))/S(2,2)*b[0][2*id1 + 0];
                Vec[1] -= -(S(0,2)*dFy(1,2*i+1)-S(2,5)*dFy(2,2*i+1))/S(2,2)*b[0][2*id1 + 1];
                Vec[2] -= -(S(1,2)*dF (0,2*i+0)-S(2,5)*dF (2,2*i+0))/S(2,2)*b[0][2*id1 + 0];
                Vec[2] -= -(S(0,2)*dF (1,2*i+1)-S(2,5)*dF (2,2*i+1))/S(2,2)*b[0][2*id1 + 1];
            }

            const long id2 = e->nodes[i]->id;

            Vec[0] -= -(-S(2,3)*d1dofx(0,i)+S(2,4)*d1dofx(1,i))*b[0][F_offset + id2]/S(2,2);
            Vec[1] -= -(-S(2,3)*d1dofy(0,i)+S(2,4)*d1dofy(1,i))*b[0][F_offset + id2]/S(2,2);
            Vec[2] -= -(-S(2,3)*d1dof (0,i)+S(2,4)*d1dof (1,i))*b[0][F_offset + id2]/S(2,2);
            Vec[3] -= (-d1dofx(0,i)-d1dofy(1,i))*b[0][F_offset + id2];
            Vec[4] -= (-d1dof(0,i))*b[0][F_offset + id2];
            Vec[5] -= ( d1dof(1,i))*b[0][F_offset + id2];

            Vec[0] -= -(S(2,3)*d1dofx(1,i)+S(2,4)*d1dofx(0,i))*b[0][F_phi_offset + id2]/S(2,2);
            Vec[1] -= -(S(2,3)*d1dofy(1,i)+S(2,4)*d1dofy(0,i))*b[0][F_phi_offset + id2]/S(2,2);
            Vec[2] -= -(S(2,3)*d1dof (1,i)+S(2,4)*d1dof (0,i))*b[0][F_phi_offset + id2]/S(2,2);
            Vec[3] -= ( d1dofx(1,i)-d1dofy(0,i))*b[0][F_phi_offset + id2];
            Vec[4] -= ( d1dof(1,i))*b[0][F_phi_offset + id2];
            Vec[5] -= ( d1dof(0,i))*b[0][F_phi_offset + id2];
        }
        for(size_t j = 0; j < 6; ++j){
            for(size_t i = 0; i < num_nodes; ++i){
                const long id1 = e->nodes[i]->u_pos[0];
                if(id1 >= 0){
                    Mat(0,j) += -(S(1,2)*dFx(0,2*i+0)-S(2,5)*dFx(2,2*i+0))/S(2,2)*b[j+1][2*id1 + 0];
                    Mat(0,j) += -(S(0,2)*dFx(1,2*i+1)-S(2,5)*dFx(2,2*i+1))/S(2,2)*b[j+1][2*id1 + 1];
                    Mat(1,j) += -(S(1,2)*dFy(0,2*i+0)-S(2,5)*dFy(2,2*i+0))/S(2,2)*b[j+1][2*id1 + 0];
                    Mat(1,j) += -(S(0,2)*dFy(1,2*i+1)-S(2,5)*dFy(2,2*i+1))/S(2,2)*b[j+1][2*id1 + 1];
                    Mat(2,j) += -(S(1,2)*dF (0,2*i+0)-S(2,5)*dF (2,2*i+0))/S(2,2)*b[j+1][2*id1 + 0];
                    Mat(2,j) += -(S(0,2)*dF (1,2*i+1)-S(2,5)*dF (2,2*i+1))/S(2,2)*b[j+1][2*id1 + 1];
                }
                const long id2 = e->nodes[i]->id;

                Mat(0,j) += -(-S(2,3)*d1dofx(0,i)+S(2,4)*d1dofx(1,i))*b[j+1][F_offset + id2]/S(2,2);
                Mat(1,j) += -(-S(2,3)*d1dofy(0,i)+S(2,4)*d1dofy(1,i))*b[j+1][F_offset + id2]/S(2,2);
                Mat(2,j) += -(-S(2,3)*d1dof (0,i)+S(2,4)*d1dof (1,i))*b[j+1][F_offset + id2]/S(2,2);
                Mat(3,j) += (-d1dofx(0,i)-d1dofy(1,i))*b[j+1][F_offset + id2];
                Mat(4,j) += (-d1dof(0,i))*b[j+1][F_offset + id2];
                Mat(5,j) += ( d1dof(1,i))*b[j+1][F_offset + id2];

                Mat(0,j) += -(S(2,3)*d1dofx(1,i)+S(2,4)*d1dofx(0,i))*b[j+1][F_phi_offset + id2]/S(2,2);
                Mat(1,j) += -(S(2,3)*d1dofy(1,i)+S(2,4)*d1dofy(0,i))*b[j+1][F_phi_offset + id2]/S(2,2);
                Mat(2,j) += -(S(2,3)*d1dof (1,i)+S(2,4)*d1dof (0,i))*b[j+1][F_phi_offset + id2]/S(2,2);
                Mat(3,j) += ( d1dofx(1,i)-d1dofy(0,i))*b[j+1][F_phi_offset + id2];
                Mat(4,j) += ( d1dof(1,i))*b[j+1][F_phi_offset + id2];
                Mat(5,j) += ( d1dof(0,i))*b[j+1][F_phi_offset + id2];
            }
        }
    }

    logger::quick_log(Mat);
    logger::quick_log("");
    logger::quick_log(Vec);
    Eigen::Vector<double, 6> Res = Mat.fullPivLu().solve(Vec);

    this->A     = Res[0];
    this->B     = Res[1];
    this->C     = Res[2];
    this->theta = Res[3];
    this->Kyz   = Res[4];
    this->Kxz   = Res[5];

    /*
     * Save phi_0 and F_0 in the global vector.
     *
     */
    for(const auto& n:boundary_nodes){
        const long id1 = n->u_pos[0];
        if(id1 < 0){
            continue;
        }
        this->F[2*n->id] += b[0][2*id1];
        this->F[2*n->id+1] += b[0][2*id1+1];
        for(size_t i = 0; i < 6; ++i){
            this->F[2*n->id] += Res[i]*b[i+1][2*id1];
            this->F[2*n->id+1] += Res[i]*b[i+1][2*id1+1];
        }
    }
    for(const auto& n:boundary_nodes){
        const long id1 = n->id;
        this->phi[n->id] += b[0][F_offset + id1];
        for(size_t i = 0; i < 6; ++i){
            this->phi[n->id] += Res[i]*b[i+1][F_offset + id1];
        }
    }
    for(const auto& n:boundary_nodes){
        const long id1 = n->id;
        this->xi[n->id] += b[0][F_phi_offset + id1];
        for(size_t i = 0; i < 6; ++i){
            this->xi[n->id] += Res[i]*b[i+1][F_phi_offset + id1];
        }
    }

    logger::quick_log("A", A);
    logger::quick_log("B", B);
    logger::quick_log("C", C);
    logger::quick_log("theta", theta);
    logger::quick_log("Kyz", Kyz);
    logger::quick_log("Kxz", Kxz);
    logger::quick_log("");
}
    
void Curvature::base_matrix_id(general_solver::MUMPSGeneral& M, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const size_t num_nodes, const size_t F_offset, const size_t F_phi_offset) const{
    Eigen::Matrix<double, 3, 3> B4;
    Eigen::Matrix<double, 3, 2> B3;
    Eigen::Matrix<double, 2, 2> B2;
    Eigen::Matrix<double, 2, 2> A{{1, 0},{0, 1}};

    std::vector<double> L4v(2*num_nodes*2*num_nodes,0);
    std::vector<double> L3v(2*num_nodes*num_nodes,0);
    std::vector<double> L3vT(2*num_nodes*num_nodes,0);
    std::vector<double> L2v(num_nodes*num_nodes,0);

    std::vector<double> L3zv(2*num_nodes*num_nodes,0);
    std::vector<double> L2zv(num_nodes*num_nodes,0);
    std::vector<double> Lzv(num_nodes*num_nodes,0);

    std::vector<long> pos_L4(num_nodes*2,0);
    std::vector<long> pos_L2(num_nodes,0);
    std::vector<long> pos_L2z(num_nodes,0);

    for(const auto& e:boundary_mesh){
        gp_Pnt c = utils::change_point(e->get_centroid(), this->rot3D);
        this->get_B_tensors_3D(e->parent, c, B4, B3, B2);
        const Eigen::MatrixXd L4 = e->L4(B4);
        const Eigen::MatrixXd L3 = e->L3(B3);
        const Eigen::MatrixXd L2 = e->L2(B2);

        const Eigen::MatrixXd L3z = e->L3z(B3);
        const Eigen::MatrixXd L2z = e->L2z(B2);

        const Eigen::MatrixXd Lz = e->diffusion_1dof(A);

        for(size_t i = 0; i < 2*num_nodes; ++i){
            for(size_t j = 0; j < 2*num_nodes; ++j){
                L4v[i*2*num_nodes + j] = L4(i,j);
            }
        }
        for(size_t i = 0; i < 2*num_nodes; ++i){
            for(size_t j = 0; j < num_nodes; ++j){
                L3v[i*num_nodes + j] = L3(i,j);
                L3vT[j*2*num_nodes + i] = L3(i,j);

                L3zv[i*num_nodes + j] = L3z(i,j);
            }
        }
        for(size_t i = 0; i < num_nodes; ++i){
            for(size_t j = 0; j < num_nodes; ++j){
                L2v[i*num_nodes + j] = L2(i,j);

                L2zv[i*num_nodes + j] = L2z(i,j);
                Lzv[i*num_nodes + j] = Lz(i,j);
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
            pos_L2z[i] = F_phi_offset + id1;
        }

        M.add_element(L4v, pos_L4);
        M.add_element(L3v, pos_L4, pos_L2);
        M.add_element(L3vT, pos_L2, pos_L4);
        M.add_element(L2v, pos_L2);

        M.add_element(L3zv, pos_L4, pos_L2z);
        M.add_element(L2zv, pos_L2, pos_L2z);
        M.add_element(Lzv, pos_L2z);
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
    this->mat->rotate_S_3D(S, this->permute_and_rotate_3D_inv);

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
         {{ B(1,1),  B(0,1), -B(1,5)},
          { B(0,1),  B(0,0), -B(0,5)},
          {-B(1,5), -B(0,5),  B(5,5)}};

    B3 = Eigen::Matrix<double, 3, 2>
         {{-B(1,3)  ,  B(1,4)  },
          {-B(0,3)  ,  B(0,4)  },
          { B(3,5)  , -B(4,5)  }};

    B2 = Eigen::Matrix<double, 2, 2>
         {{ B(3,3), -B(3,4)},
          {-B(3,4),  B(4,4)}};
}

void Curvature::GS_tri(const MeshElement* const e, const std::array<gp_Pnt, 3>& p, const std::array<gp_Pnt, 3>& px, const std::vector<std::function<double(const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)>>& fn, Eigen::VectorXd& result) const{
    gp_Vec v1(p[1], p[0]);
    gp_Vec v2(p[2], p[0]);
    result.fill(0);
    const gp_Pnt centroid(
            (p[0].X() + p[1].X() + p[2].X())/3.0,
            (p[0].Y() + p[1].Y() + p[2].Y())/3.0,
            (p[0].Z() + p[1].Z() + p[2].Z())/3.0);
    const auto S = this->get_S_3D(e, centroid);
    const double drnorm = (v1.Crossed(v2)).Magnitude()/2;
    const auto& gsi = utils::GaussLegendreTri<6>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        for(size_t i = 0; i < fn.size(); ++i){
            gp_Pnt pxi{
                it->a*px[0].X() + it->b*px[1].X() + it->c*px[2].X(),
                it->a*px[0].Y() + it->b*px[1].Y() + it->c*px[2].Y(),
                it->a*px[0].Z() + it->b*px[1].Z() + it->c*px[2].Z()
            };
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
