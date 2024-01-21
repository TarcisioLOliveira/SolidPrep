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
    this->chi.resize(full_vector_size,0);
    this->zeta.resize(full_vector_size,0);

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
    Eigen::Vector<double, 3> grad_F = e->dF_2dof_id(c, this->F);
    Eigen::Vector<double, 2> grad_phi = e->grad_1dof_id(c, this->phi);
    Eigen::Vector<double, 2> grad_xi = e->grad_1dof_id(c, this->xi);
    Eigen::Vector<double, 2> grad_zeta = e->grad_1dof_id(c, this->zeta);
    Eigen::Vector<double, 2> grad_chi = e->grad_1dof_id(c, this->chi);
    //double Fx = 0, Fy = 0;
    //for(size_t i = 0; i < 3; ++i){
    //    const auto& n = e->nodes[i];
    //    Fx += this->F[2*n->id + 0]/3;
    //    Fy += this->F[2*n->id + 1]/3;
    //}
    t_vw = -grad_phi[0] + grad_xi[1] + Kyz;
    t_uw =  grad_phi[1] + grad_xi[0] + Kxz;

    const auto S = this->get_S_3D(e->parent, c2);
    const double s_u = grad_F[1] + grad_chi[0] - grad_zeta[1];
    const double s_v = grad_F[0] - grad_chi[0] + grad_zeta[1];
    const double t_uv = -grad_F[2] + grad_chi[1] + grad_zeta[0];

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
    const size_t F_offset = 2*this->reduced_vec_len; //this->F.size();
    const size_t F_phi_offset = F_offset + this->reduced_vec_len;
    const size_t F_phi_xi_offset = F_phi_offset + this->xi.size();
    const size_t F_phi_xi_chi_offset = F_phi_xi_offset + this->reduced_vec_len;//this->chi.size();
    const size_t F_phi_xi_chi_zeta_offset = F_phi_xi_chi_offset + this->reduced_vec_len;//this->zeta.size();
    constexpr size_t MAX_B = 6;//12;
    const size_t phi_size = F_phi_xi_chi_zeta_offset + MAX_B;

    // Eigen solvers DO NOT WORK CORRECTLY HERE
    general_solver::MUMPSGeneral solver;
    solver.initialize_matrix(false, phi_size);

    this->base_matrix_id(solver, boundary_mesh, num_nodes);

    // A B C theta Kyz Kxz
    std::vector<double> b(phi_size, 0);

    /*
     * Generate the constants which multiply theta and fill the additional
     * rectangles in the matrix.
     */
    Eigen::Matrix<double, 2, 2> I{{1,0},
                                  {0,1}};

    Eigen::Matrix<double, 3, 3> B4;
    Eigen::Matrix<double, 3, 2> B3;
    Eigen::Matrix<double, 2, 2> B2;
    for(const auto& e:boundary_mesh){

        const gp_Pnt centroid = e->get_centroid();
        const gp_Pnt c = utils::change_point(centroid, this->rot3D);
        const auto S = this->get_S_3D(e->parent, c);
        this->get_B_tensors_3D(e->parent, c, B4, B3, B2);
        const Eigen::Matrix<double, 3, 3> aF
            {{S(1,2)/S(2,2), 0, 0},
             {0, S(0,2)/S(2,2), 0},
             {0, 0, -S(2,5)/S(2,2)}};
        const Eigen::Matrix<double, 2, 2> aphi
            {{-S(2,3)/S(2,2), 0},
             {0, S(2,4)/S(2,2)}};
        const Eigen::MatrixXd dFD = e->int_grad_F_D(aF, center);
        const Eigen::MatrixXd dphiD = e->int_grad_phi_D(aphi, center);

        const auto dphi = e->int_grad_phi();
        const auto dF = e->int_grad_F();
        const Eigen::MatrixXd BF = dF.transpose()*B3;
        const Eigen::MatrixXd Bphi = dphi.transpose()*B2;

        const auto N = e->source_1dof();

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->u_pos[0];
            if(id1 > -1){
                // F
                // A
                solver.add_value(2*id1 + 0, F_phi_xi_chi_zeta_offset + 0, dFD(2*i+0,0));
                solver.add_value(2*id1 + 1, F_phi_xi_chi_zeta_offset + 0, dFD(2*i+1,0));

                // B
                solver.add_value(2*id1 + 0, F_phi_xi_chi_zeta_offset + 1, dFD(2*i+0,1));
                solver.add_value(2*id1 + 1, F_phi_xi_chi_zeta_offset + 1, dFD(2*i+1,1));
                                            
                // C
                solver.add_value(2*id1 + 0, F_phi_xi_chi_zeta_offset + 2, dFD(2*i+0,2));
                solver.add_value(2*id1 + 1, F_phi_xi_chi_zeta_offset + 2, dFD(2*i+1,2));
                                            
                // theta                    
                // unused                   
                                            
                // Kyz                      
                solver.add_value(2*id1 + 0, F_phi_xi_chi_zeta_offset + 4, -BF(2*i+0,0));
                solver.add_value(2*id1 + 1, F_phi_xi_chi_zeta_offset + 4, -BF(2*i+1,0));
                                            
                // Kxz                      
                solver.add_value(2*id1 + 0, F_phi_xi_chi_zeta_offset + 5,  BF(2*i+0,1));
                solver.add_value(2*id1 + 1, F_phi_xi_chi_zeta_offset + 5,  BF(2*i+1,1));

                // phi
                 // A
                solver.add_value(F_offset + id1, F_phi_xi_chi_zeta_offset + 0, dphiD(i,0));
                                                 
                 // B                            
                solver.add_value(F_offset + id1, F_phi_xi_chi_zeta_offset + 1, dphiD(i,1));
                                                 
                 // C                            
                solver.add_value(F_offset + id1, F_phi_xi_chi_zeta_offset + 2, dphiD(i,2));
                                                 
                // theta                         
                solver.add_value(F_offset + id1, F_phi_xi_chi_zeta_offset + 3, -2*N[i]);
                                                 
                // Kyz                           
                solver.add_value(F_offset + id1, F_phi_xi_chi_zeta_offset + 4, -Bphi(i,0));
                                                 
                // Kxz                           
                solver.add_value(F_offset + id1, F_phi_xi_chi_zeta_offset + 5,  Bphi(i,1));
            }

        }
    }

    // A B square
    solver.add_value(F_phi_xi_chi_zeta_offset + 0, F_phi_xi_chi_zeta_offset + 0,  EI_vv);
    solver.add_value(F_phi_xi_chi_zeta_offset + 0, F_phi_xi_chi_zeta_offset + 1,  EI_uv);
    solver.add_value(F_phi_xi_chi_zeta_offset + 1, F_phi_xi_chi_zeta_offset + 0,  EI_uv);
    solver.add_value(F_phi_xi_chi_zeta_offset + 1, F_phi_xi_chi_zeta_offset + 1,  EI_uu);

    b[F_phi_xi_chi_zeta_offset + 0] = -M_v;
    b[F_phi_xi_chi_zeta_offset + 1] = -M_u;

    // C
    solver.add_value(F_phi_xi_chi_zeta_offset + 2, F_phi_xi_chi_zeta_offset + 2,  EA);
    b[F_phi_xi_chi_zeta_offset + 2] = -V_w;

    // theta
    b[F_phi_xi_chi_zeta_offset + 3] = -M_w;

    // t_yz
    solver.add_value(F_phi_xi_chi_zeta_offset + 4, F_phi_xi_chi_zeta_offset + 4,  Area);
    b[F_phi_xi_chi_zeta_offset + 4] = V_v;

    // t_xz
    solver.add_value(F_phi_xi_chi_zeta_offset + 5, F_phi_xi_chi_zeta_offset + 5,  Area);
    b[F_phi_xi_chi_zeta_offset + 5] = -V_u;

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

        // Kyz
        solver.add_value(F_phi_xi_chi_zeta_offset + 0, F_phi_xi_chi_zeta_offset + 4, -(S(2,3)*dx*a)/S(2,2));
        solver.add_value(F_phi_xi_chi_zeta_offset + 1, F_phi_xi_chi_zeta_offset + 4, -(S(2,3)*dy*a)/S(2,2));
        solver.add_value(F_phi_xi_chi_zeta_offset + 2, F_phi_xi_chi_zeta_offset + 4, -(S(2,3)*a)/S(2,2));

        // Kxz
        solver.add_value(F_phi_xi_chi_zeta_offset + 0, F_phi_xi_chi_zeta_offset + 5, -(S(2,4)*dx*a)/S(2,2));
        solver.add_value(F_phi_xi_chi_zeta_offset + 1, F_phi_xi_chi_zeta_offset + 5, -(S(2,4)*dy*a)/S(2,2));
        solver.add_value(F_phi_xi_chi_zeta_offset + 2, F_phi_xi_chi_zeta_offset + 5, -(S(2,4)*a)/S(2,2));

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->u_pos[0];

            if(id1 > -1){
                // F
                solver.add_value(F_phi_xi_chi_zeta_offset + 0, 2*id1 + 0, -(S(1,2)*dFx(0,2*i+0)-S(2,5)*dFx(2,2*i+0))/S(2,2));
                solver.add_value(F_phi_xi_chi_zeta_offset + 0, 2*id1 + 1, -(S(0,2)*dFx(1,2*i+1)-S(2,5)*dFx(2,2*i+1))/S(2,2));
                solver.add_value(F_phi_xi_chi_zeta_offset + 1, 2*id1 + 0, -(S(1,2)*dFy(0,2*i+0)-S(2,5)*dFy(2,2*i+0))/S(2,2));
                solver.add_value(F_phi_xi_chi_zeta_offset + 1, 2*id1 + 1, -(S(0,2)*dFy(1,2*i+1)-S(2,5)*dFy(2,2*i+1))/S(2,2));
                solver.add_value(F_phi_xi_chi_zeta_offset + 2, 2*id1 + 0, -(S(1,2)*dF (0,2*i+0)-S(2,5)*dF (2,2*i+0))/S(2,2));
                solver.add_value(F_phi_xi_chi_zeta_offset + 2, 2*id1 + 1, -(S(0,2)*dF (1,2*i+1)-S(2,5)*dF (2,2*i+1))/S(2,2));

                // chi
                solver.add_value(F_phi_xi_chi_zeta_offset + 0, F_phi_xi_offset + id1, -((-S(1,2)+S(0,2))*d1dofx(0,i)+S(2,5)*d1dofx(1,i))/S(2,2));
                solver.add_value(F_phi_xi_chi_zeta_offset + 1, F_phi_xi_offset + id1, -((-S(1,2)+S(0,2))*d1dofy(0,i)+S(2,5)*d1dofy(1,i))/S(2,2));
                solver.add_value(F_phi_xi_chi_zeta_offset + 2, F_phi_xi_offset + id1, -((-S(1,2)+S(0,2))*d1dof (0,i)+S(2,5)*d1dof (1,i))/S(2,2));

                // zeta
                solver.add_value(F_phi_xi_chi_zeta_offset + 0, F_phi_xi_chi_offset + id1, -((S(1,2)-S(0,2))*d1dofx(1,i)+S(2,5)*d1dofx(0,i))/S(2,2));
                solver.add_value(F_phi_xi_chi_zeta_offset + 1, F_phi_xi_chi_offset + id1, -((S(1,2)-S(0,2))*d1dofy(1,i)+S(2,5)*d1dofy(0,i))/S(2,2));
                solver.add_value(F_phi_xi_chi_zeta_offset + 2, F_phi_xi_chi_offset + id1, -((S(1,2)-S(0,2))*d1dof (1,i)+S(2,5)*d1dof (0,i))/S(2,2));

                // // phi
                solver.add_value(F_phi_xi_chi_zeta_offset + 0, F_offset + id1, -(-S(2,3)*d1dofx(0,i)+S(2,4)*d1dofx(1,i))/S(2,2));
                solver.add_value(F_phi_xi_chi_zeta_offset + 1, F_offset + id1, -(-S(2,3)*d1dofy(0,i)+S(2,4)*d1dofy(1,i))/S(2,2));
                solver.add_value(F_phi_xi_chi_zeta_offset + 2, F_offset + id1, -(-S(2,3)*d1dof (0,i)+S(2,4)*d1dof (1,i))/S(2,2));
                solver.add_value(F_phi_xi_chi_zeta_offset + 3, F_offset + id1, (-d1dofx(0,i)-d1dofy(1,i)));
                solver.add_value(F_phi_xi_chi_zeta_offset + 4, F_offset + id1, (-d1dof(0,i)));
                solver.add_value(F_phi_xi_chi_zeta_offset + 5, F_offset + id1, ( d1dof(1,i)));
            }

            const long id2 = e->nodes[i]->id;

            // xi
            solver.add_value(F_phi_xi_chi_zeta_offset + 0, F_phi_offset + id2, -(S(2,3)*d1dofx(1,i)+S(2,4)*d1dofx(0,i))/S(2,2));
            solver.add_value(F_phi_xi_chi_zeta_offset + 1, F_phi_offset + id2, -(S(2,3)*d1dofy(1,i)+S(2,4)*d1dofy(0,i))/S(2,2));
            solver.add_value(F_phi_xi_chi_zeta_offset + 2, F_phi_offset + id2, -(S(2,3)*d1dof (1,i)+S(2,4)*d1dof (0,i))/S(2,2));
            solver.add_value(F_phi_xi_chi_zeta_offset + 3, F_phi_offset + id2, ( d1dofx(1,i)-d1dofy(0,i)));
            solver.add_value(F_phi_xi_chi_zeta_offset + 4, F_phi_offset + id2, ( d1dof(1,i)));
            solver.add_value(F_phi_xi_chi_zeta_offset + 5, F_phi_offset + id2, ( d1dof(0,i)));
        }
    }

    //solver.print_matrix();
    //logger::quick_log("");

    solver.compute();

    double Az = 0, Bz = 0, Cz = 0;
    if(this->V_u != 0 || this->V_v != 0){

        std::vector<double> bz(phi_size, 0);

        bz[F_phi_xi_chi_zeta_offset + 0] = -V_u;
        bz[F_phi_xi_chi_zeta_offset + 1] = V_v;
        solver.solve(bz);

        Az = bz[F_phi_xi_chi_zeta_offset + 0];
        Bz = bz[F_phi_xi_chi_zeta_offset + 1];
        Cz = bz[F_phi_xi_chi_zeta_offset + 2];

        const double Kyz_z = bz[F_phi_xi_chi_zeta_offset + 4];
        const double Kxz_z = bz[F_phi_xi_chi_zeta_offset + 5];

        logger::quick_log("Az", Az);
        logger::quick_log("Bz", Bz);
        logger::quick_log("Cz", Cz);

        #pragma omp parallel for
        for(const auto& n:boundary_nodes){
            const long id1 = n->u_pos[0];
            const long id2 = n->id;
            if(id1 > -1){
                this->F[2*n->id]   = bz[2*id1];
                this->F[2*n->id+1] = bz[2*id1+1];
                this->chi[n->id]   = bz[F_phi_xi_offset + id1];
                this->zeta[n->id]  = bz[F_phi_xi_chi_offset + id1];

                this->phi[n->id]   = bz[F_offset + id1];
            }

            this->xi[n->id]  = bz[F_phi_offset + id2];
        }

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

            const Eigen::MatrixXd grad_phi = e->int_NdN(this->phi);
            const Eigen::MatrixXd grad_F = e->int_NdF(this->F);
            const Eigen::MatrixXd grad_xi = e->int_NdN(this->xi);
            const Eigen::MatrixXd grad_chi = e->int_NdN(this->chi);
            const Eigen::MatrixXd grad_zeta = e->int_NdN(this->zeta);

            for(size_t i = 0; i < num_nodes; ++i){
                const long id1 = e->nodes[i]->u_pos[0];

                const double sx  =  grad_F(i,1) + grad_chi(i,0) - grad_zeta(i,1);
                const double sy  =  grad_F(i,0) - grad_chi(i,0) + grad_zeta(i,1);
                const double txy =  grad_F(i,2) + grad_chi(i,1) + grad_zeta(i,0);
                const double tyz = -grad_phi(i,0) + grad_xi(i,1) + Kyz_z;
                const double txz =  grad_phi(i,1) + grad_xi(i,0) + Kxz_z;

                if(id1 > -1){
                    // chi
                    b[F_phi_xi_offset + id1] += txz;
                    // zeta
                    b[F_phi_xi_chi_offset + id1] += tyz;
                }

                const long id2 = e->nodes[i]->id;

                // xi
                b[F_phi_offset + id2] += (Az*Nx[i] + Bz*Ny[i] + Cz*N[i])/S(2,2);
                b[F_phi_offset + id2] -= (S(0,2)*sx + S(1,2)*sy + S(2,3)*tyz + S(2,4)*txz + S(2,5)*txy)/S(2,2);
            }
        }

        std::fill(F.begin(), F.end(), 0);
        std::fill(phi.begin(), phi.end(), 0);
        std::fill(xi.begin(), xi.end(), 0);
        std::fill(chi.begin(), chi.end(), 0);
        std::fill(zeta.begin(), zeta.end(), 0);
    }

    solver.solve(b);

    this->A     = b[F_phi_xi_chi_zeta_offset + 0];
    this->B     = b[F_phi_xi_chi_zeta_offset + 1];
    this->C     = b[F_phi_xi_chi_zeta_offset + 2];
    this->theta = b[F_phi_xi_chi_zeta_offset + 3];
    this->Kyz   = b[F_phi_xi_chi_zeta_offset + 4];
    this->Kxz   = b[F_phi_xi_chi_zeta_offset + 5];

    /*
     * Save phi_0 and F_0 in the global vector.
     *
     */
    #pragma omp parallel for
    for(const auto& n:boundary_nodes){
        const long id1 = n->u_pos[0];
        const long id2 = n->id;
        if(id1 > -1){
            this->F[2*n->id]   = b[2*id1];
            this->F[2*n->id+1] = b[2*id1+1];
            this->chi[n->id]   = b[F_phi_xi_offset + id1];
            this->zeta[n->id]  = b[F_phi_xi_chi_offset + id1];

            this->phi[n->id]   = b[F_offset + id1];
        }

        this->xi[n->id]  = b[F_phi_offset + id2];
    }

    logger::quick_log("A", A);
    logger::quick_log("B", B);
    logger::quick_log("C", C);
    logger::quick_log("theta", theta);
    logger::quick_log("Kyz", Kyz);
    logger::quick_log("Kxz", Kxz);
    logger::quick_log("");
}
    
void Curvature::base_matrix_id(general_solver::MUMPSGeneral& M, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const size_t num_nodes) const{
    const size_t F_offset = 2*this->reduced_vec_len; //this->F.size();
    const size_t F_phi_offset = F_offset + this->reduced_vec_len;
    const size_t F_phi_xi_offset = F_phi_offset + this->xi.size();
    const size_t F_phi_xi_chi_offset = F_phi_xi_offset + this->reduced_vec_len;//this->chi.size();

    Eigen::Matrix<double, 3, 3> B4;
    Eigen::Matrix<double, 3, 2> B3;
    Eigen::Matrix<double, 2, 2> B2;
    Eigen::Matrix<double, 2, 2> A{{1, 0},{0, 1}};

    std::vector<double> L4v(2*num_nodes*2*num_nodes,0);
    std::vector<double> L3v(2*num_nodes*num_nodes,0);
    std::vector<double> L3Tv(2*num_nodes*num_nodes,0);
    std::vector<double> L2v(num_nodes*num_nodes,0);

    std::vector<double> L3xiv(2*num_nodes*num_nodes,0);
    std::vector<double> L2xiv(num_nodes*num_nodes,0);

    std::vector<double> L4chiv(2*num_nodes*num_nodes,0);
    std::vector<double> L3Tchiv(num_nodes*num_nodes,0);

    std::vector<double> L4zetav(2*num_nodes*num_nodes,0);
    std::vector<double> L3Tzetav(num_nodes*num_nodes,0);

    std::vector<double> Lv(num_nodes*num_nodes,0);

    std::vector<long> pos_L4(num_nodes*2,0);
    std::vector<long> pos_L2(num_nodes,0);
    std::vector<long> pos_L2xi(num_nodes,0);
    std::vector<long> pos_L2chi(num_nodes,0);
    std::vector<long> pos_L2zeta(num_nodes,0);

    for(const auto& e:boundary_mesh){
        gp_Pnt c = utils::change_point(e->get_centroid(), this->rot3D);
        this->get_B_tensors_3D(e->parent, c, B4, B3, B2);
        const Eigen::MatrixXd L4 = e->L4(B4);
        const Eigen::MatrixXd L3 = e->L3(B3);
        const Eigen::MatrixXd L2 = e->L2(B2);

        const Eigen::MatrixXd L3xi = e->L3xi(B3);
        const Eigen::MatrixXd L2xi = e->L2xi(B2);

        const Eigen::MatrixXd L4chi = e->L4chi(B4);
        const Eigen::MatrixXd L3Tchi = e->L3Tchi(B3);

        const Eigen::MatrixXd L4zeta = e->L4zeta(B4);
        const Eigen::MatrixXd L3Tzeta = e->L3Tzeta(B3);

        const Eigen::MatrixXd Lz = e->diffusion_1dof(A);

        for(size_t i = 0; i < 2*num_nodes; ++i){
            for(size_t j = 0; j < 2*num_nodes; ++j){
                L4v[i*2*num_nodes + j] = L4(i,j);
            }
        }
        for(size_t i = 0; i < 2*num_nodes; ++i){
            for(size_t j = 0; j < num_nodes; ++j){
                L3v[i*num_nodes + j] = L3(i,j);
                L3Tv[j*2*num_nodes + i] = L3(i,j);

                L3xiv[i*num_nodes + j] = L3xi(i,j);

                L4chiv[i*num_nodes + j] = L4chi(i,j);
                L4zetav[i*num_nodes + j] = L4zeta(i,j);
            }
        }
        for(size_t i = 0; i < num_nodes; ++i){
            for(size_t j = 0; j < num_nodes; ++j){
                L2v[i*num_nodes + j] = L2(i,j);

                L2xiv[i*num_nodes + j] = L2xi(i,j);

                L3Tchiv[i*num_nodes + j] = L3Tchi(i,j);
                L3Tzetav[i*num_nodes + j] = L3Tzeta(i,j);

                Lv[i*num_nodes + j] = Lz(i,j);
            }
        }

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->u_pos[0];
            if(id1 > -1){
                pos_L4[2*i] = 2*id1;
                pos_L4[2*i+1] = 2*id1+1;
                pos_L2[i] = F_offset + id1;
                pos_L2chi[i] = F_phi_xi_offset + id1;
                pos_L2zeta[i] = F_phi_xi_chi_offset + id1;
            } else {
                pos_L4[2*i] = -1;
                pos_L4[2*i+1] = -1;
                pos_L2[i] = -1;
                pos_L2chi[i] = -1;
                pos_L2zeta[i] = -1;
            }
        }

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->id;
            pos_L2xi[i] = F_phi_offset + id1;
        }

        M.add_element(L4v, pos_L4);
        M.add_element(L3v, pos_L4, pos_L2);
        M.add_element(L3Tv, pos_L2, pos_L4);
        M.add_element(L2v, pos_L2);

        M.add_element(L3xiv, pos_L4, pos_L2xi);
        M.add_element(L2xiv, pos_L2, pos_L2xi);

        M.add_element(L4chiv, pos_L4, pos_L2chi);
        M.add_element(L3Tchiv, pos_L2, pos_L2chi);

        M.add_element(L4zetav, pos_L4, pos_L2zeta);
        M.add_element(L3Tzetav, pos_L2, pos_L2zeta);

        M.add_element(Lv, pos_L2xi);
        M.add_element(Lv, pos_L2chi);
        M.add_element(Lv, pos_L2zeta);
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
