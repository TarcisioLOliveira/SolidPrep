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

void Curvature::generate_curvature_3D(const std::vector<std::unique_ptr<MeshNode>>& boundary_nodes, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const std::vector<utils::LineBoundary>& line_bound, size_t phi_size, size_t psi_size){
    std::array<gp_Pnt, 3> points;

    const auto fn_EA =
        [this](const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EA_base_3D(e, p, px);
        };
    const auto fn_EA_u =
        [this](const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EA_u_base_3D(e, p, px);
        };
    const auto fn_EA_v =
        [this](const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EA_v_base_3D(e, p, px);
        };
    const auto fn_EA_w =
        [this](const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EA_w_base_3D(e, p, px);
        };
    const auto fn_EI_v =
        [this](const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EI_v_base_3D(e, p, px);
        };
    const auto fn_EI_u =
        [this](const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EI_u_base_3D(e, p, px);
        };
    const auto fn_EI_uv =
        [this](const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EI_uv_base_3D(e, p, px);
        };

    this->EA = this->integrate_surface_3D(boundary_mesh, fn_EA);
    this->EA_u = this->integrate_surface_3D(boundary_mesh, fn_EA_u);
    this->EA_v = this->integrate_surface_3D(boundary_mesh, fn_EA_v);
    const double EA_w = this->integrate_surface_3D(boundary_mesh, fn_EA_w);
    this->c_u = EA_u/this->EA;
    this->c_v = EA_v/this->EA;
    this->c_w = EA_w/this->EA;
    this->EI_u = this->integrate_surface_3D(boundary_mesh, fn_EI_u);
    this->EI_v = this->integrate_surface_3D(boundary_mesh, fn_EI_v);
    this->EI_uv = this->integrate_surface_3D(boundary_mesh, fn_EI_uv);
    logger::quick_log("A: ", EA);
    logger::quick_log("EA: ", EA);
    logger::quick_log("c_v: ", c_v);
    logger::quick_log("c_u: ", c_u);
    logger::quick_log("EA_u", EA_u);
    logger::quick_log("EA_v", EA_v);
    logger::quick_log("EI_u", EI_u);
    logger::quick_log("EI_v", EI_v);
    logger::quick_log("EI_uv", EI_uv);

    const Eigen::Vector<double, 2> Mv{M_u, M_v};
    const Eigen::Vector<double, 2> Vv{-V_v, V_u};
    const Eigen::Matrix<double, 2, 2> EI{{ EI_uv,  EI_v},
                                         {-EI_u, -EI_uv}};

    const Eigen::Vector<double, 2> cv = EI.fullPivLu().solve(Mv);
    this->curv_u = cv[0];
    this->curv_v = cv[1];
    const Eigen::Vector<double, 2> dcv = EI.fullPivLu().solve(Vv);
    this->dcurv_u = dcv[0];
    this->dcurv_v = dcv[1];

    this->reduced_vec_len = phi_size;
    this->F.resize(2*psi_size,0);
    this->phi.resize(psi_size,0);
    this->psi_shear.resize(psi_size,0);

    //this->C = -V_w/this->EA;

    if(this->M_u != 0 || this->M_v != 0 || this->M_w != 0){
        this->calculate_torsion(boundary_nodes, boundary_mesh);
    }
    if(this->V_u != 0 || this->V_v != 0){
        this->calculate_shear_3D(boundary_mesh, line_bound);
    }
}

void Curvature::get_shear_in_3D(const BoundaryMeshElement* e, double& s_w, double& t_uw, double& t_vw) const{
    s_w = 0;
    gp_Pnt c = e->get_centroid();
    gp_Pnt c2 = utils::change_point(c, this->rot3D);
    const double x = c.X() - this->c_u;
    const double y = c.Y() - this->c_v;

    Eigen::Vector<double, 2> grad_phi = e->grad_1dof_id(c, this->phi);
    Eigen::Vector<double, 3> grad_F = e->dF_2dof_id(c, this->F);
    t_uw =  grad_phi[1];
    t_vw = -grad_phi[0];

    double s_u = grad_F[1];
    double s_v = grad_F[0];
    double t_uv = -grad_F[2];
    if(this->V_u != 0 || this->V_v != 0){
        Eigen::Vector<double, 2> grad = e->grad_1dof_id(c, this->psi_shear);
        const auto E = this->mat->beam_E_3D(e->parent, c2, this->rot3D);
        t_vw +=  K_uv*grad[1] + E*this->dcurv_v*y*y/2;
        t_uw += -K_uw*grad[0] + E*this->dcurv_u*x*x/2;
    }

    const auto S = this->get_S_3D(e->parent, c2);

    s_w = (A*x + B*y + C)/S(2,2) - (S(0,2)*s_u + S(1,2)*s_v + S(2,3)*t_vw + S(2,4)*t_uw + S(2,5)*t_uv)/S(2,2);
}

void Curvature::calculate_torsion(const std::vector<std::unique_ptr<MeshNode>>& boundary_nodes, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh){
    const size_t num_nodes = this->elem_info->get_nodes_per_element();

    const size_t F_offset = 2*this->reduced_vec_len;
    size_t phi_size = F_offset + this->reduced_vec_len;
    const gp_Pnt center(curv_u, curv_v, 0);

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
        const Eigen::VectorXd N = MULT*2*e->source_1dof();

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
    const size_t F_phi_offset = phi_size;
    phi_size = phi_size + 4;
    solver.initialize_matrix(false, phi_size);

    this->base_matrix_upos(solver, boundary_mesh, num_nodes, F_offset);

    b.resize(phi_size);
    std::fill(b.begin(), b.end(), 0);

    /*
     * Generate the constants which multiply theta and fill the additional
     * rectangles in the matrix.
     *
     */
    double C1 = 0, Cx1 = 0, Cy1 = 0, Ct1 = 0;
    for(const auto& e:boundary_mesh){
        const auto N = e->source_1dof();
        const auto dphi = e->int_grad_phi();
        const auto dphix = e->int_grad_phi_x(center);
        const auto dphiy = e->int_grad_phi_y(center);
        const auto dF = e->int_grad_F();
        const auto dFx = e->int_grad_F_x(center);
        const auto dFy = e->int_grad_F_y(center);
        gp_Pnt c = utils::change_point(e->get_centroid(), this->rot3D);
        const auto S = this->get_S_3D(e->parent, c);

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->u_pos[0];
            if(id1 < 0){
                continue;
            }
            Cx1 += MULT*(
                    (S(1,2)*dFx(0,2*i+0)-0.5*S(2,5)*dFx(2,2*i+0))*phi_tmp[2*id1+0]+
                    (S(0,2)*dFx(1,2*i+1)-0.5*S(2,5)*dFx(2,2*i+1))*phi_tmp[2*id1+1])
                    /S(2,2);
            Cy1 += MULT*(
                    (S(1,2)*dFy(0,2*i+0)-0.5*S(2,5)*dFy(2,2*i+0))*phi_tmp[2*id1+0]+
                    (S(0,2)*dFy(1,2*i+1)-0.5*S(2,5)*dFy(2,2*i+1))*phi_tmp[2*id1+1])
                    /S(2,2);
            C1  += MULT*(
                    (S(1,2)*dF (0,2*i+0)-0.5*S(2,5)*dF (2,2*i+0))*phi_tmp[2*id1+0]+
                    (S(0,2)*dF (1,2*i+1)-0.5*S(2,5)*dF (2,2*i+1))*phi_tmp[2*id1+1])
                    /S(2,2);
            Cx1 += MULT*(-S(2,3)*dphix(0,i)+S(2,4)*dphix(1,i))*phi_tmp[F_offset + id1]/S(2,2);
            Cy1 += MULT*(-S(2,3)*dphiy(0,i)+S(2,4)*dphiy(1,i))*phi_tmp[F_offset + id1]/S(2,2);
            C1  += MULT*(-S(2,3)*dphi (0,i)+S(2,4)*dphi (1,i))*phi_tmp[F_offset + id1]/S(2,2);
            Ct1 += MULT*2*N[i]*phi_tmp[F_offset + id1];

            // A coeffs (eq 2)
            solver.add_value(F_offset + id1, F_phi_offset + 0, -(-MULT*S(2,3)*N[i]/S(2,2)));
            // B coeffs (eq 2)
            solver.add_value(F_offset + id1, F_phi_offset + 1, -(MULT*S(2,4)*N[i]/S(2,2)));

            // Cx0
            solver.add_value(F_phi_offset + 0, 2*id1 + 0, -MULT*(S(1,2)*dFx(0,2*i+0)-0.5*S(2,5)*dFx(2,2*i+0))/S(2,2));
            solver.add_value(F_phi_offset + 0, 2*id1 + 1, -MULT*(S(0,2)*dFx(1,2*i+1)-0.5*S(2,5)*dFx(2,2*i+1))/S(2,2));
            solver.add_value(F_phi_offset + 0, F_offset + id1, -MULT*(-S(2,3)*dphix(0,i)+S(2,4)*dphix(1,i))/S(2,2));
            // Cy0
            solver.add_value(F_phi_offset + 1, 2*id1 + 0, -MULT*(S(1,2)*dFy(0,2*i+0)-0.5*S(2,5)*dFy(2,2*i+0))/S(2,2));
            solver.add_value(F_phi_offset + 1, 2*id1 + 1, -MULT*(S(0,2)*dFy(1,2*i+1)-0.5*S(2,5)*dFy(2,2*i+1))/S(2,2));
            solver.add_value(F_phi_offset + 1, F_offset + id1, -MULT*(-S(2,3)*dphiy(0,i)+S(2,4)*dphiy(1,i))/S(2,2));
            // C0                                                    
            solver.add_value(F_phi_offset + 2, 2*id1 + 0, -MULT*(S(1,2)*dF (0,2*i+0)-0.5*S(2,5)*dF (2,2*i+0))/S(2,2));
            solver.add_value(F_phi_offset + 2, 2*id1 + 1, -MULT*(S(0,2)*dF (1,2*i+1)-0.5*S(2,5)*dF (2,2*i+1))/S(2,2));
            solver.add_value(F_phi_offset + 2, F_offset + id1, -MULT*(-S(2,3)*dphi (0,i)+S(2,4)*dphi (1,i))/S(2,2));
            // Ct0
            solver.add_value(F_phi_offset + 3, F_offset + id1, MULT*2*N[i]);
        }
    }

    // A B square
    solver.add_value(F_phi_offset + 0, F_phi_offset + 0, MULT*EI_v);
    solver.add_value(F_phi_offset + 0, F_phi_offset + 1, MULT*EI_uv);
    solver.add_value(F_phi_offset + 1, F_phi_offset + 0, MULT*EI_uv);
    solver.add_value(F_phi_offset + 1, F_phi_offset + 1, MULT*EI_u);

    // C
    solver.add_value(F_phi_offset + 2, F_phi_offset + 2, MULT*EA);


    // theta
    solver.add_value(F_phi_offset + 0, F_phi_offset + 3, -MULT*Cx1);
    solver.add_value(F_phi_offset + 1, F_phi_offset + 3, -MULT*Cy1);
    solver.add_value(F_phi_offset + 2, F_phi_offset + 3, -MULT*C1 );
    solver.add_value(F_phi_offset + 3, F_phi_offset + 3,  MULT*Ct1);

    // Moments
    b[F_phi_offset + 0] = -MULT*M_v;
    b[F_phi_offset + 1] = -MULT*M_u;
    b[F_phi_offset + 2] = -MULT*V_w;
    b[F_phi_offset + 3] = -MULT*M_w;

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

    solver.compute();

    solver.solve(b);
    phi_tmp = b;

    this->A     = phi_tmp[F_phi_offset + 0];
    this->B     = phi_tmp[F_phi_offset + 1];
    this->C     = phi_tmp[F_phi_offset + 2];
    this->theta = phi_tmp[F_phi_offset + 3];
    logger::quick_log("A", A);
    logger::quick_log("B", B);
    logger::quick_log("C", C);
    logger::quick_log("theta", theta);
    logger::quick_log("");

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
        const long id1 = n->u_pos[0];
        if(id1 < 0){
            continue;
        }
        this->phi[n->id] += phi_tmp[F_offset + id1];
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
        const Eigen::MatrixXd L4 = MULT*e->L4(B4);
        const Eigen::MatrixXd L3 = MULT*e->L3(B3);
        const Eigen::MatrixXd L2 = MULT*e->L2(B2);

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

void Curvature::calculate_shear_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const std::vector<utils::LineBoundary>& line_bound){

    const size_t num_nodes = this->elem_info->get_nodes_per_element();

    const size_t psi_size = this->psi_shear.size();
    Eigen::SparseMatrix<double> M = Eigen::SparseMatrix<double>(psi_size, psi_size);
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::COLAMDOrdering<int>> solver;

    Eigen::VectorXd b;
    gp_Pnt center(c_u, c_v, 0);
    b.resize(psi_size);
    b.fill(0);

    Eigen::Matrix<double, 2, 2> G{{1,0},{0,1}};

    for(const auto& e:boundary_mesh){
        gp_Pnt c = utils::change_point(e->get_centroid(), this->rot3D);
        const auto EG = this->mat->beam_EG_3D(e->parent, c, this->rot3D);
        const auto S = this->mat->S12_S13_3D(e->parent, c, this->rot3D);
        G(0,0) = MULT/EG[2]; // G_uw
        G(1,1) = MULT/EG[1]; // G_uv
        const auto N12 = e->source_1dof(S[1], EG[1], S[0], EG[2], center);
        const auto N = -MULT*EG[0]*(this->dcurv_u*N12[0] + this->dcurv_v*N12[1]);

        const auto M_e = e->diffusion_1dof(G);

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->id;
            //for(size_t j = 0; j <= i; ++j){
            for(size_t j = 0; j < num_nodes; ++j){
                const long id2 = e->nodes[j]->id;
                M.coeffRef(id1, id2) += M_e(i, j);
            }
        }

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->id;
            b[id1] += N[i];
        }
    }

    for(const auto& e:line_bound){
        gp_Pnt pc = utils::change_point(e.parent->get_centroid(), this->rot3D);
        const auto EG = this->mat->beam_EG_3D(e.parent->parent, pc, this->rot3D);
        Eigen::VectorXd N = MULT*EG[0]*e.parent->flow_1dof(dcurv_u, dcurv_v, center, e.normal, {e.edges[0]->point, e.edges[1]->point});

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e.parent->nodes[i]->id;
            b[id1] += N[i];
        }
    }

    M.makeCompressed();
    solver.compute(M);

    logger::log_assert(solver.info() == Eigen::Success, logger::ERROR, "matrix decomposition failed");
    Eigen::VectorXd psi_tmp = solver.solve(b);

    std::copy(psi_tmp.begin(), psi_tmp.end(), this->psi_shear.begin());
    Eigen::Vector<double, 2> psi_grad_int{0,0};
    for(const auto& e:boundary_mesh){
        const auto N = e->int_grad_1dof();

        for(size_t j = 0; j < 2; ++j){
            for(size_t i = 0; i < num_nodes; ++i){
                const long id1 = e->nodes[i]->id;
                psi_grad_int[j] += N(j,i)*this->psi_shear[id1];
            }
        }
    }
    this->K_uv = -(V_v + this->dcurv_u*this->EI_u/2)/psi_grad_int[1];
    this->K_uw =  (V_u - this->dcurv_v*this->EI_v/2)/psi_grad_int[0];
}

double Curvature::integrate_surface_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const std::function<double(const MeshElement* const, const gp_Pnt&, const gp_Pnt&)>& fn) const{
    double result = 0;
    if(this->elem_info->get_shape_type() == Element::Shape::TRI){
        #pragma omp parallel for reduction(+:result)
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
            result += this->GS_tri(e->parent, abs_points, rel_points, fn);
        }
    } else if(this->elem_info->get_shape_type() == Element::Shape::QUAD){

    }
    return result;
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
         {{ B(1,1)  , -B(1,5)  ,  B(0,1)/2},
          {-B(1,5)  ,  B(0,0)  , -B(0,5)/2},
          { B(0,1)/2, -B(0,5)/2,  B(5,5)/4}};

    B3 = Eigen::Matrix<double, 3, 2>
         {{-B(1,3)  , -B(1,4)  },
          {-B(0,3)  ,  B(0,4)  },
          { B(3,5)/2, -B(4,5)/2}};

    B2 = Eigen::Matrix<double, 2, 2>
         {{ B(3,3), -B(3,4)},
          {-B(3,4),  B(5,5)}};
}

double Curvature::GS_tri(const MeshElement* const e, const std::array<gp_Pnt, 3>& p, const std::array<gp_Pnt, 3>& px, const std::function<double(const MeshElement* const, const gp_Pnt&, const gp_Pnt&)>& fn) const{
    double result = 0;
    gp_Vec v1(p[1], p[0]);
    gp_Vec v2(p[2], p[0]);
    const double drnorm = (v1.Crossed(v2)).Magnitude()/2;
    const auto& gsi = utils::GaussLegendreTri<6>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
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
        result += it->w*fn(e, pi, pxi);
    }
    result *= drnorm;

    return result;
}

double Curvature::GS_quad(const MeshElement* const e, const std::array<gp_Pnt, 3>& p, const std::array<gp_Pnt, 3>& px, const std::function<double(const MeshElement* const, const gp_Pnt&, const gp_Pnt&)>& fn) const{
    double result = 0;

    return result;
}
