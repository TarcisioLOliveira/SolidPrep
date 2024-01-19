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

#ifndef CURVATURE_HPP
#define CURVATURE_HPP

#include <memory>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <gp_Dir.hxx>
#include "logger.hpp"
#include "material.hpp"
#include "element.hpp"
#include "element_factory.hpp"
#include "utils/boundary_nullifier.hpp"
#include "general_solver/mumps_general.hpp"

class Curvature{
    public:

    Curvature(const Material* mat, gp_Dir u, gp_Dir v, gp_Dir w, Eigen::Matrix<double, 2, 2> rot2D, Eigen::Matrix<double, 3, 3> rot3D, const BoundaryMeshElementFactory* elem_info, double V_u, double V_v, double V_w, double M_u, double M_v, double M_w);

    void generate_curvature_3D(const std::vector<std::unique_ptr<MeshNode>>& boundary_nodes, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, size_t reduced_vector_size, size_t full_vector_size);

    inline gp_Pnt get_center() const{
        return gp_Pnt(c_u, c_v, c_w);
    }
    void get_stress_3D(const BoundaryMeshElement* e, double& t_uw, double& t_vw, double& s_w) const;

    private:
    const Material* mat;
    const Eigen::Matrix<double, 2, 2> rot2D;
    const Eigen::Matrix<double, 3, 3> Lek_basis;
    const Eigen::Matrix<double, 3, 3> rot3D;
    const BoundaryMeshElementFactory* elem_info;
    const gp_Dir u, v, w;
    const double V_u, V_v, V_w;
    const double M_u, M_v, M_w;
    size_t reduced_vec_len;

    double Area   = 0;
    double EA     = 0;
    double EA_u   = 0;
    double EA_v   = 0;
    double EI_uu  = 0;
    double EI_vv  = 0;
    double EI_uv  = 0;
    double c_u, c_v, c_w;
    double A, B, C;
    double theta;

    double Kyz = 0;
    double Kxz = 0;

    std::vector<double> F;
    std::vector<double> phi;
    std::vector<double> xi;

    std::vector<double> permute_shear_3D;
    std::vector<double> rotate_X_Z_3D;
    std::vector<double> rotate_X_Z_3D_inv;

    std::vector<double> permute_and_rotate_3D;
    std::vector<double> permute_and_rotate_3D_inv;

    void calculate_stress_field_3D(const std::vector<std::unique_ptr<MeshNode>>& boundary_nodes, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh);

    void base_matrix_upos(general_solver::MUMPSGeneral& M, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const size_t num_nodes, const size_t F_offset) const;
    void base_matrix_id(general_solver::MUMPSGeneral& M, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const size_t num_nodes, const size_t F_offset, const size_t F_phi_offset) const;

    Eigen::VectorXd integrate_surface_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const std::vector<std::function<double(const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)>>& fn) const;
    void GS_tri(const MeshElement* const e, const std::array<gp_Pnt, 3>& p, const std::array<gp_Pnt, 3>& px, const std::vector<std::function<double(const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)>>& fn, Eigen::VectorXd& result) const;
    void GS_quad(const MeshElement* const e, const std::array<gp_Pnt, 3>& p, const std::array<gp_Pnt, 3>& px, const std::vector<std::function<double(const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px)>>& fn, Eigen::VectorXd& result) const;

    Eigen::Matrix<double, 6, 6> get_S_3D(const MeshElement* const e, const gp_Pnt& p) const;
    Eigen::Matrix<double, 6, 6> get_B_3D(const MeshElement* const e, const gp_Pnt& p) const;

    void get_B_tensors_3D(const MeshElement* const e, const gp_Pnt& p, Eigen::Matrix<double, 3, 3>& B4, Eigen::Matrix<double, 3, 2>& B3, Eigen::Matrix<double, 2, 2>& B2) const;

    inline double make_A_base_3D(const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px) const{
        (void)S;
        (void)px;
        return 1;
    }
    inline double make_EA_base_3D(const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px) const{
        (void)px;
        return 1.0/S(2,2);
    }
    inline double make_EA_u_base_3D(const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px) const{
        return px.X()/S(2,2);
    }
    inline double make_EA_v_base_3D(const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px) const{
        return px.Y()/S(2,2);
    }
    inline double make_EA_w_base_3D(const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px) const{
        return px.Z()/S(2,2);
    }
    inline double make_EI_uu_base_3D(const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px) const{
        const double dy = px.Y() - c_v;
        return dy*dy/S(2,2);
    }
    inline double make_EI_vv_base_3D(const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px) const{
        const double dx = px.X() - c_u;
        return dx*dx/S(2,2);
    }
    inline double make_EI_uv_base_3D(const Eigen::Matrix<double, 6, 6>& S, const gp_Pnt& px) const{
        const double dx = px.X() - c_u;
        const double dy = px.Y() - c_v;
        return dx*dy/S(2,2);
    }
};

#endif
