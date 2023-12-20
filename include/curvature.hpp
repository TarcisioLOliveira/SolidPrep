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

    void generate_curvature_3D(const std::vector<std::unique_ptr<MeshNode>>& boundary_nodes, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const std::vector<utils::LineBoundary>& line_bound, size_t phi_size, size_t psi_size);

    inline std::array<double, 2> get_curvatures() const{
        return {curv_u, curv_v};
    }
    inline std::array<double, 2> get_curvatures_dx() const{
        return {dcurv_u, dcurv_v};
    }
    inline gp_Pnt get_center() const{
        return gp_Pnt(c_u, c_v, c_w);
    }
    void get_shear_in_3D(const BoundaryMeshElement* e, double& s_w, double& t_uw, double& t_vw) const;

    private:
    const Material* mat;
    const Eigen::Matrix<double, 2, 2> rot2D;
    const Eigen::Matrix<double, 3, 3> Lek_basis;
    const Eigen::Matrix<double, 3, 3> rot3D;
    const BoundaryMeshElementFactory* elem_info;
    const gp_Dir u, v, w;
    const double V_u, V_v, V_w;
    const double M_u, M_v, M_w;
    const double line_mesh_size = 0.5;

    double EA;
    double EA_u;
    double EA_v;
    double EI_uu;
    double EI_vv;
    double EI_uv;
    double EI_uuu;
    double EI_vvv;
    double EI_uvv;
    double EI_uuv;
    double c_u, c_v, c_w;
    double curv_u;
    double curv_v;
    double dcurv_u;
    double dcurv_v;
    double A, B, C;

    double Kyz_1;
    double Kyz_2;
    double Kyz_3;
    double Kxz_1;
    double Kxz_2;
    double Kxz_3;

    double Lyz_vv = 0;
    double Lyz_u  = 0;
    double Lyz_uu = 0;
    double Lxz_uu = 0;
    double Lxz_v  = 0;
    double Lxz_vv = 0;

    double Syz_vv = 0;
    double Syz_uu = 0;
    double Sxz_uu = 0;
    double Sxz_vv = 0;

    double S34_uu;
    double S34_vv;
    double S34_uuu;
    double S34_vvv;
    double S34_uvv;
    double S34_uuv;
    double S35_uu;
    double S35_vv;
    double S35_uuu;
    double S35_vvv;
    double S35_uvv;
    double S35_uuv;

    typedef Eigen::SparseMatrix<double, Eigen::RowMajor> Mat;

    double K_uv, K_uw;

    double theta;
    size_t reduced_vec_len;
    std::vector<double> F;
    std::vector<double> phi;
    std::vector<double> psi_shear;

    std::vector<double> permute_shear_3D;
    std::vector<double> rotate_X_Z_3D;
    std::vector<double> rotate_X_Z_3D_inv;

    std::vector<double> permute_and_rotate_3D;
    std::vector<double> permute_and_rotate_3D_inv;

    void calculate_torsion(const std::vector<std::unique_ptr<MeshNode>>& boundary_nodes, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh);
    void calculate_shear_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const std::vector<utils::LineBoundary>& line_bound);

    void base_matrix_upos(general_solver::MUMPSGeneral& M, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const size_t num_nodes, const size_t F_offset) const;
    void base_matrix_id(general_solver::MUMPSGeneral& M, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const size_t num_nodes, const size_t F_offset) const;

    double integrate_surface_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const std::function<double(const MeshElement* const, const gp_Pnt&, const gp_Pnt& px)>& fn) const;
    double GS_tri(const MeshElement* const e, const std::array<gp_Pnt, 3>& p, const std::array<gp_Pnt, 3>& px, const std::function<double(const MeshElement* const, const gp_Pnt&, const gp_Pnt& px)>& fn) const;
    double GS_quad(const MeshElement* const e, const std::array<gp_Pnt, 3>& p, const std::array<gp_Pnt, 3>& px, const std::function<double(const MeshElement* const, const gp_Pnt&, const gp_Pnt& px)>& fn) const;

    Eigen::Matrix<double, 6, 6> get_S_3D(const MeshElement* const e, const gp_Pnt& p) const;
    Eigen::Matrix<double, 6, 6> get_B_3D(const MeshElement* const e, const gp_Pnt& p) const;

    void get_B_tensors_3D(const MeshElement* const e, const gp_Pnt& p, Eigen::Matrix<double, 3, 3>& B4, Eigen::Matrix<double, 3, 2>& B3, Eigen::Matrix<double, 2, 2>& B2) const;

    inline double make_EA_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        (void)px;
        const auto S = this->get_S_3D(e, p);
        return 1/S(2,2);
    }
    inline double make_EA_u_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        return px.X()/S(2,2);
    }
    inline double make_EA_v_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        return px.Y()/S(2,2);
    }
    inline double make_EA_w_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        return px.Z()/S(2,2);
    }
    inline double make_EI_uu_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dy = px.Y() - c_v;
        return dy*dy/S(2,2);
    }
    inline double make_EI_vv_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dx = px.X() - c_u;
        return dx*dx/S(2,2);
    }
    inline double make_EI_uv_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dx = px.X() - c_u;
        const double dy = px.Y() - c_v;
        return dx*dy/S(2,2);
    }
    inline double make_EI_uuu_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dy = px.Y() - c_v;
        return dy*dy*dy/S(2,2);
    }
    inline double make_EI_vvv_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dx = px.X() - c_u;
        return dx*dx*dx/S(2,2);
    }
    inline double make_EI_uvv_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dx = px.X() - c_u;
        const double dy = px.Y() - c_v;
        return dx*dx*dy/S(2,2);
    }
    inline double make_EI_uuv_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dx = px.X() - c_u;
        const double dy = px.Y() - c_v;
        return dx*dy*dy/S(2,2);
    }
    inline double make_S34_uu_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dy = px.Y() - c_v;
        return S(2,3)*dy*dy/(S(2,2)*S(2,2));
    }
    inline double make_S34_vv_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dx = px.X() - c_u;
        return S(2,3)*dx*dx/(S(2,2)*S(2,2));
    }
    inline double make_S35_uu_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dy = px.Y() - c_v;
        return S(2,4)*dy*dy/(S(2,2)*S(2,2));
    }
    inline double make_S35_vv_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dx = px.X() - c_u;
        return S(2,4)*dx*dx/(S(2,2)*S(2,2));
    }
    inline double make_S34_uuu_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dy = px.Y() - c_v;
        return S(2,3)*dy*dy*dy/(S(2,2)*S(2,2));
    }
    inline double make_S34_uuv_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dx = px.X() - c_u;
        const double dy = px.Y() - c_v;
        return S(2,3)*dx*dy*dy/(S(2,2)*S(2,2));
    }
    inline double make_S34_uvv_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dx = px.X() - c_u;
        const double dy = px.Y() - c_v;
        return S(2,3)*dx*dx*dy/(S(2,2)*S(2,2));
    }
    inline double make_S34_vvv_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dx = px.X() - c_u;
        return S(2,3)*dx*dx*dx/(S(2,2)*S(2,2));
    }
    inline double make_S35_uuu_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dy = px.Y() - c_v;
        return S(2,4)*dy*dy*dy/(S(2,2)*S(2,2));
    }
    inline double make_S35_uuv_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dx = px.X() - c_u;
        const double dy = px.Y() - c_v;
        return S(2,4)*dx*dy*dy/(S(2,2)*S(2,2));
    }
    inline double make_S35_uvv_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dx = px.X() - c_u;
        const double dy = px.Y() - c_v;
        return S(2,4)*dx*dx*dy/(S(2,2)*S(2,2));
    }
    inline double make_S35_vvv_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const auto S = this->get_S_3D(e, p);
        const double dx = px.X() - c_u;
        return S(2,4)*dx*dx*dx/(S(2,2)*S(2,2));
    }

    inline void get_L_xx(const gp_Pnt& p0, const gp_Pnt& p1, double& L_uu, double& L_u, double& L_v, double& L_vv) const{
        const double x0 = p0.X();
        const double y0 = p0.Y();
        const double x1 = p1.X();
        const double y1 = p1.Y();
        const double d = p0.Distance(p1);

        L_uu = d*(x0*x0 + x0*x1 + x1*x1)/3.0;
        L_u  = d*(x0+x1)/2.0;
        L_v  = d*(y0+y1)/2.0;
        L_vv = d*(y0*y0 + y0*y1 + y1*y1)/3.0;
    }

    inline void get_S_xx(const gp_Pnt& p0, const gp_Pnt& p1, double& S_uudv, double& S_vvdv, double& S_vvdu, double& S_uudu) const{
        const double x0 = p0.X();
        const double y0 = p0.Y();
        const double x1 = p1.X();
        const double y1 = p1.Y();
        if(std::abs(x1 - x0) > 1e-14){
            const double a = (y1 - y0)/(x1 - x0);
           
            if(std::abs(a) > 1e-14){ 
                S_vvdu = (y1*y1*y1 - y0*y0*y0)/(3*a);
            } else {
                S_vvdu = y0*y0*(x1 - x0);
            }
            S_uudu = (x1*x1*x1 - x0*x0*x0)/3.0;
        } else {
            S_vvdu = 0;
            S_uudu = 0;
        }
        if(std::abs(y1 - y0) > 1e-14){
            const double c = (x1 - x0)/(y1 - y0);
            
            if(std::abs(c) > 1e-14){
                S_uudv = (x1*x1*x1 - x0*x0*x0)/(3*c);
            } else {
                S_uudv = x0*x0*(y1 - y0);
            }
            S_vvdv = (y1*y1*y1 - y0*y0*y0)/3.0;
        } else {
            S_uudv = 0;
            S_vvdv = 0;
        }
    }
};

#endif
