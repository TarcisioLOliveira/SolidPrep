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
#include "material.hpp"
#include "element.hpp"
#include "element_factory.hpp"
#include "utils/boundary_nullifier.hpp"

class Curvature{
    public:

    Curvature(const Material* mat, gp_Dir u, gp_Dir v, gp_Dir w, Eigen::Matrix<double, 2, 2> rot2D, Eigen::Matrix<double, 3, 3> rot3D, const BoundaryMeshElementFactory* elem_info, double V_v, double V_w, double M_u, double M_v, double M_w);

    void generate_curvature_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const std::vector<utils::LineBoundary>& line_bound, size_t phi_size, size_t psi_size);

    inline std::array<double, 2> get_curvatures() const{
        return {curv_v, curv_w};
    }
    inline std::array<double, 2> get_curvatures_dx() const{
        return {dcurv_v, dcurv_w};
    }
    inline gp_Pnt get_center() const{
        return gp_Pnt(0, c_v, c_w);
    }
    void get_shear_in_3D(const BoundaryMeshElement* e, double& t_uv, double& t_uw) const;

    private:
    const Material* mat;
    const gp_Dir u, v, w;
    const Eigen::Matrix<double, 2, 2> rot2D;
    const Eigen::Matrix<double, 3, 3> rot3D;
    const BoundaryMeshElementFactory* elem_info;
    const double V_v, V_w;
    const double M_u, M_v, M_w;
    const double line_mesh_size = 0.5;

    double EA;
    double EI_v;
    double EI_w;
    double EI_vw;
    double c_v, c_w;
    double curv_v;
    double curv_w;
    double dcurv_v;
    double dcurv_w;

    double theta;
    std::vector<double> phi_torsion;
    std::vector<double> psi_shear;
    utils::BoundaryNullifier<1, 2> bn1;
    utils::BoundaryNullifier<2, 1> bn2;

    void calculate_torsion(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh);
    void calculate_shear_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const std::vector<utils::LineBoundary>& line_bound);

    double integrate_surface_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const std::function<double(const MeshElement* const, const gp_Pnt&, const gp_Pnt& px)>& fn) const;
    double GS_tri(const MeshElement* const e, const std::array<gp_Pnt, 3>& p, const std::array<gp_Pnt, 3>& px, const std::function<double(const MeshElement* const, const gp_Pnt&, const gp_Pnt& px)>& fn) const;

    inline double make_EA_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        (void)px;
        return this->mat->beam_E_3D(e, p, this->rot3D);
    }
    inline double make_EA_v_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        return this->mat->beam_E_3D(e, p, this->rot3D)*px.Y();
    }
    inline double make_EA_w_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        return this->mat->beam_E_3D(e, p, this->rot3D)*px.Z();
    }
    inline double make_EI_v_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const double dz = px.Z() - c_w;
        return this->mat->beam_E_3D(e, p, this->rot3D)*dz*dz;
    }
    inline double make_EI_w_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const double dy = px.Y() - c_v;
        return this->mat->beam_E_3D(e, p, this->rot3D)*dy*dy;
    }
    inline double make_EI_vw_base_3D(const MeshElement* const e, const gp_Pnt& p, const gp_Pnt& px) const{
        const double dy = px.Y() - c_v;
        const double dz = px.Z() - c_w;
        return this->mat->beam_E_3D(e, p, this->rot3D)*dy*dz;
    }
};

#endif
