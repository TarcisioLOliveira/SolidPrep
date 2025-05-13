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
#include <gp_Dir.hxx>
#include "material.hpp"
#include "element.hpp"
#include "element_factory.hpp"
#include "math/matrix.hpp"
#include "utils/gauss_legendre.hpp"

class Curvature{
    public:

    Curvature(const Material* mat, math::Matrix rot2D, math::Matrix rot3D, const BoundaryMeshElementFactory* elem_info, double V_u, double V_v, double V_w, double M_u, double M_v, double M_w);

    void generate_curvature_3D(const std::vector<std::unique_ptr<MeshNode>>& boundary_nodes, const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, size_t reduced_vector_size, size_t full_vector_size);

    inline gp_Pnt get_center() const{
        return this->center;
    }
    void get_stress_3D(const BoundaryMeshElement* e, double& t_uw, double& t_vw, double& s_w) const;

    math::Vector get_force_vector_3D(const BoundaryMeshElement* e) const;

    private:
    const Material* mat;
    const math::Matrix rot2D;
    const math::Matrix Lek_basis;
    const math::Matrix rot3D;
    const BoundaryMeshElementFactory* elem_info;
    const double V_u, V_v, V_w;
    const double M_u, M_v, M_w;
    size_t reduced_vec_len;
    size_t full_vector_size;

    double Area   = 0;
    double EA     = 0;
    double EA_u   = 0;
    double EA_v   = 0;
    double EI_uu  = 0;
    double EI_vv  = 0;
    double EI_uv  = 0;
    double c_u, c_v, c_w;
    gp_Pnt center;
    double A, B, C;
    double theta;

    double Kyz = 0;
    double Kxz = 0;

    std::vector<double> u;

    math::Matrix permute_shear_3D;
    math::Matrix rotate_X_Z_3D;
    math::Matrix rotate_X_Z_3D_inv;

    math::Matrix permute_and_rotate_3D;
    math::Matrix permute_and_rotate_3D_inv;

    void calculate_stress_field_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh);

    math::Vector integrate_surface_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const std::vector<std::function<double(const math::Matrix& S, const gp_Pnt& px)>>& fn) const;
    void GS_tri(const MeshElement* const e, const std::array<gp_Pnt, 3>& p, const std::array<gp_Pnt, 3>& px, const std::vector<std::function<double(const math::Matrix& S, const gp_Pnt& px)>>& fn, math::Vector& result) const;
    void GS_quad(const MeshElement* const e, const std::array<gp_Pnt, 4>& p, const std::array<gp_Pnt, 4>& px, const std::vector<std::function<double(const math::Matrix& S, const gp_Pnt& px)>>& fn, math::Vector& result) const;

    math::Matrix get_S_3D(const MeshElement* const e, const gp_Pnt& p) const;
    math::Matrix get_D_3D(const MeshElement* const e, const gp_Pnt& p) const;

    math::Matrix get_T_D_3D(const math::Matrix& S, const math::Matrix& D) const;
    math::Matrix get_DT_3D(const math::Matrix& S, const math::Matrix& D) const;
    math::Matrix get_TDT_3D(const math::Matrix& S, const math::Matrix& D) const;

    void get_B_tensors_3D(const MeshElement* const e, const gp_Pnt& p, math::Matrix& B4, math::Matrix& B3, math::Matrix& B2) const;

    inline double make_A_base_3D(const math::Matrix& S, const gp_Pnt& px) const{
        (void)S;
        (void)px;
        return 1;
    }
    inline double make_EA_base_3D(const math::Matrix& S, const gp_Pnt& px) const{
        (void)px;
        return 1.0/S(2,2);
    }
    inline double make_EA_u_base_3D(const math::Matrix& S, const gp_Pnt& px) const{
        return px.X()/S(2,2);
    }
    inline double make_EA_v_base_3D(const math::Matrix& S, const gp_Pnt& px) const{
        return px.Y()/S(2,2);
    }
    inline double make_EA_w_base_3D(const math::Matrix& S, const gp_Pnt& px) const{
        return px.Z()/S(2,2);
    }
    inline double make_EI_uu_base_3D(const math::Matrix& S, const gp_Pnt& px) const{
        const double dy = px.Y() - c_v;
        return dy*dy/S(2,2);
    }
    inline double make_EI_vv_base_3D(const math::Matrix& S, const gp_Pnt& px) const{
        const double dx = px.X() - c_u;
        return dx*dx/S(2,2);
    }
    inline double make_EI_uv_base_3D(const math::Matrix& S, const gp_Pnt& px) const{
        const double dx = px.X() - c_u;
        const double dy = px.Y() - c_v;
        return dx*dy/S(2,2);
    }

    inline gp_Pnt GLT_point(const utils::GLPointTri& glp, const std::array<gp_Pnt, 3>& px) const{
            return gp_Pnt{
                glp.a*px[0].X() + glp.b*px[1].X() + glp.c*px[2].X(),
                glp.a*px[0].Y() + glp.b*px[1].Y() + glp.c*px[2].Y(),
                glp.a*px[0].Z() + glp.b*px[1].Z() + glp.c*px[2].Z()
            };
    }

    inline double N_norm(double x, double y, size_t i) const{
        switch(i){
            case 0:
                return 0.25*(1-x)*(1-y);
            case 1:
                return 0.25*(1+x)*(1-y);
            case 2:
                return 0.25*(1+x)*(1+y);
            case 3:
                return 0.25*(1-x)*(1+y);
        }
        return 0;
    }
    inline math::Matrix J(const double xi, const double eta, const std::array<gp_Pnt, 4>& px) const{
        std::array<double, 4> x{
            px[0].X(),
            px[1].X(),
            px[2].X(),
            px[3].X()
        };
        std::array<double, 4> y{
            px[0].Y(),
            px[1].Y(),
            px[2].Y(),
            px[3].Y()
        };
        return math::Matrix(
        {
            eta*x[0]/4 - eta*x[1]/4 + eta*x[2]/4 - eta*x[3]/4 - x[0]/4 + x[1]/4 + x[2]/4 - x[3]/4
            ,
            eta*y[0]/4 - eta*y[1]/4 + eta*y[2]/4 - eta*y[3]/4 - y[0]/4 + y[1]/4 + y[2]/4 - y[3]/4
            ,
            x[0]*xi/4 - x[0]/4 - x[1]*xi/4 - x[1]/4 + x[2]*xi/4 + x[2]/4 - x[3]*xi/4 + x[3]/4
            ,
            xi*y[0]/4 - xi*y[1]/4 + xi*y[2]/4 - xi*y[3]/4 - y[0]/4 - y[1]/4 + y[2]/4 + y[3]/4
        }, 2, 2);
    }

    inline gp_Pnt GL_point(const utils::GLPoint& xi, const utils::GLPoint& eta, const std::array<gp_Pnt, 4>& px) const{
        double X = 0, Y = 0, Z = 0;
        for(size_t i = 0; i < px.size(); ++i){
            const double Ni = N_norm(xi.x, eta.x, i);
            const gp_Pnt& pi = px[i];
            X += Ni*pi.X();
            Y += Ni*pi.Y();
            Z += Ni*pi.Z();
        }

        return gp_Pnt(X, Y, Z);
    }
};

#endif
