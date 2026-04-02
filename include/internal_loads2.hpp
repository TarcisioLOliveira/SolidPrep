/*
 *   Copyright (C) 2026 Tarcísio Ladeia de Oliveira.
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

#ifndef INTERNAL_LOADS_2_HPP
#define INTERNAL_LOADS_2_HPP

#include <array>
#include "material.hpp"
#include "math/matrix.hpp"
#include "utils.hpp"
#include "cross_section.hpp"
#include "element_factory.hpp"
#include "utils/delayed_pointer.hpp"
#include "utils/gauss_legendre.hpp"

class Meshing;

class InternalLoads2{
    public:

    // ASSUMES CROSS SECTION IS PLANE
    // ASSUMES this->normal == this->S.normal

    InternalLoads2(CrossSection cross_section, double thickness, gp_Dir normal, gp_Dir v, gp_Dir w, utils::DelayedPointerView<Material> mat, std::array<double, 3> F, std::array<double, 3> M, utils::ProblemType type);
    InternalLoads2(InternalLoads2&&) = default;

    void initialize(const Meshing* const mesh, size_t il_id, std::vector<double>& load);

    inline math::Matrix get_B_ext(const gp_Pnt& p) const{
        const math::Vector pv({p.X(), p.Y(), p.Z()});
        if(this->type == utils::PROBLEM_TYPE_3D){
            const math::Vector p2 = this->rot3D*pv - this->center;
            const double x = p2[0];
            const double y = p2[1];
            const double z = p2[2];

            const math::Matrix M({
             0, 0, 0, 0, 2*x*z, y*z,
             0, 0, 0, 0, x*z, 2*y*z,
             x, y, 1, 0, 0, 0,
             0, 0, 0, 0, y*z, x*z,
             0, 0, 0, -y, x*x, x*y,
             0, 0, 0,  x, x*y, y*y
            }, 6, 6);

            return rot3D_D.T()*M;
        } else {
            return math::Matrix();
        }
    }
    inline math::Matrix get_N_ext(const gp_Pnt& p) const{
        const math::Vector pv({p.X(), p.Y(), p.Z()});
        if(this->type == utils::PROBLEM_TYPE_3D){
            const math::Vector p2 = this->rot3D*pv - this->center;
            const double x = p2[0];
            const double y = p2[1];
            const double z = p2[2];

            const math::Matrix M({
             0, 0, 0, -y*z, x*x*z, x*y*z,
             0, 0, 0,  x*z, x*y*z, y*y*z,
             x*z, y*z, z, 0, 0, 0
            }, 3, 6);

            return rot3D.T()*M;
        } else {
            return math::Matrix();
        }
    }

    math::Vector calculate_internal_loads(const std::vector<double>& u_ext) const;
    math::Vector calculate_constraints(const std::vector<double>& u_ext, bool print = false) const;
    math::Vector calculate_multipliers(const std::vector<double>& u_ext, bool print = false) const;
    math::Vector get_loads() const;
    void update_multipliers(const math::Vector& params_diff, double alpha, bool print = false);

    math::Matrix get_constr(const BoundaryElement* b) const;

    math::Vector get_grad(const BoundaryElement* b, const math::Vector& mult) const;
    math::Matrix get_hess(const BoundaryElement* b) const;

    inline math::Vector get_center() const{
        return this->center;
    }

    const CrossSection S;
    const double A;
    const double thickness;
    const math::Matrix rot2D;
    const math::Matrix Lek_basis;
    const math::Matrix rot3D;
    const math::Matrix rot2D_D;
    const math::Matrix rot2D_S;
    const math::Matrix rot3D_D;
    const math::Matrix rot3D_S;
    const math::Vector F;
    const math::Vector M;
    utils::DelayedPointerView<Material> mat;
    const gp_Dir normal;
    std::vector<const BoundaryElement*> submesh;
    math::Vector lambda = math::Vector(6);
    const double MU;
    math::Vector params = math::Vector(6);
    math::Vector base_strain;

    private:
    const gp_Dir v;
    const gp_Dir w;
    const utils::ProblemType type;
    math::Vector center = math::Vector(3);
    const MeshElementFactory* elem_info;
    size_t BOUNDARY_NODE_NUM;
    

    math::Vector integrate_surface_3D(const Meshing* const mesh, const std::vector<std::function<double(const math::Matrix& S, const gp_Pnt& px)>>& fn) const;
    void GS_tri(const MeshElement* const e, const std::array<gp_Pnt, 3>& p, const std::vector<std::function<double(const math::Matrix& S, const gp_Pnt& px)>>& fn, math::Vector& result) const;
    void GS_quad(const MeshElement* const e, const std::array<gp_Pnt, 4>& p, const std::vector<std::function<double(const math::Matrix& S, const gp_Pnt& px)>>& fn, math::Vector& result) const;

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
