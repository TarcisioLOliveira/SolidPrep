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

#ifndef INTERNAL_LOADS_3_HPP
#define INTERNAL_LOADS_3_HPP

#include <array>
#include <unordered_map>
#include "general_solver/mumps_general.hpp"
#include "general_solver/petsc_general_pcg.hpp"
#include "material.hpp"
#include "math/matrix.hpp"
#include "multimaterial.hpp"
#include "utils.hpp"
#include "cross_section.hpp"
#include "element_factory.hpp"
#include "utils/delayed_pointer.hpp"
#include "utils/gauss_legendre.hpp"

class Meshing;

class InternalLoads3{
    public:

    // ASSUMES CROSS SECTION IS PLANE
    // ASSUMES this->normal == this->S.normal

    InternalLoads3(CrossSection cross_section, double thickness, gp_Dir normal, gp_Dir v, gp_Dir w, utils::DelayedPointerView<Material> mat, std::array<double, 3> F, std::array<double, 3> M, utils::ProblemType type);
    InternalLoads3(InternalLoads3&&) = default;

    void initialize(const Meshing* const mesh);
    void apply_load(const std::vector<long>& node_position, std::vector<double>& load) const;

    math::Vector calculate_internal_loads(const std::vector<double>& u_ext) const;
    math::Vector get_loads() const;

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
    std::vector<double> u;
    size_t u_size = 0;

    std::unordered_map<size_t, long> F_id_map;

    general_solver::MUMPSGeneral solver;
    //general_solver::PETScGeneralPCG solver = general_solver::PETScGeneralPCG(true);

    private:
    const gp_Dir v;
    const gp_Dir w;
    const utils::ProblemType type;
    math::Vector center = math::Vector(3);
    const MeshElementFactory* elem_info;
    size_t BOUNDARY_NODE_NUM;
    MultiMaterial mmat;
    
    void calculate_force_distribution();

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
