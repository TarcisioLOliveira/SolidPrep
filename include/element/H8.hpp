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

#ifndef H8_HPP
#define H8_HPP

#include "element.hpp"
#include "material.hpp"
#include "element_factory.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <memory>
#include "element_common.hpp"
#include "math/matrix.hpp"
#include <array>

namespace element{

class H8 : public MeshElementCommon3DHex<H8>{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 5;
    static const size_t NODE_DOF       = 3;
    static const size_t NODES_PER_ELEM = 8;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;

    static const size_t BOUNDARY_NODES_PER_ELEM = 4;
    static const size_t BOUNDARY_GMSH_TYPE = 3;
    static std::unique_ptr<BoundaryMeshElementFactory> get_boundary_element_info();
    static std::unique_ptr<ContactMeshElementFactory> get_contact_element_info();

    static const spview::defs::ElementType SPVIEW_CODE = spview::defs::H8;

    enum CubeSide{
        X_MIN,
        X_MAX,
        Y_MIN,
        Y_MAX,
        Z_MIN,
        Z_MAX,
        UNKNOWN
    };

    CubeSide get_cube_side(const std::vector<gp_Pnt>& points) const;
    gp_Pnt to_surface_point(const double xi, const double eta, const CubeSide side) const;

    H8(ElementShape s);

    virtual math::Matrix get_k(const math::Matrix& D, const double t) const override;
    virtual math::Matrix get_nodal_density_gradient(gp_Pnt p) const override;
    virtual math::Matrix get_R(const math::Matrix& K, const double t, const std::vector<gp_Pnt>& points) const override;
    virtual math::Matrix get_B(const gp_Pnt& point) const override;

    virtual math::Matrix diffusion_1dof(const double t, const math::Matrix& A) const override;
    virtual math::Matrix advection_1dof(const double t, const math::Vector& v) const override;
    virtual math::Matrix absorption_1dof(const double t) const override;
    virtual math::Vector source_1dof(const double t) const override;
    virtual math::Vector flow_1dof(const double t, const MeshNode** nodes) const override;

    virtual math::Matrix get_Ni(const gp_Pnt& p) const override;

    virtual math::Matrix get_uu(const MeshElement* const e2, const std::vector<gp_Pnt>& bounds, const gp_Dir n) const override;
    virtual math::Matrix get_MnMn(const MeshElement* const e2, const std::vector<double>& u_ext, const std::vector<gp_Pnt>& bounds, const gp_Dir n) const override;

    virtual inline std::unique_ptr<MeshElementFactory> get_element_info() const override{
        return std::unique_ptr<MeshElementFactory>(new MeshElementFactoryImpl<H8>());
    }

    private:
    static const bool reg;
    virtual math::Matrix get_Nf(const double t, const std::vector<gp_Pnt>& points) const override;

    inline double N_norm(double x, double y, double z, size_t i) const{
        switch(i){
            case 0:
                return 0.125*(1-x)*(1-y)*(1-z);
            case 1:
                return 0.125*(1+x)*(1-y)*(1-z);
            case 2:
                return 0.125*(1+x)*(1+y)*(1-z);
            case 3:
                return 0.125*(1-x)*(1+y)*(1-z);
            case 4:
                return 0.125*(1-x)*(1-y)*(1+z);
            case 5:
                return 0.125*(1+x)*(1-y)*(1+z);
            case 6:
                return 0.125*(1+x)*(1+y)*(1+z);
            case 7:
                return 0.125*(1-x)*(1+y)*(1+z);
        }
        return 0;
    }
    inline double dNdx_norm(double x, double y, double z, size_t i) const{
        (void)x;
        switch(i){
            case 0:
                return -0.125*(1-y)*(1-z);
            case 1:
                return  0.125*(1-y)*(1-z);
            case 2:
                return  0.125*(1+y)*(1-z);
            case 3:
                return -0.125*(1+y)*(1-z);
            case 4:
                return -0.125*(1-y)*(1+z);
            case 5:
                return  0.125*(1-y)*(1+z);
            case 6:
                return  0.125*(1+y)*(1+z);
            case 7:
                return -0.125*(1+y)*(1+z);
        }
        return 0;
    }
    inline double dNdy_norm(double x, double y, double z, size_t i) const{
        (void)y;
        switch(i){
            case 0:
                return -0.125*(1-x)*(1-z);
            case 1:
                return -0.125*(1+x)*(1-z);
            case 2:
                return  0.125*(1+x)*(1-z);
            case 3:
                return  0.125*(1-x)*(1-z);
            case 4:
                return -0.125*(1-x)*(1+z);
            case 5:
                return -0.125*(1+x)*(1+z);
            case 6:
                return  0.125*(1+x)*(1+z);
            case 7:
                return  0.125*(1-x)*(1+z);
        }
        return 0;
    }
    inline double dNdz_norm(double x, double y, double z, size_t i) const{
        (void)z;
        switch(i){
            case 0:
                return -0.125*(1-x)*(1-y);
            case 1:
                return -0.125*(1+x)*(1-y);
            case 2:
                return -0.125*(1+x)*(1+y);
            case 3:
                return -0.125*(1-x)*(1+y);
            case 4:
                return  0.125*(1-x)*(1-y);
            case 5:
                return  0.125*(1+x)*(1-y);
            case 6:
                return  0.125*(1+x)*(1+y);
            case 7:
                return  0.125*(1-x)*(1+y);
        }
        return 0;
    }

    inline math::Vector N_mat_1dof(double x, double y, double z) const{
        return math::Vector
                {N_norm(x,y,z,0), N_norm(x,y,z,1), N_norm(x,y,z,2), N_norm(x,y,z,3),
                 N_norm(x,y,z,4), N_norm(x,y,z,5), N_norm(x,y,z,6), N_norm(x,y,z,7)};
    }

    inline math::Matrix dN_mat_1dof(double xi, double eta, double zeta,
                                                                       const std::array<double, NODES_PER_ELEM>& x, 
                                                                       const std::array<double, NODES_PER_ELEM>& y,
                                                                       const std::array<double, NODES_PER_ELEM>& z) const{
        math::Matrix Jm = J(xi, eta, zeta, x, y, z);
        Jm.invert_LU();
        math::Matrix M(
                {dNdx_norm(xi,eta,zeta,0),dNdx_norm(xi,eta,zeta,1),dNdx_norm(xi,eta,zeta,2),dNdx_norm(xi,eta,zeta,3),dNdx_norm(xi,eta,zeta,4),dNdx_norm(xi,eta,zeta,5),dNdx_norm(xi,eta,zeta,6),dNdx_norm(xi,eta,zeta,7),
                 dNdy_norm(xi,eta,zeta,0),dNdy_norm(xi,eta,zeta,1),dNdy_norm(xi,eta,zeta,2),dNdy_norm(xi,eta,zeta,3),dNdy_norm(xi,eta,zeta,4),dNdy_norm(xi,eta,zeta,5),dNdy_norm(xi,eta,zeta,6),dNdy_norm(xi,eta,zeta,7),
                 dNdz_norm(xi,eta,zeta,0),dNdz_norm(xi,eta,zeta,1),dNdz_norm(xi,eta,zeta,2),dNdz_norm(xi,eta,zeta,3),dNdz_norm(xi,eta,zeta,4),dNdz_norm(xi,eta,zeta,5),dNdz_norm(xi,eta,zeta,6),dNdz_norm(xi,eta,zeta,7)},
                NODE_DOF, NODES_PER_ELEM);

        return Jm*M;
    }

    inline math::Matrix N_mat_norm(double x, double y, double z) const{
        double Ni[NODES_PER_ELEM] = 
                {N_norm(x,y,z,0), N_norm(x,y,z,1), N_norm(x,y,z,2), N_norm(x,y,z,3),
                 N_norm(x,y,z,4), N_norm(x,y,z,5), N_norm(x,y,z,6), N_norm(x,y,z,7)};

        math::Matrix M(NODE_DOF, K_DIM);
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            M(0, 3*i  ) = Ni[i];
            M(1, 3*i+1) = Ni[i];
            M(2, 3*i+2) = Ni[i];
        }

        return M;
    }

    inline math::Vector norm_to_nat(double xi, double eta, double zeta, const std::array<double, NODES_PER_ELEM>& x, const std::array<double, NODES_PER_ELEM>& y, const std::array<double, NODES_PER_ELEM>& z) const{
        double X = 0, Y = 0, Z = 0;
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            const double Ni = N_norm(xi, eta, zeta, i);
            X += Ni*x[i];
            Y += Ni*y[i];
            Z += Ni*z[i];
        }

        return math::Vector{X, Y, Z};
    }


    inline math::Matrix J(const double xi, const double eta, const double zeta,
                                       const std::array<double, NODES_PER_ELEM>& x, 
                                       const std::array<double, NODES_PER_ELEM>& y,
                                       const std::array<double, NODES_PER_ELEM>& z) const{

        math::Matrix Jm(NODE_DOF, NODE_DOF);
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            const double Nix = dNdx_norm(xi,eta,zeta,i);
            const double Niy = dNdy_norm(xi,eta,zeta,i);
            const double Niz = dNdz_norm(xi,eta,zeta,i);

            Jm(0,0) += Nix*x[i];
            Jm(0,1) += Nix*y[i];
            Jm(0,2) += Nix*z[i];

            Jm(1,0) += Niy*x[i];
            Jm(1,1) += Niy*y[i];
            Jm(1,2) += Niy*z[i];

            Jm(2,0) += Niz*x[i];
            Jm(2,1) += Niz*y[i];
            Jm(2,2) += Niz*z[i];
        }

        return Jm;
    }

    inline math::Matrix B_mat_norm(double xi, double eta, double zeta, const std::array<double, NODES_PER_ELEM>& x, const std::array<double, NODES_PER_ELEM>& y, const std::array<double, NODES_PER_ELEM>& z) const{
        auto Ji = J(xi, eta, zeta, x, y, z);
        Ji.invert_LU();
        double Nix[NODES_PER_ELEM];
        double Niy[NODES_PER_ELEM];
        double Niz[NODES_PER_ELEM];
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            Nix[i] = dNdx_norm(xi,eta,zeta,i);
            Niy[i] = dNdy_norm(xi,eta,zeta,i);
            Niz[i] = dNdz_norm(xi,eta,zeta,i);
        }

        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            const double Nixi = Ji(0,0)*Nix[i] + Ji(0,1)*Niy[i] + Ji(0,2)*Niz[i];
            const double Niyi = Ji(1,0)*Nix[i] + Ji(1,1)*Niy[i] + Ji(1,2)*Niz[i];
            const double Nizi = Ji(2,0)*Nix[i] + Ji(2,1)*Niy[i] + Ji(2,2)*Niz[i];
            Nix[i] = Nixi;
            Niy[i] = Niyi;
            Niz[i] = Nizi;
        }

        math::Matrix Bm(S_SIZE, K_DIM);

        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            Bm(0, 3*i + 0) = Nix[i];
            Bm(1, 3*i + 1) = Niy[i];
            Bm(2, 3*i + 2) = Niz[i];
            Bm(3, 3*i + 0) = Niy[i];
            Bm(3, 3*i + 1) = Nix[i];
            Bm(4, 3*i + 0) = Niz[i];
            Bm(4, 3*i + 2) = Nix[i];
            Bm(5, 3*i + 1) = Niz[i];
            Bm(5, 3*i + 2) = Niy[i];
        }

        return Bm;
    }


    inline double N_norm_surface(double x, double y, size_t i) const{
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
    inline double dNdx_norm_surface(double x, double y, size_t i) const{
        (void)x;
        switch(i){
            case 0:
                return -0.25*(1-y);
            case 1:
                return  0.25*(1-y);
            case 2:
                return  0.25*(1+y);
            case 3:
                return -0.25*(1+y);
        }
        return 0;
    }
    inline double dNdy_norm_surface(double x, double y, size_t i) const{
        (void)y;
        switch(i){
            case 0:
                return -0.25*(1-x);
            case 1:               
                return -0.25*(1+x);
            case 2:               
                return  0.25*(1+x);
            case 3:               
                return  0.25*(1-x);
        }
        return 0;
    }
    inline math::Vector surface_to_nat(double xi, double eta, const std::array<double, BOUNDARY_NODES_PER_ELEM>& x, const std::array<double, BOUNDARY_NODES_PER_ELEM>& y, const std::array<double, BOUNDARY_NODES_PER_ELEM>& z) const{
        double X = 0, Y = 0, Z = 0;
        for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
            const double Ni = N_norm_surface(xi, eta, i);
            X += Ni*x[i];
            Y += Ni*y[i];
            Z += Ni*z[i];
        }

        return math::Vector{X, Y, Z};
    }
    inline double surface_drnorm(double xi, double eta, const std::array<double, BOUNDARY_NODES_PER_ELEM>& x, const std::array<double, BOUNDARY_NODES_PER_ELEM>& y, const std::array<double, BOUNDARY_NODES_PER_ELEM>& z) const{
        double dXdx = 0, dYdx = 0, dZdx = 0;
        double dXdy = 0, dYdy = 0, dZdy = 0;
        for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
            const double Nix = dNdx_norm_surface(xi, eta, i);
            const double Niy = dNdy_norm_surface(xi, eta, i);

            dXdx += Nix*x[i];
            dYdx += Nix*y[i];
            dZdx += Nix*z[i];

            dXdy += Niy*x[i];
            dYdy += Niy*y[i];
            dZdy += Niy*z[i];
        }
        gp_Vec rx{dXdx, dYdx, dZdx};
        gp_Vec ry{dXdy, dYdy, dZdy};

        return (rx.Crossed(ry)).Magnitude();
    }
};

}

#endif
