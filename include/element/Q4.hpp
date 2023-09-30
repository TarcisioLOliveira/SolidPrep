/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
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

#ifndef Q4_HPP
#define Q4_HPP

#include "element.hpp"
#include "material.hpp"
#include "utils.hpp"
#include "element_factory.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <memory>
#include "element_common.hpp"
#include <array>

namespace element{

class Q4 : public MeshElementCommon2DQuad<Q4>{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 3;
    static const size_t NODE_DOF       = 2;
    static const size_t NODES_PER_ELEM = 4;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;

    static const size_t BOUNDARY_NODES_PER_ELEM = 2;
    static const size_t BOUNDARY_GMSH_TYPE = 1;

    static const spview::defs::ElementType SPVIEW_CODE = spview::defs::Q4;

    Q4(ElementShape s);

    virtual std::vector<double> get_k(const std::vector<double>& D, const double t) const override;
    virtual std::vector<double> get_nodal_density_gradient(gp_Pnt p) const override;
    virtual std::vector<double> get_R(const std::vector<double>& K, const double t, const std::vector<gp_Pnt>& points) const override;
    virtual std::vector<double> get_Rf(const std::vector<double>& S, const std::vector<double>& F, const gp_Pnt& C, const double t, const std::vector<gp_Pnt>& points) const override;

    virtual Eigen::MatrixXd diffusion_1dof(const double t, const std::vector<double>& A) const override;
    virtual Eigen::MatrixXd advection_1dof(const double t, const std::vector<double>& v) const override;
    virtual Eigen::MatrixXd absorption_1dof(const double t) const override;
    virtual Eigen::VectorXd source_1dof(const double t) const override;

    virtual inline std::unique_ptr<MeshElementFactory> get_element_info() const override{
        return std::unique_ptr<MeshElementFactory>(new MeshElementFactoryImpl<Q4>());
    }

    private:
    virtual std::vector<double> get_DB(const std::vector<double>& D, const gp_Pnt& point) const override;
    virtual std::vector<double> get_Nf(const double t, const std::vector<gp_Pnt>& points) const override;
    virtual std::vector<double> get_B(const gp_Pnt& point) const override;

    double a[NODES_PER_ELEM], b[NODES_PER_ELEM], c[NODES_PER_ELEM], d[NODES_PER_ELEM], A;

    inline double N(double x, double y, size_t i) const{
        return a[i] + b[i]*x + c[i]*y + d[i]*x*y;
    }

    inline Eigen::Vector<double, NODES_PER_ELEM> N_mat_1dof(double x, double y) const{
        return Eigen::Vector<double, NODES_PER_ELEM>
                {N(x,y,0), N(x,y,1), N(x,y,2), N(x,y,3)};
    }

    inline double dNdx(double x, double y, size_t i) const{
        (void)x;
        return b[i] + d[i]*y;
    }
    inline double dNdy(double x, double y, size_t i) const{
        (void)y;
        return c[i] + d[i]*x;
    }
    inline Eigen::Matrix<double, NODE_DOF, NODES_PER_ELEM> dN_mat_1dof(double x, double y) const{
        return Eigen::Matrix<double, NODE_DOF, NODES_PER_ELEM>
                {{dNdx(x,y,0),dNdx(x,y,1),dNdx(x,y,2),dNdx(x,y,3)},
                 {dNdy(x,y,0),dNdy(x,y,1),dNdy(x,y,2),dNdy(x,y,3)}};
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
    inline double dNdx_norm(double x, double y, size_t i) const{
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
    inline double dNdy_norm(double x, double y, size_t i) const{
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

    
    inline Eigen::Matrix<double, NODE_DOF, K_DIM> N_mat(double x, double y) const{
        double Ni[4] = {N(x,y,0), N(x,y,1), N(x,y,2), N(x,y,3)};

        return Eigen::Matrix<double, NODE_DOF, K_DIM>
                {{Ni[0], 0, Ni[1], 0, Ni[2], 0, Ni[3], 0},
                 {0, Ni[0], 0, Ni[1], 0, Ni[2], 0, Ni[3]}};
    }

    inline Eigen::Matrix<double, NODE_DOF, K_DIM> N_mat_norm(double x, double y) const{
        double Ni[4] = {N_norm(x,y,0), N_norm(x,y,1), N_norm(x,y,2), N_norm(x,y,3)};

        return Eigen::Matrix<double, NODE_DOF, K_DIM>
                {{Ni[0], 0, Ni[1], 0, Ni[2], 0, Ni[3], 0},
                 {0, Ni[0], 0, Ni[1], 0, Ni[2], 0, Ni[3]}};
    }

    inline Eigen::Vector<double, NODE_DOF> norm_to_nat(double xi, double eta, const std::array<double, NODES_PER_ELEM>& x, const std::array<double, NODES_PER_ELEM>& y) const{
        double X = 0, Y = 0;
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            const double Ni = N_norm(xi, eta, i);
            X += Ni*x[i];
            Y += Ni*y[i];
        }

        return Eigen::Vector<double, NODE_DOF>{X, Y};
    }


    inline Eigen::Matrix<double, NODE_DOF, NODE_DOF> J(const double xi, const double eta, 
                                       const std::array<double, NODES_PER_ELEM>& x, const std::array<double, NODES_PER_ELEM>& y) const{

        return Eigen::Matrix<double, NODE_DOF, NODE_DOF>
        {{
            eta*x[0]/4 - eta*x[1]/4 + eta*x[2]/4 - eta*x[3]/4 - x[0]/4 + x[1]/4 + x[2]/4 - x[3]/4
            ,
            eta*y[0]/4 - eta*y[1]/4 + eta*y[2]/4 - eta*y[3]/4 - y[0]/4 + y[1]/4 + y[2]/4 - y[3]/4
            },{
            x[0]*xi/4 - x[0]/4 - x[1]*xi/4 - x[1]/4 + x[2]*xi/4 + x[2]/4 - x[3]*xi/4 + x[3]/4
            ,
            xi*y[0]/4 - xi*y[1]/4 + xi*y[2]/4 - xi*y[3]/4 - y[0]/4 - y[1]/4 + y[2]/4 + y[3]/4
        }};
    }
    inline Eigen::Matrix<double, S_SIZE, K_DIM> B_mat_nat(double X, double Y) const{
        const double Nix[NODES_PER_ELEM] = {dNdx(X,Y,0), dNdx(X,Y,1), dNdx(X,Y,2), dNdx(X,Y,3)};
        const double Niy[NODES_PER_ELEM] = {dNdy(X,Y,0), dNdy(X,Y,1), dNdy(X,Y,2), dNdy(X,Y,3)};

        return Eigen::Matrix<double, S_SIZE, K_DIM>
                {{Nix[0], 0, Nix[1], 0, Nix[2], 0, Nix[3], 0},
                 {0, Niy[0], 0, Niy[1], 0, Niy[2], 0, Niy[3]},
                 {Niy[0], Nix[0], Niy[1], Nix[1], Niy[2], Nix[2], Niy[3], Nix[3]}};
    }

    inline Eigen::Matrix<double, S_SIZE, K_DIM> B_mat_norm(double xi, double eta, const std::array<double, NODES_PER_ELEM>& x, const std::array<double, NODES_PER_ELEM>& y) const{
        auto Ji = J(xi, eta, x, y);
        Ji = Ji.inverse();
        double Nix[NODES_PER_ELEM] = {dNdx_norm(xi,eta,0), dNdx_norm(xi,eta,1), dNdx_norm(xi,eta,2), dNdx_norm(xi,eta,3)};
        double Niy[NODES_PER_ELEM] = {dNdy_norm(xi,eta,0), dNdy_norm(xi,eta,1), dNdy_norm(xi,eta,2), dNdy_norm(xi,eta,3)};

        for(size_t i = 0; i < NODE_DOF*NODE_DOF; ++i){
            double Nixi = Ji(0,0)*Nix[i] + Ji(0,1)*Niy[i];
            double Niyi = Ji(1,0)*Nix[i] + Ji(1,1)*Niy[i];
            Nix[i] = Nixi;
            Niy[i] = Niyi;
        }

        return Eigen::Matrix<double, S_SIZE, K_DIM>
                {{Nix[0], 0, Nix[1], 0, Nix[2], 0, Nix[3], 0},
                 {0, Niy[0], 0, Niy[1], 0, Niy[2], 0, Niy[3]},
                 {Niy[0], Nix[0], Niy[1], Nix[1], Niy[2], Nix[2], Niy[3], Nix[3]}};
    }
};

}

#endif
