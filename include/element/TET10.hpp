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

#ifndef TET10_HPP
#define TET10_HPP

#include <Eigen/Core>
#include "element.hpp"
#include "material.hpp"
#include <vector>
#include "utils.hpp"
#include "element_factory.hpp"
#include "element_common.hpp"
#include <array>

namespace element{

class TET10 : public MeshElementCommon3DTet<TET10>{
    public:
    static const size_t ORDER          = 2;
    static const size_t GMSH_TYPE      = 11;
    static const size_t NODE_DOF       = 3;
    static const size_t NODES_PER_ELEM = 10;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;

    static const size_t BOUNDARY_NODES_PER_ELEM = 6;
    static const size_t BOUNDARY_GMSH_TYPE = 9;
    static std::unique_ptr<BoundaryMeshElementFactory> get_boundary_element_info();

    static const spview::defs::ElementType SPVIEW_CODE = spview::defs::TET10;

    TET10(ElementShape s);

    virtual std::vector<double> get_k(const std::vector<double>& D, const double t) const override;
    virtual std::vector<double> get_nodal_density_gradient(gp_Pnt p) const override;
    virtual std::vector<double> get_R(const std::vector<double>& K, const double t, const std::vector<gp_Pnt>& points) const override;
    virtual std::vector<double> get_Rf(const std::vector<double>& S, const std::vector<double>& F, const gp_Pnt& C, const double t, const std::vector<gp_Pnt>& points) const override;
    virtual std::vector<double> get_B(const gp_Pnt& point) const override;

    virtual Eigen::MatrixXd diffusion_1dof(const double t, const std::vector<double>& A) const override;
    virtual Eigen::MatrixXd advection_1dof(const double t, const std::vector<double>& v) const override;
    virtual Eigen::MatrixXd absorption_1dof(const double t) const override;
    virtual Eigen::VectorXd source_1dof(const double t) const override;
    virtual Eigen::VectorXd flow_1dof(const double t, const MeshNode** nodes) const override;

    virtual inline std::unique_ptr<MeshElementFactory> get_element_info() const override{
        return std::unique_ptr<MeshElementFactory>(new MeshElementFactoryImpl<TET10>());
    }

    private:
    virtual std::vector<double> get_DB(const std::vector<double>& D, const gp_Pnt& point) const override;
    virtual std::vector<double> get_Nf(const double t, const std::vector<gp_Pnt>& points) const override;

    inline gp_Pnt GL_point(double c1, double c2, double c3, double c4) const{
        return gp_Pnt(
            c1*this->nodes[0]->point.X() + c2*this->nodes[1]->point.X() + c3*this->nodes[2]->point.X() + c4*this->nodes[3]->point.X(),
            c1*this->nodes[0]->point.Y() + c2*this->nodes[1]->point.Y() + c3*this->nodes[2]->point.Y() + c4*this->nodes[3]->point.Y(),
            c1*this->nodes[0]->point.Z() + c2*this->nodes[1]->point.Z() + c3*this->nodes[2]->point.Z() + c4*this->nodes[3]->point.Z()
        );
    }

    inline gp_Pnt GL_point_2D(double c1, double c2, double c3, const std::vector<gp_Pnt>& p) const{
        return gp_Pnt(
            c1*p[0].X() + c2*p[1].X() + c3*p[2].X(),
            c1*p[0].Y() + c2*p[1].Y() + c3*p[2].Y(),
            c1*p[0].Z() + c2*p[1].Z() + c3*p[2].Z()
        );
    }

    void get_coeffs();

    double a[4], b[4], c[4], d[4], V;
    static constexpr double n = 2;

    inline double L(const gp_Pnt& p, size_t i) const{
        return (a[i] + b[i]*p.X() + c[i]*p.Y() + d[i]*p.Z())/(6*V);
    }
    inline std::array<double, 4> La(const gp_Pnt& p) const{
        return {L(p,0), L(p,1), L(p,2), L(p,3)};
    }
    inline double Na(double Li, size_t a) const{
        switch(a){
            case 0:
                return 1;
            case 1:
                return 2*Li;
            case 2:
                return Li*(2*Li-1);
        }
        return 1;
    }
    inline double dNadx(double Li, size_t i, size_t a) const{
        switch(a){
            case 0:
                return 1;
            case 1:
                return 2*b[i]/(6*V);
            case 2:
                return (4*Li*b[i] - b[i])/(6*V);
        }
        return 1;
    }
    inline double dNady(double Li, size_t i, size_t a) const{
        switch(a){
            case 0:
                return 1;
            case 1:
                return 2*c[i]/(6*V);
            case 2:
                return (4*Li*c[i] - c[i])/(6*V);
        }
        return 1;
    }
    inline double dNadz(double Li, size_t i, size_t a) const{
        switch(a){
            case 0:
                return 1;
            case 1:
                return 2*d[i]/(6*V);
            case 2:
                return (4*Li*d[i] - d[i])/(6*V);
        }
        return 1;
    }
    inline double N(const std::array<double, 4>& Li, size_t i) const{
        switch(i){
            case 0:
                return Na(Li[0], 2);
            case 1:
                return Na(Li[1], 2);
            case 2:
                return Na(Li[2], 2);
            case 3:
                return Na(Li[3], 2);
            case 4:
                return Na(Li[0], 1)*Na(Li[1], 1);
            case 5:
                return Na(Li[1], 1)*Na(Li[2], 1);
            case 6:
                return Na(Li[0], 1)*Na(Li[2], 1);
            case 7:
                return Na(Li[0], 1)*Na(Li[3], 1);
            case 8:
                return Na(Li[2], 1)*Na(Li[3], 1);
            case 9:
                return Na(Li[1], 1)*Na(Li[3], 1);
        }
        return 0;
    }
    inline double dNdx(const std::array<double, 4>& Li, size_t i) const{
        switch(i){
            case 0:
                return dNadx(Li[0], 0, 2);
            case 1:                  
                return dNadx(Li[1], 1, 2);
            case 2:                  
                return dNadx(Li[2], 2, 2);
            case 3:                  
                return dNadx(Li[3], 3, 2);
            case 4:                  
                return dNadx(Li[0], 0, 1)*Na(Li[1], 1) + Na(Li[0], 1)*dNadx(Li[1], 1, 1);
            case 5:                                                                 
                return dNadx(Li[1], 1, 1)*Na(Li[2], 1) + Na(Li[1], 1)*dNadx(Li[2], 2, 1);
            case 6:                                                                 
                return dNadx(Li[0], 0, 1)*Na(Li[2], 1) + Na(Li[0], 1)*dNadx(Li[2], 2, 1);
            case 7:                                                                 
                return dNadx(Li[0], 0, 1)*Na(Li[3], 1) + Na(Li[0], 1)*dNadx(Li[3], 3, 1);
            case 8:                                                                 
                return dNadx(Li[2], 2, 1)*Na(Li[3], 1) + Na(Li[2], 1)*dNadx(Li[3], 3, 1);
            case 9:                                                                 
                return dNadx(Li[1], 1, 1)*Na(Li[3], 1) + Na(Li[1], 1)*dNadx(Li[3], 3, 1);
        }
        return 0;
    }
    inline double dNdy(const std::array<double, 4>& Li, size_t i) const{
        switch(i){
            case 0:
                return dNady(Li[0], 0, 2);
            case 1:                  
                return dNady(Li[1], 1, 2);
            case 2:                  
                return dNady(Li[2], 2, 2);
            case 3:                  
                return dNady(Li[3], 3, 2);
            case 4:                  
                return dNady(Li[0], 0, 1)*Na(Li[1], 1) + Na(Li[0], 1)*dNady(Li[1], 1, 1);
            case 5:                                                                 
                return dNady(Li[1], 1, 1)*Na(Li[2], 1) + Na(Li[1], 1)*dNady(Li[2], 2, 1);
            case 6:                                                                 
                return dNady(Li[0], 0, 1)*Na(Li[2], 1) + Na(Li[0], 1)*dNady(Li[2], 2, 1);
            case 7:                                                                 
                return dNady(Li[0], 0, 1)*Na(Li[3], 1) + Na(Li[0], 1)*dNady(Li[3], 3, 1);
            case 8:                                                                 
                return dNady(Li[2], 2, 1)*Na(Li[3], 1) + Na(Li[2], 1)*dNady(Li[3], 3, 1);
            case 9:                                                                 
                return dNady(Li[1], 1, 1)*Na(Li[3], 1) + Na(Li[1], 1)*dNady(Li[3], 3, 1);
        }
        return 0;
    }
    inline double dNdz(const std::array<double, 4>& Li, size_t i) const{
        switch(i){
            case 0:
                return dNadz(Li[0], 0, 2);
            case 1:                  
                return dNadz(Li[1], 1, 2);
            case 2:                  
                return dNadz(Li[2], 2, 2);
            case 3:                  
                return dNadz(Li[3], 3, 2);
            case 4:                  
                return dNadz(Li[0], 0, 1)*Na(Li[1], 1) + Na(Li[0], 1)*dNadz(Li[1], 1, 1);
            case 5:                                                                 
                return dNadz(Li[1], 1, 1)*Na(Li[2], 1) + Na(Li[1], 1)*dNadz(Li[2], 2, 1);
            case 6:                                                                 
                return dNadz(Li[0], 0, 1)*Na(Li[2], 1) + Na(Li[0], 1)*dNadz(Li[2], 2, 1);
            case 7:                                                                 
                return dNadz(Li[0], 0, 1)*Na(Li[3], 1) + Na(Li[0], 1)*dNadz(Li[3], 3, 1);
            case 8:                                                                 
                return dNadz(Li[2], 2, 1)*Na(Li[3], 1) + Na(Li[2], 1)*dNadz(Li[3], 3, 1);
            case 9:                                                                 
                return dNadz(Li[1], 1, 1)*Na(Li[3], 1) + Na(Li[1], 1)*dNadz(Li[3], 3, 1);
        }
        return 0;
    }

    inline Eigen::Matrix<double, DIM, K_DIM> N_mat(const gp_Pnt& p) const{
        const auto Lv = La(p);

        Eigen::Matrix<double, DIM, K_DIM> NN;
        NN.fill(0);
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            const double Ni = N(Lv, i);
            NN(0,3*i+0) = Ni;
            NN(1,3*i+1) = Ni;
            NN(2,3*i+2) = Ni;
        }

        return NN;
    }
    inline Eigen::Matrix<double, S_SIZE, K_DIM> B_mat(const gp_Pnt& p) const{
        const auto Lv = La(p);

        Eigen::Matrix<double, S_SIZE, K_DIM> dNN;
        dNN.fill(0);
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            const double Nx = dNdx(Lv, i);
            const double Ny = dNdy(Lv, i);
            const double Nz = dNdz(Lv, i);
            dNN(0,3*i+0) = Nx;
            dNN(1,3*i+1) = Ny;
            dNN(2,3*i+2) = Nz;

            dNN(3,3*i+0) = Ny;
            dNN(3,3*i+1) = Nx;

            dNN(4,3*i+0) = Nz;
            dNN(4,3*i+2) = Nx;

            dNN(5,3*i+1) = Nz;
            dNN(5,3*i+2) = Ny;
        }

        return dNN;
    }

    inline Eigen::Vector<double, NODES_PER_ELEM> N_mat_1dof(const gp_Pnt& p) const{
        const auto Lv = La(p);
        Eigen::Vector<double, NODES_PER_ELEM> NN;
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            const double Ni = N(Lv, i);
            NN[i] = Ni;
        }

        return NN;
    }
    inline Eigen::Matrix<double, NODE_DOF, NODES_PER_ELEM> dN_mat_1dof(const gp_Pnt& p) const{
        const auto Lv = La(p);
        Eigen::Matrix<double, NODE_DOF, NODES_PER_ELEM> dNN;
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            dNN(0,i) = dNdx(Lv, i);
            dNN(1,i) = dNdy(Lv, i);
            dNN(2,i) = dNdz(Lv, i);
        }

        return dNN;
    }
    inline double dNdx_norm_surface(double x, double y, size_t i) const{
        (void)x;
        (void)y;
        // x, y, 1 - x - y
        switch(i){
            case 0:
                return  1;
            case 1:
                return  0;
            case 2:
                return -1;
        }
        return 0;
    }
    inline double dNdy_norm_surface(double x, double y, size_t i) const{
        (void)x;
        (void)y;
        // x, y, 1 - x - y
        switch(i){
            case 0:
                return 0;
            case 1:               
                return 1;
            case 2:               
                return -1;
        }
        return 0;
    }

    inline Eigen::Vector<double, NODE_DOF> surface_to_nat(double xi, double eta, const double A[3], const double B[3], const double C[3], const std::array<double, BOUNDARY_NODES_PER_ELEM>& x, const std::array<double, BOUNDARY_NODES_PER_ELEM>& y, const std::array<double, BOUNDARY_NODES_PER_ELEM>& z) const{
        double X = 0, Y = 0, Z = 0;
        for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
            const double Ni = A[i] + B[i]*xi + C[i]*eta;
            X += Ni*x[i];
            Y += Ni*y[i];
            Z += Ni*z[i];
        }

        return Eigen::Vector<double, NODE_DOF>{X, Y, Z};
    }
};

}

#endif
