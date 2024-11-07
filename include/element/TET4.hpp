/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
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

#ifndef TET4_HPP
#define TET4_HPP

#include <Eigen/Core>
#include <vector>
#include <array>
#include "element.hpp"
#include "material.hpp"
#include "utils.hpp"
#include "element_factory.hpp"
#include "element_common.hpp"

namespace element{

class TET4 : public MeshElementCommon3DTet<TET4>{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 4;
    static const size_t NODE_DOF       = 3;
    static const size_t NODES_PER_ELEM = 4;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;

    static const size_t BOUNDARY_NODES_PER_ELEM = 3;
    static const size_t BOUNDARY_GMSH_TYPE = 2;
    static std::unique_ptr<BoundaryMeshElementFactory> get_boundary_element_info();
    static std::unique_ptr<ContactMeshElementFactory> get_contact_element_info();

    static const spview::defs::ElementType SPVIEW_CODE = spview::defs::TET4;

    TET4(ElementShape s);

    virtual std::vector<double> get_k(const std::vector<double>& D, const double t) const override;
    virtual std::vector<double> get_nodal_density_gradient(gp_Pnt p) const override;
    virtual std::vector<double> get_R(const std::vector<double>& K, const double t, const std::vector<gp_Pnt>& points) const override;
    virtual std::vector<double> get_Rf(const std::vector<double>& S, const std::vector<double>& F, const gp_Pnt& C, const double t, const std::vector<gp_Pnt>& points) const override;
    virtual std::vector<double> get_B(const gp_Pnt& point) const override;

    virtual std::vector<double> get_dk_sh(const std::vector<double>& D, const double t, const size_t n, const size_t dof) const override;
    virtual void calculate_coefficients() override;

    virtual Eigen::MatrixXd diffusion_1dof(const double t, const std::vector<double>& A) const override;
    virtual Eigen::MatrixXd advection_1dof(const double t, const std::vector<double>& v) const override;
    virtual Eigen::MatrixXd absorption_1dof(const double t) const override;
    virtual Eigen::VectorXd source_1dof(const double t) const override;
    virtual Eigen::VectorXd flow_1dof(const double t, const MeshNode** nodes) const override;

    virtual std::vector<double> get_Ni(const gp_Pnt& p) const override;

    virtual inline std::unique_ptr<MeshElementFactory> get_element_info() const override{
        return std::unique_ptr<MeshElementFactory>(new MeshElementFactoryImpl<TET4>());
    }

    private:
    virtual std::vector<double> get_DB(const std::vector<double>& D, const gp_Pnt& point) const override;
    virtual std::vector<double> get_Nf(const double t, const std::vector<gp_Pnt>& points) const override;

    typedef std::array<double, NODES_PER_ELEM*NODES_PER_ELEM> CoeffMat;

    CoeffMat C;
    double V;

    inline double N(double x, double y, double z, size_t i, const CoeffMat& M) const{
        const double* const a = M.data();
        const double* const b = a + NODES_PER_ELEM;
        const double* const c = b + NODES_PER_ELEM;
        const double* const d = c + NODES_PER_ELEM;
        return (a[i] + b[i]*x + c[i]*y + d[i]*z)/(6*V);
    }

    inline double N_norm(double x, double y, double z, size_t i) const{
        switch(i){
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
            case 3:
                return 1 - x - y - z;
        }
        return 0;
    }

    inline Eigen::Matrix<double, DIM, K_DIM> N_mat(double x, double y, double z, const CoeffMat& M) const{
        const double Ni[NODES_PER_ELEM] = {N(x, y, z, 0, M), N(x, y, z, 1, M), N(x, y, z, 2, M), N(x, y, z, 3, M)};

        return Eigen::Matrix<double, DIM, K_DIM>{{Ni[0], 0, 0, Ni[1], 0, 0, Ni[2], 0, 0, Ni[3], 0, 0},
                                                 {0, Ni[0], 0, 0, Ni[1], 0, 0, Ni[2], 0, 0, Ni[3], 0},
                                                 {0, 0, Ni[0], 0, 0, Ni[1], 0, 0, Ni[2], 0, 0, Ni[3]}};
    }

    inline Eigen::Vector<double, 4> N_mat_1dof(double x, double y, double z, const CoeffMat& M) const{
        return Eigen::Vector<double, 4>(N(x, y, z, 0, M), N(x, y, z, 1, M), N(x, y, z, 2, M), N(x, y, z, 3, M));
    }
    inline Eigen::Matrix<double, 3, 4> dN_mat_1dof(const CoeffMat& M) const{
        const double* const a = M.data();
        const double* const b = a + NODES_PER_ELEM;
        const double* const c = b + NODES_PER_ELEM;
        const double* const d = c + NODES_PER_ELEM;
        return Eigen::Matrix<double, 3, 4>{{b[0], b[1], b[2], b[3]},
                                           {c[0], c[1], c[2], c[3]},
                                           {d[0], d[1], d[2], d[3]}}/(6*V);
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


    Eigen::Matrix<double, S_SIZE, K_DIM> B(const CoeffMat& M) const{
        const double* const a = M.data();
        const double* const b = a + NODES_PER_ELEM;
        const double* const c = b + NODES_PER_ELEM;
        const double* const d = c + NODES_PER_ELEM;
        Eigen::Matrix<double, S_SIZE, K_DIM>  B{
            {b[0]/(6*V), 0, 0, b[1]/(6*V), 0, 0, b[2]/(6*V), 0, 0, b[3]/(6*V), 0, 0},
            {0, c[0]/(6*V), 0, 0, c[1]/(6*V), 0, 0, c[2]/(6*V), 0, 0, c[3]/(6*V), 0},
            {0, 0, d[0]/(6*V), 0, 0, d[1]/(6*V), 0, 0, d[2]/(6*V), 0, 0, d[3]/(6*V)},
            {c[0]/(6*V), b[0]/(6*V), 0, c[1]/(6*V), b[1]/(6*V), 0, c[2]/(6*V), b[2]/(6*V), 0, c[3]/(6*V), b[3]/(6*V), 0},
            {d[0]/(6*V), 0, b[0]/(6*V), d[1]/(6*V), 0, b[1]/(6*V), d[2]/(6*V), 0, b[2]/(6*V), d[3]/(6*V), 0, b[3]/(6*V)},
            {0, d[0]/(6*V), c[0]/(6*V), 0, d[1]/(6*V), c[1]/(6*V), 0, d[2]/(6*V), c[2]/(6*V), 0, d[3]/(6*V), c[3]/(6*V)}
        }; 
        return B;
    }

    CoeffMat get_C_derivative(const size_t n, const size_t dof) const;
};

}

#endif
