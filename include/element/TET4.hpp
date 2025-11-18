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

#include "element.hpp"
#include "element_factory.hpp"
#include "element_common.hpp"
#include "math/matrix.hpp"

namespace element{

class TET4 : public MeshElementCommon3DTet<TET4>{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 4;
    static const size_t NODE_DOF       = 3;
    static const size_t NODES_PER_ELEM = 4;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;
    static const size_t INTEG_ORDER    = 1;

    static const size_t BOUNDARY_NODES_PER_ELEM = 3;
    static const size_t BOUNDARY_GMSH_TYPE = 2;
    static std::unique_ptr<BoundaryMeshElementFactory> get_boundary_element_info();
    static std::unique_ptr<ContactMeshElementFactory> get_contact_element_info();

    static const spview::defs::ElementType SPVIEW_CODE = spview::defs::TET4;

    TET4(ElementShape s);

    virtual math::Matrix get_k(const math::Matrix& D, const double t) const override;
    virtual math::Matrix get_nodal_density_gradient(gp_Pnt p) const override;
    virtual math::Matrix get_R(const math::Matrix& K, const double t, const std::vector<gp_Pnt>& points) const override;
    virtual math::Matrix get_B(const gp_Pnt& point) const override;

    virtual math::Matrix get_dk_sh(const math::Matrix& D, const double t, const size_t n, const size_t dof) const override;
    virtual math::Matrix get_dB_sh(const gp_Pnt& p, const size_t n, const size_t dof) const override;
    virtual math::Matrix get_dN_sh(const gp_Pnt& p, const size_t n, const size_t dof) const override;
    virtual void calculate_coefficients() override;

    virtual math::Matrix diffusion_1dof(const double t, const math::Matrix& A) const override;
    virtual math::Matrix advection_1dof(const double t, const math::Vector& v) const override;
    virtual math::Matrix absorption_1dof(const double t) const override;
    virtual math::Matrix robin_1dof(const double t, const std::vector<gp_Pnt>& points) const override;
    virtual math::Vector source_1dof(const double t) const override;
    virtual math::Vector flow_1dof(const double t, const MeshNode** nodes) const override;

    virtual math::Matrix get_Ni(const gp_Pnt& p) const override;
    virtual math::Vector get_Ni_1dof(const gp_Pnt& p) const override;
    virtual math::Matrix get_dNi(const gp_Pnt& p) const override{
        (void) p;
        return this->dN_mat_dim_dof(this->C);
    }

    virtual inline std::unique_ptr<MeshElementFactory> get_element_info() const override{
        return std::unique_ptr<MeshElementFactory>(new MeshElementFactoryImpl<TET4>());
    }

    private:
    static const bool reg;
    virtual math::Matrix get_Nf(const double t, const std::vector<gp_Pnt>& points) const override;

    math::Matrix C;
    double V;

    inline gp_Pnt GL_point_tri(double c1, double c2, double c3, const std::vector<gp_Pnt>& p) const{
        return gp_Pnt(
            c1*p[0].X() + c2*p[1].X() + c3*p[2].X(),
            c1*p[0].Y() + c2*p[1].Y() + c3*p[2].Y(),
            c1*p[0].Z() + c2*p[1].Z() + c3*p[2].Z()
        );
    }

    inline gp_Pnt GL_point(double c1, double c2, double c3, double c4) const{
        return gp_Pnt(
            c1*this->nodes[0]->point.X() + c2*this->nodes[1]->point.X() + c3*this->nodes[2]->point.X() + c4*this->nodes[3]->point.X(),
            c1*this->nodes[0]->point.Y() + c2*this->nodes[1]->point.Y() + c3*this->nodes[2]->point.Y() + c4*this->nodes[3]->point.Y(),
            c1*this->nodes[0]->point.Z() + c2*this->nodes[1]->point.Z() + c3*this->nodes[2]->point.Z() + c4*this->nodes[3]->point.Z()
        );
    }

    inline double N(double x, double y, double z, size_t i, const math::Matrix& M) const{
        const double* const a = M.data();
        const double* const b = a + NODES_PER_ELEM;
        const double* const c = b + NODES_PER_ELEM;
        const double* const d = c + NODES_PER_ELEM;
        return a[i] + b[i]*x + c[i]*y + d[i]*z;
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

    inline math::Matrix N_mat(double x, double y, double z, const math::Matrix& M) const{
        const double Ni[NODES_PER_ELEM] = {N(x, y, z, 0, M), N(x, y, z, 1, M), N(x, y, z, 2, M), N(x, y, z, 3, M)};

        return math::Matrix({Ni[0], 0, 0, Ni[1], 0, 0, Ni[2], 0, 0, Ni[3], 0, 0,
                             0, Ni[0], 0, 0, Ni[1], 0, 0, Ni[2], 0, 0, Ni[3], 0,
                             0, 0, Ni[0], 0, 0, Ni[1], 0, 0, Ni[2], 0, 0, Ni[3]}, DIM, K_DIM);
    }

    inline math::Vector N_mat_1dof(double x, double y, double z, const math::Matrix& M) const{
        return math::Vector({N(x, y, z, 0, M), N(x, y, z, 1, M), N(x, y, z, 2, M), N(x, y, z, 3, M)});
    }
    inline math::Matrix dN_mat_1dof(const math::Matrix& M) const{
        const double* const a = M.data();
        const double* const b = a + NODES_PER_ELEM;
        const double* const c = b + NODES_PER_ELEM;
        const double* const d = c + NODES_PER_ELEM;
        return math::Matrix({b[0], b[1], b[2], b[3],
                             c[0], c[1], c[2], c[3],
                             d[0], d[1], d[2], d[3]}, 3, 4);
    }
    inline math::Matrix dN_mat_dim_dof(const math::Matrix& M) const{
        const double* const a = M.data();
        const double* const b = a + NODES_PER_ELEM;
        const double* const c = b + NODES_PER_ELEM;
        const double* const d = c + NODES_PER_ELEM;
        math::Matrix dN(DIM*DIM, DIM*NODES_PER_ELEM, 0);
        for(size_t i = 0; i < DIM; ++i){
            for(size_t j = 0; j < NODES_PER_ELEM; ++j){
                dN(DIM*i + 0, j*DIM + i) = b[j];
                dN(DIM*i + 1, j*DIM + i) = c[j];
                dN(DIM*i + 2, j*DIM + i) = d[j];
            }
        }
        return dN;
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

    inline math::Vector surface_to_nat(double xi, double eta, const double A[3], const double B[3], const double C[3], const std::array<double, BOUNDARY_NODES_PER_ELEM>& x, const std::array<double, BOUNDARY_NODES_PER_ELEM>& y, const std::array<double, BOUNDARY_NODES_PER_ELEM>& z) const{
        double X = 0, Y = 0, Z = 0;
        for(size_t i = 0; i < BOUNDARY_NODES_PER_ELEM; ++i){
            const double Ni = A[i] + B[i]*xi + C[i]*eta;
            X += Ni*x[i];
            Y += Ni*y[i];
            Z += Ni*z[i];
        }

        return math::Vector{X, Y, Z};
    }


    inline math::Matrix B(const math::Matrix& M) const{
        const double* const a = M.data();
        const double* const b = a + NODES_PER_ELEM;
        const double* const c = b + NODES_PER_ELEM;
        const double* const d = c + NODES_PER_ELEM;
        math::Matrix B({
            b[0], 0, 0, b[1], 0, 0, b[2], 0, 0, b[3], 0, 0,
            0, c[0], 0, 0, c[1], 0, 0, c[2], 0, 0, c[3], 0,
            0, 0, d[0], 0, 0, d[1], 0, 0, d[2], 0, 0, d[3],
            c[0], b[0], 0, c[1], b[1], 0, c[2], b[2], 0, c[3], b[3], 0,
            d[0], 0, b[0], d[1], 0, b[1], d[2], 0, b[2], d[3], 0, b[3],
            0, d[0], c[0], 0, d[1], c[1], 0, d[2], c[2], 0, d[3], c[3]
        }, S_SIZE, K_DIM); 
        return B;
    }

    math::Matrix get_C_derivative(const size_t n, const size_t dof) const;
};

}

#endif
