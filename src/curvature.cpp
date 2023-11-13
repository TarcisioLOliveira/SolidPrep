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

#include "curvature.hpp"
#include <Eigen/Dense>
#include "logger.hpp"
#include "utils/gauss_legendre.hpp"

Curvature::Curvature(const Material* mat, gp_Dir u, gp_Dir v, gp_Dir w, Eigen::Matrix<double, 2, 2> rot2D, Eigen::Matrix<double, 3, 3> rot3D, const BoundaryMeshElementFactory* elem_info, double V_v, double V_w, double M_u, double M_v, double M_w):
    mat(mat), u(u), v(v), w(w),
    rot2D(rot2D), rot3D(rot3D), 
    elem_info(elem_info),
    V_v(V_v), V_w(V_w), 
    M_u(M_u), 
    M_v(M_v),
    M_w(M_w)
{

}

void Curvature::generate_curvature_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, size_t phi_size, size_t psi_size){
    std::array<gp_Pnt, 3> points;

    const auto fn_EA =
        [this](const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EA_base_3D(p, px);
        };
    const auto fn_EA_v =
        [this](const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EA_v_base_3D(p, px);
        };
    const auto fn_EA_w =
        [this](const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EA_w_base_3D(p, px);
        };
    const auto fn_EI_v =
        [this](const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EI_v_base_3D(p, px);
        };
    const auto fn_EI_w =
        [this](const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EI_w_base_3D(p, px);
        };
    const auto fn_EI_vw =
        [this](const gp_Pnt& p, const gp_Pnt& px)->double{
            return this->make_EI_vw_base_3D(p, px);
        };

    this->EA = this->integrate_surface_3D(boundary_mesh, fn_EA);
    const double EA_v = this->integrate_surface_3D(boundary_mesh, fn_EA_v);
    const double EA_w = this->integrate_surface_3D(boundary_mesh, fn_EA_w);
    this->c_v = EA_v/this->EA;
    this->c_w = EA_w/this->EA;
    this->EI_v = this->integrate_surface_3D(boundary_mesh, fn_EI_v);
    this->EI_w = this->integrate_surface_3D(boundary_mesh, fn_EI_w);
    this->EI_vw = this->integrate_surface_3D(boundary_mesh, fn_EI_vw);
    logger::quick_log("A: ", EA/180000);
    logger::quick_log("EA: ", EA);
    logger::quick_log("c_v: ", c_v);
    logger::quick_log("c_w: ", c_w);
    logger::quick_log("EI_v", EI_v/180000);
    logger::quick_log("EI_w", EI_w/180000);
    logger::quick_log("EI_vw", EI_vw/180000);

    const Eigen::Vector<double, 2> Mv{M_v, M_w};
    const Eigen::Vector<double, 2> Vv{V_v, V_w};
    const Eigen::Matrix<double, 2, 2> EI{{ EI_vw,  EI_w},
                                         {-EI_v, -EI_vw}};

    const Eigen::Vector<double, 2> cv = EI.fullPivLu().solve(Mv);
    this->curv_v = cv[0];
    this->curv_w = cv[1];
    const Eigen::Vector<double, 2> dcv = EI.fullPivLu().solve(Vv);
    this->dcurv_v = dcv[0];
    this->dcurv_w = dcv[1];

    this->phi_torsion.resize(phi_size,0);
    this->psi_shear.resize(psi_size,0);

    if(this->M_u != 0){
        this->calculate_torsion(boundary_mesh);
    }
    if(this->V_v != 0 || this->V_w != 0){
        this->calculate_shear_3D(boundary_mesh);
    }
}

void Curvature::get_shear_in_3D(const BoundaryMeshElement* e, double& t_uv, double& t_uw) const{
    auto c = e->get_centroid();

    Eigen::Vector<double, 2> grad = this->theta*e->grad_1dof(c, this->phi_torsion);
    t_uv = -grad[1];
    t_uw =  grad[0];
    grad = e->grad_1dof(c, this->psi_shear);
    const auto EG = this->mat->beam_EG_2D(c, this->u);
    t_uv += EG[1]*(grad[0])/2;
    t_uw += EG[2]*(grad[1])/2;
}

void Curvature::calculate_torsion(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh){
    const size_t num_nodes = this->elem_info->get_nodes_per_element();

    const size_t phi_size = this->phi_torsion.size();
    Eigen::SparseMatrix<double> M = Eigen::SparseMatrix<double>(phi_size, phi_size);
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::COLAMDOrdering<int>> solver;

    Eigen::VectorXd b;
    b.resize(phi_size);
    b.fill(0);

    Eigen::Matrix<double, 2, 2> G{{0,0},{0,0}};

    for(const auto& e:boundary_mesh){
        const auto N = 2*e->source_1dof();
        gp_Pnt c = utils::change_point(e->get_centroid(), this->rot3D);
        const auto EG = this->mat->beam_EG_3D(c, this->u);
        G(0,0) = 1.0/EG[2]; // G_uw
        G(1,1) = 1.0/EG[1]; // G_uv

        const auto M_e = e->diffusion_1dof(G);

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->u_pos[0];
            if(id1 < 0){
                continue;
            }
            for(size_t j = 0; j < num_nodes; ++j){
                const long id2 = e->nodes[j]->u_pos[0];
                if(id2 < 0){
                    continue;
                }
                if(std::abs(M_e(i,j)) > 1e-14){
                    M.coeffRef(id1, id2) += M_e(i, j);
                }
            }
        }

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->u_pos[0];
            if(id1 < 0){
                continue;
            }
            b[id1] += N[i];
        }
    }

    M.makeCompressed();
    solver.analyzePattern(M);

    solver.factorize(M);
    logger::log_assert(solver.info() == Eigen::Success, logger::ERROR, "matrix decomposition failed");
    Eigen::VectorXd phi_tmp = solver.solve(b);

    std::copy(phi_tmp.begin(), phi_tmp.end(), this->phi_torsion.begin());
    double phi_int = 0;
    for(const auto& e:boundary_mesh){
        const auto N = e->source_1dof();

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->u_pos[0];
            if(id1 < 0){
                continue;
            }
            phi_int += N[i]*this->phi_torsion[id1];
        }
    }
    this->theta = -this->M_u/(2*phi_int);
}

void Curvature::calculate_shear_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh){
    const size_t num_nodes = this->elem_info->get_nodes_per_element();

    const size_t psi_size = this->psi_shear.size();
    Eigen::SparseMatrix<double> M = Eigen::SparseMatrix<double>(psi_size, psi_size);
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::COLAMDOrdering<int>> solver;

    Eigen::VectorXd b;
    b.resize(psi_size);
    b.fill(0);

    Eigen::Matrix<double, 2, 2> G{{0,0},{0,0}};

    for(const auto& e:boundary_mesh){
        gp_Pnt c = utils::change_point(e->get_centroid(), this->rot3D);
        const auto EG = this->mat->beam_EG_3D(c, this->u);
        G(0,0) = EG[1]/2; // G_uv
        G(1,1) = EG[2]/2; // G_uw
        Eigen::Vector<double, 3> v1{0, dcurv_v, dcurv_w};
        const auto N1 = EG[0]*e->source_1dof(v1);
        Eigen::Vector<double, 2> v2{EG[1]*curv_v, EG[2]*curv_w};
        const auto N2 = e->source_grad_1dof(v2);
        const auto N = N1;// - N2;

        const auto M_e = e->diffusion_1dof(G);

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->id;
            if(id1 < 0){
                continue;
            }
            for(size_t j = 0; j < num_nodes; ++j){
                const long id2 = e->nodes[j]->id;
                if(id2 < 0){
                    continue;
                }
                if(std::abs(M_e(i,j)) > 1e-14){
                    M.coeffRef(id1, id2) += M_e(i, j);
                }
            }
        }

        for(size_t i = 0; i < num_nodes; ++i){
            const long id1 = e->nodes[i]->id;
            if(id1 < 0){
                continue;
            }
            b[id1] += N[i];
        }
    }

    M.makeCompressed();
    solver.analyzePattern(M);

    solver.factorize(M);
    logger::log_assert(solver.info() == Eigen::Success, logger::ERROR, "matrix decomposition failed");
    Eigen::VectorXd psi_tmp = solver.solve(b);

    std::copy(psi_tmp.begin(), psi_tmp.end(), this->psi_shear.begin());
}

double Curvature::integrate_surface_3D(const std::vector<std::unique_ptr<BoundaryMeshElement>>& boundary_mesh, const std::function<double(const gp_Pnt&, const gp_Pnt&)>& fn) const{
    double result = 0;
    if(this->elem_info->get_shape_type() == Element::Shape::TRI){
        #pragma omp parallel for reduction(+:result)
        for(size_t i = 0; i < boundary_mesh.size(); ++i){
            Eigen::Matrix<double, 3, 3> points{{0,0,0},{0,0,0},{0,0,0}};
            const auto& e = boundary_mesh[i];
            std::array<gp_Pnt, 3> rel_points;
            for(size_t x = 0; x < 3; ++x){
                rel_points[x] = e->nodes[x]->point;
                for(size_t y = 0; y < 3; ++y){
                    points(y, x) = rel_points[x].Coord(y+1);
                }
            }
            Eigen::Matrix<double, 3, 3> rotd_p = this->rot3D*points;
            std::array<gp_Pnt, 3> abs_points{
                gp_Pnt(rotd_p(0, 0), rotd_p(1, 0), rotd_p(2, 0)),
                gp_Pnt(rotd_p(0, 1), rotd_p(1, 1), rotd_p(2, 1)),
                gp_Pnt(rotd_p(0, 2), rotd_p(1, 2), rotd_p(2, 2))
            };
            result += this->GS_tri(abs_points, rel_points, fn);
        }
    } else if(this->elem_info->get_shape_type() == Element::Shape::QUAD){

    }
    return result;
}

double Curvature::GS_tri(const std::array<gp_Pnt, 3>& p, const std::array<gp_Pnt, 3>& px, const std::function<double(const gp_Pnt&, const gp_Pnt&)>& fn) const{
    double result = 0;
    gp_Vec v1(p[1], p[0]);
    gp_Vec v2(p[2], p[0]);
    const double drnorm = (v1.Crossed(v2)).Magnitude()/2;
    const auto& gsi = utils::GaussLegendreTri<6>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        gp_Pnt pi{
            it->a*p[0].X() + it->b*p[1].X() + it->c*p[2].X(),
            it->a*p[0].Y() + it->b*p[1].Y() + it->c*p[2].Y(),
            it->a*p[0].Z() + it->b*p[1].Z() + it->c*p[2].Z()
        };
        gp_Pnt pxi{
            it->a*px[0].X() + it->b*px[1].X() + it->c*px[2].X(),
            it->a*px[0].Y() + it->b*px[1].Y() + it->c*px[2].Y(),
            it->a*px[0].Z() + it->b*px[1].Z() + it->c*px[2].Z()
        };
        result += it->w*fn(pi, pxi);
    }
    result *= drnorm;

    return result;
}
