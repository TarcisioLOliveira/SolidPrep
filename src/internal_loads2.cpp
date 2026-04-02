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

#include "internal_loads2.hpp"
#include "math/matrix.hpp"
#include "meshing.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include "utils/basis_tensor.hpp"

inline math::Matrix Lek_rot3D(math::Matrix L, math::Matrix R){
    return R*L;//L*R;
}

InternalLoads2::InternalLoads2(CrossSection cross_section, double thickness, gp_Dir normal, gp_Dir v, gp_Dir w, utils::DelayedPointerView<Material> mat, std::array<double, 3> F, std::array<double, 3> M, utils::ProblemType type):
    S(std::move(cross_section)), 
    A(S.get_area()),
    thickness(thickness),
    rot2D({normal.X(), v.X(),
           normal.Y(), v.Y()}, 2, 2),
    Lek_basis(
        {0, 0, 1,
         0, 1, 0,
         -1, 0, 0}, 3, 3),
    rot3D(Lek_rot3D(Lek_basis,
        math::Matrix(
        {normal.X(), v.X(), w.X(),
         normal.Y(), v.Y(), w.Y(),
         normal.Z(), v.Z(), w.Z()}, 3, 3))),
    rot2D_D(utils::basis_tensor_2D(rot2D)),
    rot2D_S(utils::basis_tensor_2D_inv_T(rot2D)),
    rot3D_D(utils::basis_tensor_3D(rot3D)),
    rot3D_S(utils::basis_tensor_3D_inv_T(rot3D)),
    F({F[0], F[1], F[2]}), M({M[0], M[1], M[2]}),
    mat(mat),
    normal(normal),
    MU(1),
    v(v), w(w),
    type(type){

}


void InternalLoads2::initialize(const Meshing* const mesh, size_t il_id, std::vector<double>& load){
    size_t bound_nodes_per_elem = mesh->elem_info->get_boundary_nodes_per_element();

    const auto& boundary_elements = mesh->boundary_elements;

    this->elem_info = mesh->elem_info;
    this->BOUNDARY_NODE_NUM = elem_info->get_boundary_nodes_per_element();

    // Get associated boundary elements

    std::vector<bool> apply_spring(boundary_elements.size());
    size_t num_elems = 0;
    // Fill apply_spring concurrently
    // (CrossSection::is_inside() is pretty slow, unfortunately
    #pragma omp parallel for reduction(+:num_elems)
    for(size_t i = 0; i < boundary_elements.size(); ++i){
        const auto& e = boundary_elements[i];
        apply_spring[i] = this->S.is_inside(e.get_centroid(bound_nodes_per_elem));
        if(apply_spring[i]){
            ++num_elems;
        }
    }

    this->submesh.resize(num_elems);
    size_t cur_elem = 0;
    ElementShape sh;
    
    sh.nodes.resize(bound_nodes_per_elem);
    auto boundary_elem_info = mesh->elem_info->get_boundary_element_info();
    for(size_t i = 0; i < boundary_elements.size(); ++i){
        if(apply_spring[i]){
            this->submesh[cur_elem] = &boundary_elements[i];

            ++cur_elem;
        }
    }

    const auto fn_A =
        [this](const math::Matrix& S, const gp_Pnt& px)->double{
            return this->make_A_base_3D(S, px);
        };
    const auto fn_EA =
        [this](const math::Matrix& S, const gp_Pnt& px)->double{
            return this->make_EA_base_3D(S, px);
        };
    const auto fn_EA_u =
        [this](const math::Matrix& S, const gp_Pnt& px)->double{
            return this->make_EA_u_base_3D(S, px);
        };
    const auto fn_EA_v =
        [this](const math::Matrix& S, const gp_Pnt& px)->double{
            return this->make_EA_v_base_3D(S, px);
        };
    const auto fn_EA_w =
        [this](const math::Matrix& S, const gp_Pnt& px)->double{
            return this->make_EA_w_base_3D(S, px);
        };

    // Calculate neutral axis position
    const auto result1 = this->integrate_surface_3D(mesh, {fn_A, fn_EA, fn_EA_u, fn_EA_v, fn_EA_w});
    double A    = result1[0];
    double EA   = result1[1];
    double EA_u = result1[2];
    double EA_v = result1[3];
    double EA_w = result1[4];
    this->center[0] = EA_u/EA;
    this->center[1] = EA_v/EA;
    this->center[2] = EA_w/EA;

    this->center = this->rot3D.T()*this->center;

    logger::quick_log("A: ", A);
    logger::quick_log("EA: ", EA);
    logger::quick_log("EA_u", EA_u);
    logger::quick_log("EA_v", EA_v);
    logger::quick_log("EA_w", EA_w);
    logger::quick_log("c_u: ", center[0]);
    logger::quick_log("c_v: ", center[1]);


    if(this->type == utils::PROBLEM_TYPE_2D){

    } else {
        //const size_t bnum = this->elem_info->get_boundary_nodes_per_element();
        //const size_t node_num = this->elem_info->get_nodes_per_element();
        //const size_t dof = this->elem_info->get_dof_per_node();
        //std::vector<gp_Pnt> bounds(bnum);

        //math::Matrix M(6,6);
        //for(const auto& b:this->submesh){
        //    for(size_t i = 0; i < bnum; ++i){
        //        bounds[i] = b->nodes[i]->point;
        //    }
        //    const auto D = this->mat->stiffness_3D(b->parent, b->parent->get_centroid());
        //    M += b->parent->il2_strain_matrix(center, this->rot3D, rot3D_D, D, bounds);
        //}

        //math::LU lu(M);
        math::Vector V = this->get_loads();
        //const size_t il_size = mesh->internal_loads2->size();
        for(size_t i = 0; i < 6; ++i){
            //load[mesh->il2_offset + 6*(il_size + il_id) + i] = V[i];
            load[mesh->il2_offset + 6*il_id + i] = V[i];
        }

        //lu.solve(V);
        //this->base_strain = V;
        //logger::quick_log("Base strain:");
        //logger::quick_log(this->base_strain);

        //for(const auto& b:this->submesh){
        //    for(size_t i = 0; i < bnum; ++i){
        //        bounds[i] = b->nodes[i]->point;
        //    }
        //    const auto D = this->mat->stiffness_3D(b->parent, b->parent->get_centroid());
        //    const auto fe = b->parent->il2_strain_load(center, rot3D, rot3D_D, D, base_strain, bounds);
        //    for(size_t i = 0; i < node_num; ++i){
        //        const size_t id = b->parent->nodes[i]->id;
        //        for(size_t j = 0; j < dof; ++j){
        //            const long p = mesh->node_positions[0][dof*id + j];
        //            if(p >= 0){
        //                load[p] -= fe[dof*i + j];
        //            }
        //        }
        //    }
        //}
    }
}

math::Vector InternalLoads2::calculate_internal_loads(const std::vector<double>& u_ext) const{
    math::Vector result(6);
    const size_t bnum = this->elem_info->get_boundary_nodes_per_element();
    std::vector<gp_Pnt> bounds(bnum);

    if(this->type == utils::PROBLEM_TYPE_2D){
        for(const auto& e:this->submesh){
            for(size_t i = 0; i < bnum; ++i){
                bounds[i] = e->nodes[i]->point;
            }
            const auto D = this->mat->stiffness_2D(e->parent, e->parent->get_centroid());
            result += e->parent->calculate_internal_loads(this->center, this->rot2D, this->rot2D_D, D, this->normal, u_ext, params, bounds);
        }
    } else {
        for(const auto& e:this->submesh){
            for(size_t i = 0; i < bnum; ++i){
                bounds[i] = e->nodes[i]->point;
            }
            const auto D = this->mat->stiffness_3D(e->parent, e->parent->get_centroid());
            result += e->parent->calculate_internal_loads(this->center, this->rot3D, this->rot3D_D, D, this->normal, u_ext, params, bounds);
        }
    }

    return result;
}

math::Vector InternalLoads2::calculate_constraints(const std::vector<double>& u_ext, bool print) const{
    auto result = this->calculate_internal_loads(u_ext);

    if(print){
        logger::quick_log("Internal loads:");
        logger::quick_log("My Mx Fz Mz Fx Fy");
        logger::quick_log(result);
    }

    result -= this->get_loads();

    if(print){
        logger::quick_log("IL constraints:");
        logger::quick_log(result);
    }

    return result;
}

math::Vector InternalLoads2::calculate_multipliers(const std::vector<double>& u_ext, bool print) const{
    auto c = this->calculate_constraints(u_ext, print);
    const math::Vector result = this->lambda + this->MU*c;
    if(print){
        logger::quick_log("IL multipliers:");
        logger::quick_log("My Mx Fz Mz Fx Fy");
        logger::quick_log(result);
    }

    return result;
}
math::Vector InternalLoads2::get_loads() const{

    const auto f = this->Lek_basis.T()*this->F;
    const auto m = this->Lek_basis.T()*this->M;

    return {
        m[1],
        m[0],
        f[2],
        m[2],
        f[0],
        f[1]
    };
}
void InternalLoads2::update_multipliers(const math::Vector& params_diff, double alpha, bool print){
    this->params += alpha*params_diff;
    if(print){
        logger::quick_log("IL params:");
        logger::quick_log("A B C theta A' B':");
        //logger::quick_log("C Kxz Kyz A B theta");
        logger::quick_log(this->params);
    }
}

math::Matrix InternalLoads2::get_constr(const BoundaryElement* b) const{
    std::vector<gp_Pnt> bounds(BOUNDARY_NODE_NUM);
    for(size_t i = 0; i < BOUNDARY_NODE_NUM; ++i){
        bounds[i] = b->nodes[i]->point;
    }
    if(this->type == utils::PROBLEM_TYPE_2D){
        const auto D = this->mat->stiffness_2D(b->parent, b->get_centroid(BOUNDARY_NODE_NUM));
        return b->parent->constraint_matrix(this->center, this->rot2D, this->rot2D_D, D, bounds);
    } else {
        const auto D = this->mat->stiffness_3D(b->parent, b->get_centroid(BOUNDARY_NODE_NUM));
        return b->parent->constraint_matrix(this->center, this->rot3D, this->rot3D_D, D, bounds);
    }
}
math::Vector InternalLoads2::get_grad(const BoundaryElement* b, const math::Vector& mult) const{
    std::vector<gp_Pnt> bounds(BOUNDARY_NODE_NUM);
    for(size_t i = 0; i < BOUNDARY_NODE_NUM; ++i){
        bounds[i] = b->nodes[i]->point;
    }
    if(this->type == utils::PROBLEM_TYPE_2D){
        const auto D = this->mat->stiffness_2D(b->parent, b->get_centroid(BOUNDARY_NODE_NUM));
        return b->parent->internal_loads_grad(mult, this->center, this->rot2D, this->rot2D_D, D, bounds);
    } else {
        const auto D = this->mat->stiffness_3D(b->parent, b->get_centroid(BOUNDARY_NODE_NUM));
        return b->parent->internal_loads_grad(mult, this->center, this->rot3D, this->rot3D_D, D, bounds);
    }
}
math::Matrix InternalLoads2::get_hess(const BoundaryElement* b) const{
    std::vector<gp_Pnt> bounds(BOUNDARY_NODE_NUM);
    for(size_t i = 0; i < BOUNDARY_NODE_NUM; ++i){
        bounds[i] = b->nodes[i]->point;
    }
    if(this->type == utils::PROBLEM_TYPE_2D){
        const auto D = this->mat->stiffness_2D(b->parent, b->get_centroid(BOUNDARY_NODE_NUM));
        auto hess =  b->parent->internal_loads_hess(center, rot2D, rot2D_D, D, this->normal, bounds);
        //const auto load = b->parent->il2_param_load(center, rot2D, rot2D_D, D, normal, bounds);
        //const size_t K_DIM = load.get_H();
        //const size_t S_SIZE = load.get_W();
        //for(size_t i = 0; i < K_DIM; ++i){
        //    for(size_t j = 0; j < S_SIZE; ++j){
        //       hess(i, j + K_DIM) = load(i, j);
        //       hess(j + K_DIM, i) = load(i, j);
        //    }
        //}

        return hess; 
    } else {
        const auto D = this->mat->stiffness_3D(b->parent, b->get_centroid(BOUNDARY_NODE_NUM));
        auto hess = b->parent->internal_loads_hess(center, rot3D, rot3D_D, D, this->normal, bounds);
        //const auto load = b->parent->il2_param_load(center, rot3D, rot3D_D, D, normal, bounds);
        //const size_t K_DIM = load.get_H();
        //const size_t S_SIZE = load.get_W();
        //for(size_t i = 0; i < K_DIM; ++i){
        //    for(size_t j = 0; j < S_SIZE; ++j){
        //       hess(i, j + K_DIM) = load(i, j);
        //       hess(j + K_DIM, i) = load(i, j);
        //    }
        //} 

        return hess; 
    }
}

class VectorWrapper{
    public:
    math::Vector vec;
    VectorWrapper() = default;
    VectorWrapper(size_t s):vec(s){}

    VectorWrapper& operator+=(const VectorWrapper& v){
        if(this->vec.get_N() == 0){
            this->vec = math::Vector(v.vec.get_N());
        }
        this->vec += v.vec;

        return *this;
    }
};

math::Vector InternalLoads2::integrate_surface_3D(const Meshing* const mesh, const std::vector<std::function<double(const math::Matrix& S, const gp_Pnt& px)>>& fn) const{
    VectorWrapper result(fn.size());
    #pragma omp declare reduction(vecsum : VectorWrapper : omp_out += omp_in) 
    if(mesh->elem_info->get_shape_type() == Element::Shape::TRI){
        constexpr size_t N = 3;
        #pragma omp parallel
        {
            VectorWrapper result_tmp(fn.size());
            #pragma omp for reduction(vecsum:result)
            for(size_t i = 0; i < submesh.size(); ++i){
                const auto& e = submesh[i];
                std::array<gp_Pnt, N> abs_points;
                for(size_t x = 0; x < N; ++x){
                    abs_points[x] = e->nodes[x]->point;
                }
                this->GS_tri(e->parent, abs_points, fn, result_tmp.vec);
                result += result_tmp;
            }
        }
    } else if(mesh->elem_info->get_shape_type() == Element::Shape::QUAD){
        constexpr size_t N = 4;
        #pragma omp parallel
        {
            VectorWrapper result_tmp(fn.size());
            #pragma omp for reduction(vecsum:result)
            for(size_t i = 0; i < submesh.size(); ++i){
                const auto& e = submesh[i];
                std::array<gp_Pnt, N> abs_points;
                for(size_t x = 0; x < N; ++x){
                    abs_points[x] = e->nodes[x]->point;
                }
                this->GS_quad(e->parent, abs_points, fn, result_tmp.vec);
                result += result_tmp;
            }
        }
    }
    return result.vec;
}

void InternalLoads2::GS_tri(const MeshElement* const e, const std::array<gp_Pnt, 3>& p, const std::vector<std::function<double(const math::Matrix& S, const gp_Pnt& px)>>& fn, math::Vector& result) const{
    gp_Vec v1(p[1], p[0]);
    gp_Vec v2(p[2], p[0]);
    result.fill(0);
    const gp_Pnt centroid(
            (p[0].X() + p[1].X() + p[2].X())/3.0,
            (p[0].Y() + p[1].Y() + p[2].Y())/3.0,
            (p[0].Z() + p[1].Z() + p[2].Z())/3.0);
    math::Matrix S;
    if(this->type == utils::PROBLEM_TYPE_2D){
        S = this->mat->stiffness_inverse_2D(e, centroid);
        this->mat->rotate_S(S, this->rot2D_S);
    } else {
        S = this->mat->stiffness_inverse_3D(e, centroid);
        this->mat->rotate_S(S, this->rot3D_S);
    }
    const double drnorm = (v1.Crossed(v2)).Magnitude()/2;
    const auto& gsi = utils::GaussLegendreTri<2>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        for(size_t i = 0; i < fn.size(); ++i){
            gp_Pnt pxi(this->GLT_point(*it, p));
            result[i] += it->w*fn[i](S, pxi);
        }
    }
    result *= drnorm;
}

void InternalLoads2::GS_quad(const MeshElement* const e, const std::array<gp_Pnt, 4>& p, const std::vector<std::function<double(const math::Matrix& S, const gp_Pnt& px)>>& fn, math::Vector& result) const{
    result.fill(0);
    const gp_Pnt centroid(
            (p[0].X() + p[1].X() + p[2].X() + p[3].X())/4.0,
            (p[0].Y() + p[1].Y() + p[2].Y() + p[3].Y())/4.0,
            (p[0].Z() + p[1].Z() + p[2].Z() + p[3].Z())/4.0);
    math::Matrix S;
    if(this->type == utils::PROBLEM_TYPE_2D){
        S = this->mat->stiffness_inverse_2D(e, centroid);
        this->mat->rotate_S(S, this->rot2D_S);
    } else {
        S = this->mat->stiffness_inverse_3D(e, centroid);
        this->mat->rotate_S(S, this->rot3D_S);
    }
    const auto& gsi = utils::GaussLegendre<2>::get();
    for(auto xi = gsi.begin(); xi < gsi.end(); ++xi){
        for(auto eta = gsi.begin(); eta < gsi.end(); ++eta){
            for(size_t i = 0; i < fn.size(); ++i){
                gp_Pnt pxi(this->GL_point(*xi, *eta, p));
                const double detJ = std::abs(this->J(xi->x, eta->x, p).determinant());
                result[i] += xi->w*eta->w*detJ*fn[i](S, pxi);
            }
        }
    }
}
