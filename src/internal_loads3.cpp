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

#include <algorithm>
#include <set>
#include "internal_loads3.hpp"
#include "math/matrix.hpp"
#include "math/slicer.hpp"
#include "meshing.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include "utils/basis_tensor.hpp"

inline math::Matrix Lek_rot3D(math::Matrix L, math::Matrix R){
    return R*L;//L*R;
}

InternalLoads3::InternalLoads3(CrossSection cross_section, double thickness, gp_Dir normal, gp_Dir v, gp_Dir w, utils::DelayedPointerView<Material> mat, std::array<double, 3> F, std::array<double, 3> M, utils::ProblemType type):
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
    rot2D_D(utils::basis_tensor_2D(rot2D.T())),
    rot2D_S(utils::basis_tensor_2D_inv_T(rot2D.T())),
    rot3D_D(utils::basis_tensor_3D(rot3D.T())),
    rot3D_S(utils::basis_tensor_3D_inv_T(rot3D.T())),
    F({F[0], F[1], F[2]}), M({M[0], M[1], M[2]}),
    mat(mat),
    normal(normal),
    v(v), w(w),
    type(type),
    mmat({mat}, type, false)
{

}


void InternalLoads3::initialize(const Meshing* const mesh){
    size_t bound_nodes_per_elem = mesh->elem_info->get_boundary_nodes_per_element();

    const auto& boundary_elements = mesh->boundary_elements;

    this->elem_info = mesh->elem_info;
    this->BOUNDARY_NODE_NUM = elem_info->get_boundary_nodes_per_element();
    const size_t node_num = elem_info->get_nodes_per_element();
    const size_t dof = elem_info->get_dof_per_node();
    const size_t s_size = elem_info->get_D_dimension();
    //const size_t dim = (this->type == utils::PROBLEM_TYPE_2D) ? 2 : 3;

    // Get associated boundary elements
    if(this->submesh.size() == 0){
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
                for(size_t n = 0; n < node_num; ++n){
                    this->F_id_map[boundary_elements[i].parent->nodes[n]->id] = -1;
                }
                for(size_t n = 0; n < BOUNDARY_NODE_NUM; ++n){
                    this->F_id_map[boundary_elements[i].nodes[n]->id] = 0;
                }

                ++cur_elem;
            }
        }

        size_t cur_id = 0;
        for(auto& n:F_id_map){
            if(n.second == 0){
                n.second = cur_id;
                ++cur_id;
            }
        }

        this->u_size = (cur_id*dof + s_size)*2;
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

    this->calculate_force_distribution();
}

void InternalLoads3::calculate_force_distribution(){
    const size_t K_DIM = this->elem_info->get_k_dimension();
    const size_t NODE_DOF = this->elem_info->get_dof_per_node();
    const size_t NODES_PER_ELEM = this->elem_info->get_nodes_per_element();
    const size_t S_SIZE = elem_info->get_D_dimension();
    const size_t DIM = (this->type == utils::PROBLEM_TYPE_2D) ? 2 : 3;

    if(this->u.size() != u_size){
        this->u.resize(u_size, 0); 
        this->solver.initialize_matrix(false, u_size);
    } else {
        std::fill(this->u.begin(), this->u.end(), 0.0);
        this->solver.make_zero();
    }

    const size_t V_offset = this->u.size() - 2*S_SIZE;
    const size_t f_offset = V_offset/2;

    std::vector<long> posi(K_DIM);
    std::vector<long> posf(K_DIM);
    std::vector<long> pos_M1(S_SIZE);
    std::vector<long> pos_M2(S_SIZE);

    std::vector<gp_Pnt> bounds(BOUNDARY_NODE_NUM);

    for(size_t i = 0; i < S_SIZE; ++i){
        pos_M1[i] = V_offset + i;
        pos_M2[i] = V_offset + S_SIZE + i;
    }

    math::Matrix D;
    math::Matrix Fe;
    math::Matrix M1;
    math::Matrix M2;

    long max_u = 0;
    long max_id = 0;

    const auto I = math::Matrix::identity(DIM);

    for(const auto& b:this->submesh){
        const gp_Pnt centroid = b->get_centroid(this->BOUNDARY_NODE_NUM);
        const gp_Pnt c = utils::change_point(centroid, this->rot3D);
        for(size_t i = 0; i < BOUNDARY_NODE_NUM; ++i){
            bounds[i] = b->nodes[i]->point;
        }

        if(this->type == utils::PROBLEM_TYPE_2D){
            D = this->mat->stiffness_2D(b->parent, c);
            //Fe = -b->parent->get_R(this->rot2D, this->thickness, bounds);
            //M1 = b->parent->constraint_matrix_alt(this->center, this->rot2D, bounds);
            Fe = b->parent->internal_loads_hess(this->center, this->rot2D, this->rot2D_D, D, this->normal, bounds);
            M1 = b->parent->constraint_matrix(this->center, this->rot2D, this->rot2D_D, this->mmat, bounds);
            M2 = b->parent->constraint_matrix(this->center, this->rot2D, this->rot2D_D, this->mmat, bounds);
        } else {
            D = this->mat->stiffness_3D(b->parent, c);
            //Fe = -b->parent->get_R(this->rot3D, this->thickness, bounds);
            //M1 = b->parent->constraint_matrix_alt(this->center, this->rot3D, bounds);
            Fe = b->parent->internal_loads_hess(this->center, this->rot3D, this->rot3D_D, D, this->normal, bounds);
            M1 = b->parent->constraint_matrix(this->center, this->rot3D, this->rot3D_D, this->mmat, bounds);
            M2 = b->parent->constraint_matrix(this->center, this->rot3D, this->rot3D_D, this->mmat, bounds);
        }

        //Fe = -b->parent->get_R(I, this->thickness, bounds);
        const auto Ke = b->parent->get_k_face(D, this->thickness, bounds) + b->parent->get_R(I, this->thickness, bounds);

        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            const auto f_id = this->F_id_map[b->parent->nodes[i]->id];
            if(f_id >= 0){
                for(size_t j = 0; j < NODE_DOF; ++j){
                    const long id = NODE_DOF*f_id + j;
                    posi[i*NODE_DOF + j] = id;
                    posf[i*NODE_DOF + j] = id + f_offset;
                }
            } else {
                for(size_t j = 0; j < NODE_DOF; ++j){
                    posi[i*NODE_DOF + j] = -1;
                    posf[i*NODE_DOF + j] = -1;
                }
            }
        }

        auto max_tmp = *std::max_element(posi.begin(), posi.end());
        max_u = std::max(max_u, max_tmp);
        max_tmp = *std::max_element(posf.begin(), posf.end());
        max_id = std::max(max_id, max_tmp);

        this->solver.add_element(Ke, posi);
        this->solver.add_element(Fe, posi, posf);
        this->solver.add_element(M1, pos_M1, posf);
        this->solver.add_element(M2, pos_M2, posi);
    }

    this->solver.compute();

    const auto V = this->get_loads();

    for(size_t i = 0; i < S_SIZE; ++i){
        this->u[V_offset + i] = V[i];
        this->u[V_offset + S_SIZE + i] = V[i];
    }

    this->solver.solve(this->u);
}

void InternalLoads3::apply_load(const std::vector<long>& node_position, std::vector<double>& load) const{
    const size_t K_DIM = this->elem_info->get_k_dimension();
    const size_t NODE_DOF = this->elem_info->get_dof_per_node();
    const size_t NODES_PER_ELEM = this->elem_info->get_nodes_per_element();
    //const size_t S_SIZE = elem_info->get_D_dimension();
    //const size_t DIM = (this->type == utils::PROBLEM_TYPE_2D) ? 2 : 3;

    math::Vector ue(K_DIM);
    std::vector<gp_Pnt> bounds(BOUNDARY_NODE_NUM);
    math::Matrix D;
    math::Matrix Fe;
    math::Matrix M1;
    math::Matrix M2;
    if(this->type == utils::PROBLEM_TYPE_2D){

    } else {
        for(const auto& b:this->submesh){
            const gp_Pnt centroid = b->get_centroid(this->BOUNDARY_NODE_NUM);
            const gp_Pnt c = utils::change_point(centroid, this->rot3D);
            for(size_t i = 0; i < BOUNDARY_NODE_NUM; ++i){
                bounds[i] = b->nodes[i]->point;
            }

            if(this->type == utils::PROBLEM_TYPE_2D){
                D = this->mat->stiffness_2D(b->parent, c);
            } else {
                D = this->mat->stiffness_3D(b->parent, c);
            }

            for(size_t i = 0; i < NODES_PER_ELEM; ++i){
                const auto f_id = this->F_id_map.at(b->parent->nodes[i]->id);
                if(f_id >= 0){
                    for(size_t j = 0; j < NODE_DOF; ++j){
                        const long id = NODE_DOF*f_id + j;
                        ue[i*NODE_DOF + j] = this->u[id];
                    }
                } else {
                    for(size_t j = 0; j < NODE_DOF; ++j){
                        ue[i*NODE_DOF + j] = 0;
                    }
                }
            }

            const auto fe = b->parent->stress_to_force(D, this->normal, ue, this->thickness, bounds);

            for(size_t i = 0; i < NODES_PER_ELEM; ++i){
                for(size_t j = 0; j < NODE_DOF; ++j){
                    const auto u_id = b->parent->nodes[i]->u_pos[j];
                    const auto u_pos = node_position[u_id];
                    if(u_pos >= 0){
                        load[u_pos] += fe[i*NODE_DOF + j];
                    }
                }
            }

        }
    }
}

math::Vector InternalLoads3::calculate_internal_loads(const std::vector<double>& u_ext) const{
    math::Vector result(6);
    const size_t bnum = this->elem_info->get_boundary_nodes_per_element();
    std::vector<gp_Pnt> bounds(bnum);

    if(this->type == utils::PROBLEM_TYPE_2D){
        for(const auto& e:this->submesh){
            for(size_t i = 0; i < bnum; ++i){
                bounds[i] = e->nodes[i]->point;
            }
            const auto D = this->mat->stiffness_2D(e->parent, e->parent->get_centroid());
            result += e->parent->calculate_internal_loads(this->center, this->rot2D, this->rot2D_D, this->mmat, this->normal, u_ext, bounds);
        }
    } else {
        for(const auto& e:this->submesh){
            for(size_t i = 0; i < bnum; ++i){
                bounds[i] = e->nodes[i]->point;
            }
            const auto D = this->mat->stiffness_3D(e->parent, e->parent->get_centroid());
            result += e->parent->calculate_internal_loads(this->center, this->rot3D, this->rot3D_D, this->mmat, this->normal, u_ext, bounds);
        }
    }

    return result;
}

math::Vector InternalLoads3::get_loads() const{

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

math::Vector InternalLoads3::integrate_surface_3D(const Meshing* const mesh, const std::vector<std::function<double(const math::Matrix& S, const gp_Pnt& px)>>& fn) const{
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

void InternalLoads3::GS_tri(const MeshElement* const e, const std::array<gp_Pnt, 3>& p, const std::vector<std::function<double(const math::Matrix& S, const gp_Pnt& px)>>& fn, math::Vector& result) const{
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

void InternalLoads3::GS_quad(const MeshElement* const e, const std::array<gp_Pnt, 4>& p, const std::vector<std::function<double(const math::Matrix& S, const gp_Pnt& px)>>& fn, math::Vector& result) const{
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
