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

#include <map>
#include <set>
#include "internal_loads2.hpp"
#include "visualization.hpp"
#include "math/matrix.hpp"
#include "meshing.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include "utils/basis_tensor.hpp"
#include "element_common.hpp" // Get MAX_GN

inline math::Matrix Lek_rot3D(math::Matrix L, math::Matrix R){
    return R*L;//L*R;
    //return L*R;
}

const static size_t MAX_GN = (element::MAX_GN+1)/2;

InternalLoads2::InternalLoads2(CrossSection cross_section, double thickness, gp_Dir normal, gp_Dir v, gp_Dir w, utils::DelayedPointerView<Material> mat, std::array<double, 3> F, std::array<double, 3> M, bool use_M, std::array<double, 3> P_last, bool inverted, bool calculate, utils::ProblemType type):
    S(std::move(cross_section)), 
    A(S.get_area()),
    thickness(thickness),
    R0({normal.X(), v.X(), w.X(),
         normal.Y(), v.Y(), w.Y(),
         normal.Z(), v.Z(), w.Z()}, 3, 3),
    rot2D({normal.X(), v.X(),
           normal.Y(), v.Y()}, 2, 2),
    Lek_basis(
        { 0, 0, 1,
          0, 1, 0,
         -1, 0, 0}, 3, 3),
    rot3D(Lek_rot3D(Lek_basis, R0)),
    rot2D_D(utils::basis_tensor_2D(rot2D)),
    rot2D_S(utils::basis_tensor_2D_inv_T(rot2D)),
    rot3D_D(utils::basis_tensor_3D(rot3D)),
    rot3D_S(utils::basis_tensor_3D_inv_T(rot3D)),
    F(3), M(3),
    use_M(use_M),
    mat(mat),
    normal(normal),
    v(v), w(w),
    type(type),
    F_orig(R0*math::Vector({F[0], F[1], F[2]})),
    M_orig(R0*math::Vector({M[0], M[1], M[2]})),
    P_prev({P_last[0], P_last[1], P_last[2]}),
    inverted_calc(inverted),
    calculate_from_last_point(calculate),
    mmat({mat}, type, false)
{
    this->F = rot3D.T()*this->F_orig;
}

class Edge{
    public:
    Edge(const Node* n1, const Node* n2, const BoundaryElement* e):
        n1( (n1->id < n2->id) ? n1 : n2), 
        n2( (n1->id < n2->id) ? n2 : n1),
        e(e)
    {}

    bool operator<(const Edge& e) const{
        if(n1->id < e.n1->id){
            return true;
        } else if(n1->id == e.n1->id){
            return n2->id < e.n2->id;
        }
        return false;
    }

    const Node* n1;
    const Node* n2;
    const BoundaryElement* e;
};

void InternalLoads2::initialize_views(Visualization* v){
        this->force_view = v->add_view("Internal Loads", spview::defs::ViewType::VECTOR, spview::defs::DataType::OTHER);
}

void InternalLoads2::initialize(const Meshing* const mesh, size_t& max_dofs){
    size_t bound_nodes_per_elem = mesh->elem_info->get_boundary_nodes_per_element();

    const auto& boundary_elements = mesh->boundary_elements;

    this->elem_info = mesh->elem_info;
    this->BOUNDARY_NODE_NUM = elem_info->get_boundary_nodes_per_element();
    const size_t node_num = elem_info->get_nodes_per_element();
    const size_t dof = elem_info->get_dof_per_node();
    const size_t s_size = elem_info->get_D_dimension();
    const size_t dim = (this->type == utils::PROBLEM_TYPE_2D) ? 2 : 3;
    const auto elem_type = elem_info->get_shape_type();
    const size_t vert_node_num = (elem_type == Element::Shape::TRI) ? 3 : 4;
    const size_t order = elem_info->get_element_order();

    this->u_dof = mesh->node_list.size()*dof;

    this->max_constr_num = s_size;
    if(this->use_M){
        this->il_constr_num = this->max_constr_num;
    } else {
        this->il_constr_num = dim;
    }

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
        this->prism_mesh.resize(num_elems);
        std::set<const MeshNode*> bnodes;
        for(size_t i = 0; i < boundary_elements.size(); ++i){
            if(apply_spring[i]){
                for(size_t j = 0; j < BOUNDARY_NODE_NUM; ++j){
                    bnodes.insert(boundary_elements[i].nodes[j]);
                }
            }
        }
        this->prism_nodes.reserve(num_elems);
        const auto RT = this->rot3D.T();
        std::map<size_t, MeshNode*> bound_node_mapping;
        for(auto& n:bnodes){
            gp_Pnt p;
            for(size_t i = 0; i < 3; ++i){
                double xi = 0;
                for(size_t j = 0; j < 3; ++j){
                    xi += RT(i,j)*n->point.Coord(1+j);
                }
                p.SetCoord(1+i, xi);
            }
            this->prism_nodes.push_back(std::make_unique<MeshNode>(p, n->id, 0));
            bound_node_mapping[n->id] = this->prism_nodes.back().get();
        }
        size_t cur_elem = 0;
        
        auto boundary_elem_info = mesh->elem_info->get_boundary_element_info();
        this->prism_elem_info = mesh->elem_info->get_prism_element_info();
        std::vector<MeshNode*> new_nodes(BOUNDARY_NODE_NUM);
        for(size_t i = 0; i < boundary_elements.size(); ++i){
            if(apply_spring[i]){
                this->submesh[cur_elem] = &boundary_elements[i];
                for(size_t n = 0; n < node_num; ++n){
                    this->F_id_map[boundary_elements[i].parent->nodes[n]->id] = -1;
                }
                for(size_t n = 0; n < BOUNDARY_NODE_NUM; ++n){
                    this->F_id_map[boundary_elements[i].nodes[n]->id] = 0;
                    new_nodes[n] = bound_node_mapping[boundary_elements[i].nodes[n]->id];
                }
                this->prism_mesh[cur_elem].reset(prism_elem_info->make_element(new_nodes, boundary_elements[i].parent));

                ++cur_elem;
            }
        }

        std::map<Edge, size_t> edge_count;
        for(const auto& e:this->submesh){
            for(size_t i = 0; i < vert_node_num; ++i){
                const size_t j = (i + 1) % BOUNDARY_NODE_NUM;
                const Edge edge(e->nodes[i], e->nodes[j], e);
                auto it = edge_count.find(edge);
                if(it != edge_count.end()){
                    it->second += 1;
                } else {
                    edge_count[edge] = 1;
                }
            }
        }

        std::set<const Node*> line_nodes;
        size_t bound_count = 0;
        for(const auto& e:edge_count){
            if(e.second == 1){
                line_nodes.insert(e.first.n1);
                line_nodes.insert(e.first.n2);
                ++bound_count;
            }
        }
        if(order > 1){
            for(const auto& e:edge_count){
                if(e.second == 1){
                    const auto p1 = e.first.n1->point;
                    const auto p2 = e.first.n2->point;
                    const gp_Vec v1(p1, p2);
                    for(size_t i = vert_node_num; i < BOUNDARY_NODE_NUM; ++i){
                        const auto p3 = e.first.e->nodes[i]->point;
                        const gp_Vec v2(p1, p3);
                        if((v1.Crossed(v2)).Magnitude() < 1e-10){
                            line_nodes.insert(e.first.e->nodes[i]);
                        }
                    }
                }
            }
        }
        this->edge_mesh.resize(bound_count);
        size_t ei = 0;
        for(const auto& e:edge_count){
            if(e.second == 1){
                auto& emi = edge_mesh[ei];
                emi.bounds.resize(2);
                emi.bounds[0] = e.first.n1->point;
                emi.bounds[1] = e.first.n2->point;
                emi.parent = e.first.e;

                gp_Vec v1(e.first.e->nodes[1]->point, e.first.e->nodes[0]->point);
                gp_Vec v2(e.first.e->nodes[2]->point, e.first.e->nodes[0]->point);
                gp_Ax1 ax(e.first.e->nodes[0]->point, v1.Crossed(v2));
                gp_Vec d(emi.bounds[1], emi.bounds[0]);
                d.Rotate(ax, M_PI/2.0);

                gp_Pnt bc = e.first.e->get_centroid(this->BOUNDARY_NODE_NUM);
                gp_Pnt lc = emi.bounds[0];
                lc.BaryCenter(1, emi.bounds[1], 1);
                gp_Vec bclc(bc, lc);

                double dir = bclc.Dot(d);
                if(dir > 0){
                    d.Reverse();
                }

                math::Vector dv({d.X(), d.Y(), d.Z()});
                dv.normalize();
                emi.transformed_normal = this->rot3D.T()*dv;

                ++ei;
            }
        }
        //logger::quick_log("edge mesh size", edge_mesh.size());


        size_t cur_id = 0;
        for(auto& n:F_id_map){
            if(n.second == 0){
                n.second = cur_id;
                ++cur_id;
            }
        }
        this->node_pos.resize(cur_id*this->f_node_dof, 0);

        //for(auto& n:line_nodes){
        //    const auto id = F_id_map[n->id];
        //    if(id >= 0){

        //        // psi
        //        this->node_pos[id*this->f_node_dof + 4] = -1;
        //        this->node_pos[id*this->f_node_dof + 5] = -1;
        //    }
        //}
        //for(const auto& e:edge_count){
        //    if(e.second == 1){
        //        const auto po1 = e.first.n1->point;
        //        const auto po2 = e.first.n2->point;
        //        math::Vector v1({po1.X(), po1.Y(), po1.Z()});
        //        math::Vector v2({po2.X(), po2.Y(), po2.Z()});
        //        math::Vector p1(rot3D.T()*v1);
        //        math::Vector p2(rot3D.T()*v2);
        //        const auto id1 = F_id_map[e.first.n1->id];
        //        const auto id2 = F_id_map[e.first.n2->id];
        //        // tau_1
        //        if(std::abs(p1[1] - p2[1]) > 1e-10){
        //            this->node_pos[id1*this->f_node_dof + 4] = -1;
        //            this->node_pos[id2*this->f_node_dof + 4] = -1;
        //        }
        //        // tau_2
        //        if(std::abs(p1[0] - p2[0]) > 1e-10){
        //            this->node_pos[id1*this->f_node_dof + 5] = -1;
        //            this->node_pos[id2*this->f_node_dof + 5] = -1;
        //        }
        //        //this->node_pos[id1*this->f_node_dof + 0] = -1;
        //        //this->node_pos[id2*this->f_node_dof + 0] = -1;
        //        //this->node_pos[id1*this->f_node_dof + 1] = -1;
        //        //this->node_pos[id2*this->f_node_dof + 1] = -1;
        //        //this->node_pos[id1*this->f_node_dof + 2] = -1;
        //        //this->node_pos[id2*this->f_node_dof + 2] = -1;
        //        //this->node_pos[id1*this->f_node_dof + 3] = -1;
        //        //this->node_pos[id2*this->f_node_dof + 3] = -1;
        //        //this->node_pos[id1*this->f_node_dof + 4] = -1;
        //        //this->node_pos[id2*this->f_node_dof + 4] = -1;
        //        //this->node_pos[id1*this->f_node_dof + 5] = -1;
        //        //this->node_pos[id2*this->f_node_dof + 5] = -1;
        //    }
        //}

        //for(auto& n:F_id_map){
        //    const auto id = n.second;
        //    if(id >= 0){
        //        this->node_pos[id*this->f_node_dof + 0] = -1;
        //        this->node_pos[id*this->f_node_dof + 1] = -1;
        //        this->node_pos[id*this->f_node_dof + 3] = -1;
        //    }
        //}
        long dof_id = 0;
        for(size_t i = 0; i < this->node_pos.size(); ++i){
            auto& pi = this->node_pos[i];
            if(pi > -1){
                pi = dof_id + max_dofs;
                ++dof_id;
            }
        }

        for(const auto& l:edge_mesh){
            for(size_t n = 0; n < node_num; ++n){
                this->line_id_map[l.parent->parent->nodes[n]->id] = -1;
            }
        }
        for(const auto& n:line_nodes){
            this->line_id_map[n->id] = 0;
        }
        size_t line_id = 0;
        for(auto& n:line_id_map){
            if(n.second == 0){
                n.second = line_id;
                ++line_id;
            }
        }


        this->il2_offset = max_dofs;
        //this->F_total_dof = cur_id*dim;
        this->F_total_dof = dof_id;
        //this->xc_offset = this->ds_offset + cur_id*dim;
        
        //this->ds_offset = this->il2_offset + this->F_total_dof;
        //this->constraint_offset = this->ds_offset + cur_id*dim;

        //this->ds_offset = this->il2_offset + this->F_total_dof;
        //this->constraint_offset = this->ds_offset + this->F_total_dof;
                                                  
        //this->bc_offset = this->il2_offset + this->F_total_dof;
        //this->constraint_offset = this->bc_offset + line_id*dim;

        //this->f2_offset = this->il2_offset + this->F_total_dof;
        //this->constraint_offset = this->f2_offset + this->F_total_dof;

        this->constraint_offset = this->il2_offset;// + 3*this->F_total_dof;

        //this->xc_offset = this->constraint_offset + il_constr_num;

        //if(use_M){
        //    this->n_offset = this->xc_offset + dim;
        //} else {
        //    this->n_offset = this->xc_offset;
        //}
        this->xc_offset = this->constraint_offset + 3*this->F_total_dof;// + il_constr_num;
        this->n_offset = this->xc_offset;
    }
    //if(use_M){
    //    max_dofs = this->n_offset + 1;//2;
    //} else {
    //    max_dofs = this->n_offset;// + 1;
    //}
    max_dofs = this->n_offset;// + 1;

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

    logger::quick_log("A", A);
    logger::quick_log("EA", EA);
    logger::quick_log("EA_u", EA_u);
    logger::quick_log("EA_v", EA_v);
    logger::quick_log("EA_w", EA_w);
    logger::quick_log(this->center);

    {
        // 0 to 1
        math::Matrix C({
             1,  1,  1,
             1,  0,  0,
             1,  2,  3}, 3, 3);

        this->V = math::Matrix(3, 3, 0);
        math::LU Cl(C);
        math::Vector r(3);

        const auto dp = this->rot3D.T()*(this->P_prev - this->center);
        const math::Vector nv({
                this->normal.X(),
                this->normal.Y(),
                this->normal.Z()
                });
        const auto Rnv = this->rot3D.T()*nv;
        for(size_t i = 0; i < dim; ++i){
            r[0] = dp[i];
            r[1] = 2.5*std::abs(r[0])*Rnv[i];
            r[2] = 2.5*std::abs(r[0])*Rnv[i];

            Cl.solve(r);

            for(size_t j = 0; j < 3; ++j){
                V(i,j) = r[j];
            }
        }
    }
    logger::quick_log(V);
    {
        // 0 to 1
        auto Nzt = math::Matrix(4, 4, 0);
        for(size_t i = 0; i < 4; ++i){
            const double t = i/3.0;
            for(size_t j = 0; j < 4; ++j){
                Nzt(i,j) = std::pow(t, j);
            }
        }
        Nzt.invert_LU();
        this->Nz = Nzt.T();
    }

    double pV = 0;
    for(const auto& p:this->prism_mesh){
        pV += p->prism_volume(this->V);
    }
    logger::quick_log("PRISM VOLUME", pV);

    if(this->use_M){
        this->D = 0;
        //for(size_t i = 0; i < dim; ++i){
        //    this->D += this->center[i]*this->normal.Coord(1+i);
        //}

        //math::Matrix Mc_tmp(3, 3, 0);

        //const math::Matrix M_conv(
        //        {1, 0, 0,
        //         0, 1, 0}, 2, 3);
        //const auto R_conv = M_conv*(this->rot3D.T());

        //// Contribution from the integrals
        //for(size_t i = 0; i < 3; ++i){
        //    Mc_tmp(0,i) = -F[1]*R_conv(0,i) - F[2]*R_conv(1,i);
        //    Mc_tmp(1,i) = -F[0]*R_conv(0,i);
        //    Mc_tmp(2,i) = -F[0]*R_conv(1,i);
        //}


        math::Matrix Mc_tmp(
                {-F[1],  F[0], 0,
                     0, -F[2], 0,
                 -F[2],     0, 0}, 3, 3);

        // here rot3D transforms x_c, not F
        Mc_tmp = Mc_tmp*this->rot3D.T();

        // The second line does need to have inverted signs due to get_loads()
        const math::Matrix M_skew(
            {
                      0, -F_orig[2],  F_orig[1],
             -F_orig[2],          0,  F_orig[0],
             -F_orig[1],  F_orig[0],     0
            }, 3, 3);

        //logger::quick_log(Mc_tmp);
        //logger::quick_log(this->rot3D.T()*M_skew);
        //exit(0);
        

        // Calculated using the global coordinate system
        //if(this->calculate_from_last_point){
        //    if(!this->inverted_calc){
        //        //this->M = -((this->center - this->P_prev).cross(this->F) + this->M_orig);
        //        this->M = this->rot3D.T()*(this->M_orig - (-this->P_prev).cross(this->F_orig));
        //        // This is (apparently) correct due to the signs used in M_skew
        //        // and in get_loads()
        //        Mc_tmp += this->rot3D.T()*M_skew;
        //    } else {
        //        //this->M = ((this->center - this->P_prev).cross(this->F) + this->M_orig);
        //        this->M = this->rot3D.T()*(-(this->M_orig - (-this->P_prev).cross(this->F_orig)));
        //        // This is (apparently) correct due to the signs used in M_skew
        //        // and in get_loads()
        //        Mc_tmp -= this->rot3D.T()*M_skew;
        //    }
        //} else {
        //    this->M = this->rot3D.T()*this->M_orig;
        //}
        this->M = this->rot3D.T()*this->M_orig;

        //this->center = this->rot3D.T()*this->center;

        //this->Mc = std::move(Mc_tmp);

        this->assume_xc.resize(dim, true);
        //this->assume_xc.resize(dim, false);
        //for(size_t i = 0; i < dim; ++i){
        //    bool all_zero = true;
        //    for(size_t j = 0; j < dim; ++j){
        //        if(std::abs(this->Mc(j,i)) > 1e-10){
        //            all_zero = false;
        //            break;
        //        }
        //    }
        //    if(all_zero && std::abs(this->normal.Coord(1+i)) < 1e-10){
        //        this->assume_xc[i] = true;
        //        --max_dofs;
        //        //--this->constraint_offset;
        //        --this->n_offset;
        //    }
        //}
    }

    auto Mv = this->integrate_vv_mat(mesh);
    math::LU Mlu(Mv);
    this->vv = this->get_loads();
    Mlu.solve(this->vv);

    logger::quick_log("il2_offset", il2_offset);
    //logger::quick_log("bc_offset", bc_offset);
    logger::quick_log("constraint_offset", constraint_offset);
    logger::quick_log("xc_offset", xc_offset);
    //logger::quick_log("n_offset", n_offset);
    //logger::quick_log("assume_xc");
    //logger::quick_log(assume_xc);

    logger::quick_log("A: ", A);
    logger::quick_log("EA: ", EA);
    logger::quick_log("EA_u", EA_u);
    logger::quick_log("EA_v", EA_v);
    logger::quick_log("EA_w", EA_w);
    logger::quick_log("Material center");
    logger::quick_log(this->center);
    logger::quick_log("F");
    logger::quick_log(this->F);
    if(this->use_M){
        logger::quick_log("M");
        logger::quick_log(this->M);

        //logger::quick_log("Neutral axis matrix");
        //logger::quick_log(Mc);
    }

}

void InternalLoads2::display_views(const std::vector<double>& u) const{
    std::vector<double> uf(this->u_dof, 0); 

    const size_t dof      = elem_info->get_dof_per_node();
    const size_t node_num = this->prism_elem_info->get_nodes_per_element();

    //for(const auto& b:submesh){
    //    std::vector<long> u_pos(dof*node_num);
    //    std::vector<long> f_pos(f_node_dof*node_num);
    //    for(size_t i = 0; i < node_num; ++i){
    //        const auto& n = b->parent->nodes[i];
    //        const auto f_id = F_id_map.at(n->id);
    //        for(size_t j = 0; j < dof; ++j){
    //            const size_t p = n->id*dof + j;
    //            u_pos[i*dof + j] = p;
    //        }
    //        if(f_id >= 0){
    //            for(size_t j = 0; j < f_node_dof; ++j){
    //                const size_t fp = f_id*f_node_dof + j;
    //                f_pos[i*f_node_dof + j] = node_pos[fp];
    //            }
    //        } else {
    //            for(size_t j = 0; j < f_node_dof; ++j){
    //                f_pos[i*f_node_dof + j] = -1;
    //            }
    //        }
    //    }
    //    const auto F = get_F_mat(b);
    //    
    //    for(size_t i = 0; i < u_pos.size(); ++i){
    //        for(size_t j = 0; j < f_pos.size(); ++j){
    //            if(f_pos[j] > -1){
    //                uf[u_pos[i]] -= F(i,j)*u[f_pos[j]];
    //            }
    //        }
    //    }
    //}
    
    /*
    double energy = 0;
    for(const auto& b:submesh){
        math::Vector f_pos(f_node_dof*node_num);
        for(size_t i = 0; i < node_num; ++i){
            const auto& n = b->parent->nodes[i];
            const auto f_id = F_id_map.at(n->id);
            if(f_id >= 0){
                for(size_t j = 0; j < f_node_dof; ++j){
                    const size_t fp = f_id*f_node_dof + j;
                    f_pos[i*f_node_dof + j] = u[node_pos[fp]];
                }
            } else {
                for(size_t j = 0; j < f_node_dof; ++j){
                    f_pos[i*f_node_dof + j] = -1;
                }
            }
        }
        const auto F = get_K_mat(b);

        energy += 0.5*(f_pos.T()*F*f_pos);
    }
    logger::quick_log("IL2 strain energy", energy);
    
    const auto& n = this->normal;
    const math::Matrix M({
            n.X(), 0, 0, n.Y(), n.Z(), 0,
            0, n.Y(), 0, n.X(), 0, n.Z(),
            0, 0, n.Z(), 0, n.X(), n.Y()}, 3, 6);
    for(const auto& b:submesh){
        math::Vector f_pos(f_node_dof);
        for(size_t i = 0; i < this->BOUNDARY_NODE_NUM; ++i){
            const auto& n = b->nodes[i];
            const auto f_id = F_id_map.at(n->id);
            if(f_id >= 0){
                for(size_t j = 0; j < f_node_dof; ++j){
                    const size_t fp = f_id*f_node_dof + j;
                    f_pos[j] = u[node_pos[fp]];
                }
                const auto fi = M*this->rot3D_D.T()*f_pos;
                for(size_t j = 0; j < dof; ++j){
                    const size_t p = n->id*dof + j;
                    uf[p] = fi[j];
                }
            }
        }
    }

    this->force_view->update_view(uf);
    */
    for(size_t ei = 0; ei < this->submesh.size(); ++ei){
        const auto& b = this->submesh[ei];
        const auto& p = this->prism_mesh[ei];
        std::vector<long> u_pos(dof*node_num);
        std::vector<long> f_pos(dof*node_num);
        for(size_t i = 0; i < node_num; ++i){
            const auto& n = b->nodes[i];
            const auto f_id = F_id_map.at(n->id);
            for(size_t j = 0; j < dof; ++j){
                const size_t p = n->id*dof + j;
                u_pos[i*dof + j] = p;
            }
            if(f_id >= 0){
                for(size_t j = 0; j < f_node_dof; ++j){
                    const size_t fp = f_id*f_node_dof + j;
                    f_pos[i*f_node_dof + j] = node_pos[fp] + 2*this->F_total_dof;
                }
            } else {
                for(size_t j = 0; j < f_node_dof; ++j){
                    f_pos[i*f_node_dof + j] = -1;
                }
            }
        }
        //const auto F = get_F_mat(b);
        const auto FE = this->get_F2(p.get());
        
        for(size_t i = 0; i < u_pos.size(); ++i){
            if(f_pos[i] > -1){
                uf[u_pos[i]] += FE[i];
            }
        }
    }
    this->force_view->update_view(uf);
}

void InternalLoads2::apply_constraints(std::vector<double>& load) const{
    //if(this->type == utils::PROBLEM_TYPE_2D){

    //} else {
    //    math::Vector V = this->get_loads();
    //    for(size_t i = 0; i < this->il_constr_num; ++i){
    //        load[this->constraint_offset + i] = V[i];
    //        //load[this->constraint_offset + 6 + i] = V[i];
    //    }
    //    //if(this->use_M){
    //    //    load[this->n_offset] = this->D;
    //    //}
    //}

    
    const size_t dof      = this->elem_info->get_dof_per_node();
    const size_t node_num = this->prism_elem_info->get_nodes_per_element();
    std::vector<long> f3_pos(dof*node_num);
    math::Vector fsum(3, 0);
    for(size_t ei = 0; ei < this->submesh.size(); ++ei){
        //const auto& b = this->submesh[ei];
        const auto& p = this->prism_mesh[ei];
        for(size_t i = 0; i < node_num; ++i){
            const auto& n = p->nodes[i];
            const auto f_id = this->F_id_map.at(n->id);
            if(f_id >= 0){
                for(size_t j = 0; j < dof; ++j){
                    const size_t fp = f_id*dof + j;
                    f3_pos[i*this->f_node_dof + j] = this->node_pos[fp] + 2*this->F_total_dof;
                }
            } else {
                for(size_t j = 0; j < dof; ++j){
                    f3_pos[i*this->f_node_dof + j] = -1;
                }
            }
        }

        const auto FE = this->get_F2(p.get());
        for(size_t i = 0; i < f3_pos.size(); ++i){
            if(f3_pos[i] >= 0){
                load[f3_pos[i]] += FE[i];
                fsum( i % 3 ) += FE[i];
            }
        }
    }
    logger::quick_log("APPLIED FORCE SUM");
    logger::quick_log(fsum);
}

math::Matrix InternalLoads2::get_F_mat(const BoundaryElement* b) const{
    std::vector<gp_Pnt> bounds(BOUNDARY_NODE_NUM);
    for(size_t i = 0; i < BOUNDARY_NODE_NUM; ++i){
        bounds[i] = b->nodes[i]->point;
    }

    if(this->type == utils::PROBLEM_TYPE_2D){
        //auto NN = -b->parent->get_R(this->rot2D, this->thickness, bounds);
        auto NN = -b->parent->il2_R2(this->rot2D, this->rot3D_D, this->normal, this->thickness, bounds);

        return NN; 
    } else {
        //auto NN = -b->parent->get_R(this->rot3D, this->thickness, bounds);
        // Use of rot3D_S follows TCT Ting
        auto NN = -b->parent->il2_R2(this->rot3D, this->rot3D_S, this->normal, this->thickness, bounds);

        return NN; 
    }
}

math::Matrix InternalLoads2::get_K_mat(const BoundaryElement* b) const{
    std::vector<gp_Pnt> bounds(BOUNDARY_NODE_NUM);
    for(size_t i = 0; i < BOUNDARY_NODE_NUM; ++i){
        bounds[i] = b->nodes[i]->point;
    }

    if(this->type == utils::PROBLEM_TYPE_2D){
        auto NN = b->parent->il2_F(this->rot2D, this->rot2D_S, this->mmat, bounds, this->thickness);

        return NN; 
    } else {
        auto NN = b->parent->il2_F(this->rot3D, this->rot3D_S, this->mmat, bounds, this->thickness);
        //auto NN = b->parent->il2_F_alt(this->rot3D, bounds, this->thickness);

        return NN; 
    }
}
math::Matrix InternalLoads2::get_K2(const PrismMeshElement* b) const{
    std::vector<gp_Pnt> bounds(BOUNDARY_NODE_NUM);
    for(size_t i = 0; i < BOUNDARY_NODE_NUM; ++i){
        bounds[i] = b->nodes[i]->point;
    }

    if(this->type == utils::PROBLEM_TYPE_2D){
        //auto NN = b->parent->il2_F(this->rot2D, this->rot2D_S, this->mmat, bounds, this->thickness);

        //return NN; 
        return math::Matrix();
    } else {
        auto NN = b->il2_k_ext(this->rot3D, this->rot3D_D, this->mmat, this->V, this->Nz);

        return NN; 
    }
}

math::Vector InternalLoads2::get_F2(const PrismMeshElement* b) const{
    std::vector<gp_Pnt> bounds(BOUNDARY_NODE_NUM);
    for(size_t i = 0; i < BOUNDARY_NODE_NUM; ++i){
        bounds[i] = b->nodes[i]->point;
    }

    if(this->type == utils::PROBLEM_TYPE_2D){
        //auto NN = b->parent->il2_F(this->rot2D, this->rot2D_S, this->mmat, bounds, this->thickness);

        //return NN; 

        return {};
    } else {
        auto NN = b->il2_f_ext(this->center, this->normal, this->rot3D, this->rot3D_D, this->mmat, this->vv);
        //auto NN = b->parent->il2_F_alt(this->rot3D, bounds, this->thickness);

        return NN; 
    }
}

math::Matrix InternalLoads2::get_bc(const EdgeElement* b) const{

    if(this->type == utils::PROBLEM_TYPE_2D){
        auto NN = b->parent->parent->il2_bc(this->rot2D, b->transformed_normal, b->bounds);

        return NN; 
    } else {
        //auto NN = b->parent->parent->il2_bc(this->rot3D, b->transformed_normal, b->bounds);
        auto NN = b->parent->parent->il2_bc2(this->rot3D, b->transformed_normal, b->bounds);

        return NN; 
    }
}

math::Matrix InternalLoads2::get_t1t2(const BoundaryElement* b) const{
    std::vector<gp_Pnt> bounds(BOUNDARY_NODE_NUM);
    for(size_t i = 0; i < BOUNDARY_NODE_NUM; ++i){
        bounds[i] = b->nodes[i]->point;
    }

    if(this->type == utils::PROBLEM_TYPE_2D){
        auto NN = b->parent->il2_t1t2(this->rot2D, thickness, bounds);

        return NN; 
    } else {
        auto NN = b->parent->il2_t1t2(this->rot3D, thickness, bounds);

        return NN; 
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
            result += e->parent->calculate_internal_loads(this->rot2D, this->rot2D_D, this->mmat, u_ext, bounds);
        }
    } else {
        for(const auto& e:this->submesh){
            for(size_t i = 0; i < bnum; ++i){
                bounds[i] = e->nodes[i]->point;
            }
            result += e->parent->calculate_internal_loads(this->rot3D, this->rot3D_D, this->mmat, u_ext, bounds);
        }
    }
    //-m[1],
    // m[0],
    // f[2],
    // m[2],
    // f[0],
    // f[1],
    logger::quick_log("IL2 resulting F");
    logger::quick_log(-result[2], result[5], result[4]);
    logger::quick_log("IL2 resulting M");
    logger::quick_log(-result[3], result[0], result[1]);
    logger::quick_log("(there are approximate results, seeing as they disagree with reactions)");


    return result;
}

void InternalLoads2::check_results(const std::vector<double>& u) const{
    if(this->use_M){
        (void)u;
        //const size_t node_num = elem_info->get_nodes_per_element();
        //const size_t s_size = elem_info->get_D_dimension();
        //const size_t dim = (this->type == utils::PROBLEM_TYPE_2D) ? 2 : 3;

        //math::Vector xc(dim);
        //for(size_t i = 0; i < dim; ++i){
        //    if(!this->assume_xc[i]){
        //        xc[i] = u[this->xc_offset + i];
        //    } else {
        //        logger::quick_log("IL2: assuming neutral axis equal to material center for coord", i);
        //        xc[i] = this->center[i];
        //    }
        //}
        //logger::quick_log("IL2 neutral axis");
        //logger::quick_log(xc);
        if(this->calculate_from_last_point){
            logger::quick_log("IL2 resulting moments");
            if(!this->inverted_calc){
                //logger::quick_log((this->M_orig - (xc - this->P_prev).cross(this->F_orig)));
                logger::quick_log((this->M_orig - (- this->P_prev).cross(this->F_orig)));
            } else {
                logger::quick_log(-(this->M_orig - (- this->P_prev).cross(this->F_orig)));
                //logger::quick_log(-(this->M_orig - (xc - this->P_prev).cross(this->F_orig)));
            }
        }
    }
}

math::Vector InternalLoads2::calculate_constraints(const std::vector<double>& u_ext, bool print) const{
    auto result = this->calculate_internal_loads(u_ext);

    if(print){
        logger::quick_log("Internal loads:");
        logger::quick_log("My Mz Fx Mx Fz Fy");
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

    return math::Vector(6);
}
math::Vector InternalLoads2::get_loads() const{

    const auto f = this->F;
    const auto m = this->M;

    if(this->use_M){
        return {
            -m[1],
             m[0],
             f[2],
             m[2],
             f[0],
             f[1],
            //    0,
            //-f[2],
            // f[1]
        };
    } else {
        return {
             f[2],
             f[0],
             f[1]
        };
    }
}

math::Matrix InternalLoads2::get_constr(const BoundaryElement* b) const{
    std::vector<gp_Pnt> bounds(BOUNDARY_NODE_NUM);
    for(size_t i = 0; i < BOUNDARY_NODE_NUM; ++i){
        bounds[i] = b->nodes[i]->point;
    }
    return b->parent->constraint_matrix_alt(this->center, this->rot3D, bounds, !this->calculate_from_last_point, this->thickness);
}
math::Matrix InternalLoads2::get_constr2(const BoundaryElement* b) const{
    std::vector<gp_Pnt> bounds(BOUNDARY_NODE_NUM);
    for(size_t i = 0; i < BOUNDARY_NODE_NUM; ++i){
        bounds[i] = b->nodes[i]->point;
    }
    const auto D = this->mat->stiffness_3D(b->parent, b->parent->get_centroid());
    return b->parent->constraint_matrix(this->center, this->rot3D, this->rot3D_S, this->mmat, bounds, this->calculate_from_last_point);
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
math::Matrix InternalLoads2::integrate_vv_mat(const Meshing* const mesh) const{
    math::Matrix result(6, 6, 0);
    if(mesh->elem_info->get_shape_type() == Element::Shape::TRI){
        constexpr size_t N = 3;
        {
            math::Matrix result_tmp(6, 6, 0);
            for(size_t i = 0; i < submesh.size(); ++i){
                const auto& e = submesh[i];
                std::array<gp_Pnt, N> abs_points;
                for(size_t x = 0; x < N; ++x){
                    abs_points[x] = e->nodes[x]->point;
                }
                this->GS_vv_mat(e->parent, abs_points, result_tmp);
                result += result_tmp;
                result_tmp.fill(0);
            }
        }
    } else if(mesh->elem_info->get_shape_type() == Element::Shape::QUAD){
        //constexpr size_t N = 4;
        //#pragma omp parallel
        //{
        //    VectorWrapper result_tmp(fn.size());
        //    #pragma omp for reduction(vecsum:result)
        //    for(size_t i = 0; i < submesh.size(); ++i){
        //        const auto& e = submesh[i];
        //        std::array<gp_Pnt, N> abs_points;
        //        for(size_t x = 0; x < N; ++x){
        //            abs_points[x] = e->nodes[x]->point;
        //        }
        //        this->GS_quad(e->parent, abs_points, fn, result_tmp.vec);
        //        result += result_tmp;
        //    }
        //}
    }
    return result;
}

void InternalLoads2::GS_vv_mat(const MeshElement* const e, const std::array<gp_Pnt, 3>& p, math::Matrix& result) const{
    gp_Vec v1(p[1], p[0]);
    gp_Vec v2(p[2], p[0]);
    result.fill(0);
    math::Matrix D;
    const double drnorm = (v1.Crossed(v2)).Magnitude()/2;
    const auto& gsi = utils::GaussLegendreTri<MAX_GN>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        gp_Pnt pxi(this->GLT_point(*it, p));
        if(this->type == utils::PROBLEM_TYPE_2D){
            D = this->mat->stiffness_2D(e, pxi);
            this->mat->rotate_D(D, this->rot2D_D);
        } else {
            D = this->mat->stiffness_3D(e, pxi);
            this->mat->rotate_D(D, this->rot3D_D);
        }
        result += it->w*this->make_vv_mat_3D(D, pxi);
    }
    result *= drnorm;
}

void InternalLoads2::GS_tri(const MeshElement* const e, const std::array<gp_Pnt, 3>& p, const std::vector<std::function<double(const math::Matrix& S, const gp_Pnt& px)>>& fn, math::Vector& result) const{
    gp_Vec v1(p[1], p[0]);
    gp_Vec v2(p[2], p[0]);
    result.fill(0);
    math::Matrix S;
    const double drnorm = (v1.Crossed(v2)).Magnitude()/2;
    const auto& gsi = utils::GaussLegendreTri<MAX_GN>::get();
    for(auto it = gsi.begin(); it < gsi.end(); ++it){
        for(size_t i = 0; i < fn.size(); ++i){
            gp_Pnt pxi(this->GLT_point(*it, p));
            if(this->type == utils::PROBLEM_TYPE_2D){
                S = this->mat->stiffness_inverse_2D(e, pxi);
                this->mat->rotate_S(S, this->rot2D_S);
            } else {
                S = this->mat->stiffness_inverse_3D(e, pxi);
                this->mat->rotate_S(S, this->rot3D_S);
            }
            result[i] += it->w*fn[i](S, pxi);
        }
    }
    result *= drnorm;
}

void InternalLoads2::GS_quad(const MeshElement* const e, const std::array<gp_Pnt, 4>& p, const std::vector<std::function<double(const math::Matrix& S, const gp_Pnt& px)>>& fn, math::Vector& result) const{
    result.fill(0);
    math::Matrix S;
    const auto& gsi = utils::GaussLegendre<MAX_GN>::get();
    for(auto xi = gsi.begin(); xi < gsi.end(); ++xi){
        for(auto eta = gsi.begin(); eta < gsi.end(); ++eta){
            for(size_t i = 0; i < fn.size(); ++i){
                gp_Pnt pxi(this->GL_point(*xi, *eta, p));
                if(this->type == utils::PROBLEM_TYPE_2D){
                    S = this->mat->stiffness_inverse_2D(e, pxi);
                    this->mat->rotate_S(S, this->rot2D_S);
                } else {
                    S = this->mat->stiffness_inverse_3D(e, pxi);
                    this->mat->rotate_S(S, this->rot3D_S);
                }
                const double detJ = std::abs(this->J(xi->x, eta->x, p).determinant());
                result[i] += xi->w*eta->w*detJ*fn[i](S, pxi);
            }
        }
    }
}
