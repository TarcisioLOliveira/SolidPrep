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

#include "internal_loads.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include "meshing.hpp"
#include <memory>
#include <set>

inline Eigen::Matrix<double, 3, 3> Lek_rot3D(Eigen::Matrix<double, 3, 3> L, Eigen::Matrix<double, 3, 3> R){
    return L*R;
}

InternalLoads::InternalLoads(CrossSection cross_section, double thickness, gp_Dir normal, gp_Dir v, gp_Dir w, Material* mat, std::array<double, 3> F, std::array<double, 3> M, MeshElementFactory* elem, BoundaryMeshElementFactory* bound_elem, utils::ProblemType type):
    S(std::move(cross_section)), 
    A(S.get_area()),
    thickness(thickness),
    rot2D{{normal.X(), v.X()},
          {normal.Y(), v.Y()}},
    Lek_basis
        {{0, 0, -1},
         {0, 1, 0},
         {1, 0, 0}},
    rot3D(Lek_rot3D(Lek_basis,
        Eigen::Matrix<double, 3, 3>
        {{normal.X(), v.X(), w.X()},
         {normal.Y(), v.Y(), w.Y()},
         {normal.Z(), v.Z(), w.Z()}})),
    F(F), M(M),
    mat(mat),
    normal(normal),
    elem_info(elem),
    boundary_elem_info(),
    v(v), w(w),
    type(type){

    this->boundary_elem_info = bound_elem;

    this->curvature = std::make_unique<Curvature>(mat, rot2D, rot3D, this->boundary_elem_info, F[0], F[1], F[2], M[0], M[1], M[2]);
}

void InternalLoads::calculate_curvature(std::vector<BoundaryElement>& boundary_elements){
    this->generate_mesh(boundary_elements);
    this->curvature->generate_curvature_3D(this->boundary_nodes, this->boundary_mesh, this->phi_size, this->boundary_nodes.size());

}

void InternalLoads::apply_load_2D(const std::vector<long>& node_positions, std::vector<double>& load_vector) const{
    auto is_between_points = [](gp_Pnt p1, gp_Pnt p2, gp_Pnt p)->bool{
        gp_Mat M(1, p1.X(), p1.Y(), 1, p2.X(), p2.Y(), 1, p.X(), p.Y());
        bool in_line = std::abs(M.Determinant()) < Precision::Confusion();
        bool within_bounds = p.Distance(p1) - p1.Distance(p2) < Precision::Confusion() &&
                             p.Distance(p2) - p1.Distance(p2) < Precision::Confusion();

        return in_line && within_bounds;
    };

    const size_t bound_nodes_per_elem = this->elem_info->get_boundary_nodes_per_element();
    const size_t nodes_per_elem = this->elem_info->get_nodes_per_element();
    const size_t dof = this->elem_info->get_dof_per_node();

    const size_t N = 2;
    Eigen::Vector<double, N> F0{this->F[0]/A, this->F[1]/A};
    Eigen::Vector<double, N> Fr = this->rot2D*F0;
    std::vector<double> F{Fr[0], Fr[1]};
    Eigen::Matrix<double, N, N> S{{0, 0},
                                  {0, 0}};
    std::vector<double> Sn(N*N);

    gp_Dir Snormal = this->S.get_normal();
    double Ssize = this->S.get_dimension();
    gp_Dir line_dir = Snormal.Rotated(gp_Ax1(center, gp_Dir(0,0,1)), M_PI/2);
    gp_Pnt p1 = center.Translated( 0.5*Ssize*line_dir);
    gp_Pnt p2 = center.Translated(-0.5*Ssize*line_dir);
    for(size_t j = 0; j < this->submesh.size(); ++j){
        const auto& e = this->submesh[j];
        const gp_Pnt c = e->get_centroid(bound_nodes_per_elem);
        const auto E = this->mat->beam_E_2D(e->parent, c, this->rot2D);
        S(0,1) = E*this->curv[0];
        Eigen::Matrix<double, N, N> Srot = this->rot2D*S*this->rot2D.transpose();
        for(size_t row = 0; row < N; ++row){
            for(size_t col = 0; col < N; ++col){
                Sn[row*N + col] = Srot(row,col);
            }
        }
        std::vector<gp_Pnt> list;
        list.reserve(bound_nodes_per_elem);
        for(size_t i = 0; i < bound_nodes_per_elem; ++i){
            auto n = e->nodes[i];
            if(is_between_points(p1, p2, n->point)){
                list.push_back(n->point);
            }
        }

        std::vector<double> Rf;
        if(list.size() >= 2){
             Rf = submesh[j]->parent->get_Rf(Sn, F, this->center, this->thickness, list);
        } else if(list.size() > 1) {
            for(size_t i = 0; i < bound_nodes_per_elem; ++i){
                size_t j = (i+1)%bound_nodes_per_elem;
                auto n1 = e->nodes[i]->point;
                auto n2 = e->nodes[j]->point;
                if(p1.IsEqual(n1, Precision::Confusion()) || p1.IsEqual(n2, Precision::Confusion()) ||
                   p2.IsEqual(n1, Precision::Confusion()) || p2.IsEqual(n2, Precision::Confusion())){
                    continue;
                }
                if(is_between_points(n1, n2, p1)){
                    list.push_back(p1);
                    break;
                } else if(is_between_points(n1, n2, p2)){
                    list.push_back(p2);
                    break;
                }
            }
            Rf = submesh[j]->parent->get_Rf(Sn, F, this->center, this->thickness, list);
        }
        if(Rf.size() > 0){
            //logger::quick_log(fe);
            //for(auto& ff:fe){
            //    F += ff;
            //}
            for(size_t i = 0; i < nodes_per_elem; ++i){
                for(size_t j = 0; j < dof; ++j){
                    auto& n = e->parent->nodes[i];
                    const auto p = node_positions[n->u_pos[j]];
                    if(p >= 0){
                        load_vector[p] += Rf[i*dof+j];
                    }
                }
            }
        }
    }
}
void InternalLoads::apply_load_3D(const std::vector<long>& node_positions, std::vector<double>& load_vector) const{
    const size_t bound_nodes_per_elem = this->elem_info->get_boundary_nodes_per_element();
    const size_t nodes_per_elem = this->elem_info->get_nodes_per_element();
    const size_t dof = this->elem_info->get_dof_per_node();

    Eigen::Vector<double, 3> FF{0, 0, 0};
    for(size_t j = 0; j < submesh.size(); ++j){
        const auto& e = submesh[j];
        const auto& b = this->boundary_mesh[j];
        std::vector<gp_Pnt> points(bound_nodes_per_elem);
        for(size_t i = 0; i < bound_nodes_per_elem; ++i){
            points[i] = e->nodes[i]->point;
        }

        const auto Rf = this->curvature->get_force_vector_3D(b.get());
        for(size_t i = 0; i < nodes_per_elem; ++i){
            for(size_t j = 0; j < dof; ++j){
                const auto n = e->parent->nodes[i];
                const auto p = node_positions[n->u_pos[j]];
                FF[j] += Rf[i*dof + j];
                if(p >= 0){
                    load_vector[p] += Rf[i*dof+j];
                }
            }
        }
    }
    logger::quick_log("Calculated:", FF[0], FF[1], FF[2]);
}

class PointSort{
    public:
    bool operator()(const std::unique_ptr<MeshNode>& n1, const std::unique_ptr<MeshNode>& n2) const{
        if(n1->point.X() < n2->point.X()){
            return true;
        } else if(n1->point.X() == n2->point.X() && n1->point.Y() < n2->point.Y()){
            return true;
        } else if(n1->point.X() == n2->point.X() && n1->point.Y() == n2->point.Y() && n1->point.Z() < n2->point.Z()){
            return true;
        }
        return false;
    }
};

struct NodeComp{
    NodeComp(gp_Pnt p): p{p} {}
    bool operator()(const std::unique_ptr<MeshNode>& o) const{ 
        return o->point.IsEqual(p, Precision::Confusion());
    }
    const gp_Pnt p;
};

void InternalLoads::generate_mesh(const std::vector<BoundaryElement>& boundary_elements){
    size_t bound_nodes_per_elem = this->elem_info->get_boundary_nodes_per_element();

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

    // Copy and transform necessary nodes
    std::set<std::unique_ptr<MeshNode>, PointSort> nodes;
    for(size_t i = 0; i < boundary_elements.size(); ++i){
        if(apply_spring[i]){
            const auto& b = boundary_elements[i];
            for(size_t j = 0; j < bound_nodes_per_elem; ++j){
                const auto& n = b.nodes[j];
                gp_Pnt p = utils::change_point(n->point, this->rot3D.transpose());
                std::unique_ptr<MeshNode> node(std::make_unique<MeshNode>(p, 0, 1));
                node->u_pos[0] = 0;
                nodes.insert(std::move(node));
            }
        }
    }

    this->boundary_nodes.resize(nodes.size());
    for(size_t i = 0; i < this->boundary_nodes.size(); ++i){
        auto nit = nodes.begin();
        this->boundary_nodes[i] = std::move(nodes.extract(nit).value());
    }
    nodes.clear();

    // Generate the proper mesh using the associated boundary elements
    this->submesh.resize(num_elems);
    this->boundary_mesh.resize(num_elems);
    size_t cur_elem = 0;
    ElementShape sh;
    if(this->bound_elem_higher_order == nullptr){
        sh.nodes.resize(bound_nodes_per_elem);
        for(size_t i = 0; i < boundary_elements.size(); ++i){
            if(apply_spring[i]){
                this->submesh[cur_elem] = &boundary_elements[i];
                const auto& b = boundary_elements[i];
                for(size_t j = 0; j < bound_nodes_per_elem; ++j){
                    const auto& n = b.nodes[j];
                    gp_Pnt p = utils::change_point(n->point, this->rot3D.transpose());
                    MeshNode* nn = std::find_if(this->boundary_nodes.begin(), this->boundary_nodes.end(), NodeComp(p))->get();
                    sh.nodes[j] = nn;
                }
                this->boundary_mesh[cur_elem].reset(this->boundary_elem_info->make_element(sh, boundary_elements[i].parent));
                ++cur_elem;
            }
        }
    } else {
        std::vector<BoundaryElement> applied;
        const size_t orig_bound_nodes = this->elem_info->get_boundary_nodes_per_element();
        applied.reserve(num_elems);
        std::vector<const MeshNode*> nodes(orig_bound_nodes);
        for(size_t i = 0; i < boundary_elements.size(); ++i){
            if(apply_spring[i]){
                const auto& b = boundary_elements[i];
                for(size_t j = 0; j < orig_bound_nodes; ++j){
                    nodes[j] = b.nodes[j];
                }
                applied.emplace_back(nodes, b.parent, b.normal, b.geom_id);
                this->submesh[cur_elem] = &boundary_elements[i];
                ++cur_elem;
            }
        }
        cur_elem = 0;

        auto new_elems = this->increase_element_order(applied);
        bound_nodes_per_elem = this->boundary_elem_info->get_nodes_per_element();
        sh.nodes.resize(bound_nodes_per_elem);

        for(size_t i = 0; i < new_elems.size(); ++i){
            const auto& b = new_elems[i];
            for(size_t j = 0; j < bound_nodes_per_elem/2; ++j){
                const auto& n = b.nodes[j];
                gp_Pnt p = utils::change_point(n->point, this->rot3D.transpose());
                MeshNode* nn = std::find_if(this->boundary_nodes.begin(), this->boundary_nodes.end(), NodeComp(p))->get();
                sh.nodes[j] = nn;
            }
            for(size_t j = bound_nodes_per_elem/2; j < bound_nodes_per_elem; ++j){
                const auto& n = b.nodes[j];
                gp_Pnt p = n->point;
                MeshNode* nn = std::find_if(this->boundary_nodes.begin(), this->boundary_nodes.end(), NodeComp(p))->get();
                sh.nodes[j] = nn;
            }
            this->boundary_mesh[cur_elem].reset(this->boundary_elem_info->make_element(sh, boundary_elements[i].parent));
            ++cur_elem;
        }
    }

    long npos = 0;
    long id = 0;
    for(auto& l:this->line_nodes){
        l->u_pos[0] = -1;
    }
    for(size_t i = 0; i < this->boundary_nodes.size(); ++i){
        if(this->boundary_nodes[i]->u_pos[0] > -1){
            this->boundary_nodes[i]->u_pos[0] = npos;
            ++npos;
        }
        this->boundary_nodes[i]->id = id;
        ++id;
    }
    this->phi_size = npos;
    nodes.clear();
}

class ElementEdge{
    public:
    std::array<const Node*, 2> vertices;
    std::vector<const BoundaryElement*> parents;

    bool is_connected_to(const ElementEdge& l) const{
        return (this->vertices[0]->point.IsEqual(l.vertices[0]->point, Precision::Confusion()) || this->vertices[1]->point.IsEqual(l.vertices[1]->point, Precision::Confusion())) ||
               (this->vertices[1]->point.IsEqual(l.vertices[0]->point, Precision::Confusion()) || this->vertices[0]->point.IsEqual(l.vertices[1]->point, Precision::Confusion()));
    }

    bool operator==(const ElementEdge& l) const{
        return (this->vertices[0]->point.IsEqual(l.vertices[0]->point, Precision::Confusion()) && this->vertices[1]->point.IsEqual(l.vertices[1]->point, Precision::Confusion())) ||
               (this->vertices[1]->point.IsEqual(l.vertices[0]->point, Precision::Confusion()) && this->vertices[0]->point.IsEqual(l.vertices[1]->point, Precision::Confusion()));
    }
};

std::vector<BoundaryElement> InternalLoads::increase_element_order(const std::vector<BoundaryElement>& boundary_elements){
    std::vector<BoundaryElement> new_elements;
    new_elements.reserve(boundary_elements.size());
    std::list<ElementEdge> edges;
    const size_t N = (this->elem_info->get_shape_type() == Element::Shape::TRI) ? 3 : 4;

    for(auto& e : boundary_elements){
        for(size_t i = 0; i < N; ++i){
            const size_t j = (i+1) % N;
            const auto& n1 = e.nodes[i];
            const auto& n2 = e.nodes[j];
            // TODO: detect inner boundaries
            ElementEdge d{{n1, n2}, {&e}};
            auto it = std::find(edges.begin(), edges.end(), d);
            if(it == edges.end()){
                edges.push_back(d);
            } else {
                it->parents.push_back(&e);
            }
        }
    }

    std::vector<std::vector<const MeshNode*>> new_nodes(boundary_elements.size());
    std::map<const BoundaryElement*, size_t> pointer_map;
    for(size_t i = 0; i < boundary_elements.size(); ++i){
        auto& vec = new_nodes[i];
        vec.reserve(2*N);
        pointer_map[&boundary_elements[i]] = i;
        for(size_t j = 0; j < N; ++j){
            vec.push_back(boundary_elements[i].nodes[j]);
        }
    }

    const size_t prev_size = this->boundary_nodes.size();
    this->boundary_nodes.resize(prev_size + edges.size());
    auto it = edges.begin();
    for(size_t i = 0; i < edges.size(); ++i){
        gp_Pnt p(it->vertices[0]->point);
        p.BaryCenter(1, it->vertices[1]->point, 1);
        p = utils::change_point(p, this->rot3D.transpose());
        std::unique_ptr<MeshNode> node(std::make_unique<MeshNode>(p, 0, 1));
        this->boundary_nodes[prev_size + i] = std::move(node);
        for(auto& p:it->parents){
            new_nodes[pointer_map[p]].push_back(this->boundary_nodes[prev_size + i].get());
        }
        ++it;
    }
    for(size_t i = 0; i < boundary_elements.size(); ++i){
        auto& b = boundary_elements[i];
        new_elements.emplace_back(new_nodes[i], b.parent, b.normal, b.geom_id);
    }
    
    return new_elements;
}
