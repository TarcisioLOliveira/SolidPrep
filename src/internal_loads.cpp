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
#include "boundary_element/BTRI6.hpp"
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

    if(bound_elem->get_element_order() >= 2){
        this->boundary_elem_info = bound_elem;
    } else if(elem_info->get_shape_type() == Element::Shape::TRI){
        this->bound_elem_higher_order.reset(static_cast<BoundaryMeshElementFactory*>(
                    new BoundaryMeshElementFactoryImpl<boundary_element::BTRI6>()
                ));
        this->boundary_elem_info = this->bound_elem_higher_order.get();
    } else {
        // TODO
        logger::log_assert(false, logger::ERROR, "BQ8 not implemented");
    }

    this->curvature = std::make_unique<Curvature>(mat, normal, v, w, rot2D, rot3D, this->boundary_elem_info, F[1], F[2], M[0], M[1], M[2]);
}

void InternalLoads::calculate_curvature(std::vector<BoundaryElement>& boundary_elements){
    this->generate_mesh(boundary_elements);
    this->curvature->generate_curvature_3D(this->boundary_nodes, this->boundary_mesh, this->line_bounds, this->phi_size, this->boundary_nodes.size());

    this->curv = this->curvature->get_curvatures();
    this->center = this->curvature->get_center();
}

void InternalLoads::apply_load_2D(std::vector<double>& load_vector) const{
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
             Rf = submesh[j]->parent->get_Rf(Sn, F, center, this->thickness, list);
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
                    auto n = e->parent->nodes[i];
                    if(n->u_pos[j] >= 0){
                        load_vector[n->u_pos[j]] += Rf[i*dof+j];
                    }
                }
            }
        }
    }
}
void InternalLoads::apply_load_3D(std::vector<double>& load_vector) const{
    const size_t bound_nodes_per_elem = this->elem_info->get_boundary_nodes_per_element();
    const size_t nodes_per_elem = this->elem_info->get_nodes_per_element();
    const size_t dof = this->elem_info->get_dof_per_node();

    const size_t N = 3;
    Eigen::Matrix<double, N, N> S{{0, 0, 0},
                                  {0, 0, 0},
                                  {0, 0, 0}};
    std::vector<double> Sn(N*N);

    double t_yz = 0, t_xz = 0;

    std::vector<double> FF(3, 0);
    for(size_t j = 0; j < submesh.size(); ++j){
        const auto& e = submesh[j];
        const auto& b = this->boundary_mesh[j];
        std::vector<gp_Pnt> points(bound_nodes_per_elem);
        for(size_t i = 0; i < bound_nodes_per_elem; ++i){
            points[i] = e->nodes[i]->point;
        }

        const gp_Pnt c = e->get_centroid(bound_nodes_per_elem);
        const auto E = this->mat->beam_E_3D(e->parent, c, this->rot3D);
        S(2,0) = E*this->curv[0];
        S(2,1) = E*this->curv[1];
        this->curvature->get_shear_in_3D(b.get(), t_xz, t_yz);
        Eigen::Vector<double, N> F0{t_xz, t_yz, this->F[0]/A};
        Eigen::Vector<double, N> Fr = this->rot3D*F0;
        std::vector<double> F{Fr[0], Fr[1], Fr[2]};
        
        Eigen::Matrix<double, N, N> Srot = this->rot3D*S*this->rot3D.transpose();
        for(size_t row = 0; row < N; ++row){
            for(size_t col = 0; col < N; ++col){
                Sn[row*N + col] = Srot(row,col);
            }
        }

        const auto Rf = e->parent->get_Rf(Sn, F, this->center, this->thickness, points);
        //logger::quick_log(fe);
        for(size_t i = 0; i < nodes_per_elem; ++i){
            for(size_t j = 0; j < dof; ++j){
                FF[j] += Rf[i*dof + j];
            }
        }
        for(size_t i = 0; i < nodes_per_elem; ++i){
            for(size_t j = 0; j < dof; ++j){
                const auto n = e->parent->nodes[i];
                if(n->u_pos[j] >= 0){
                    load_vector[n->u_pos[j]] += Rf[i*dof+j];
                }
            }
        }
    }
    logger::quick_log(FF);
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
                applied.emplace_back(nodes, b.parent, b.normal);
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

    this->generate_boundary();

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
    nodes.clear();
    this->phi_size = npos;
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
        new_elements.emplace_back(new_nodes[i], b.parent, b.normal);
    }
    
    return new_elements;
}

void InternalLoads::generate_boundary() {
    const auto colinear =
        [](const gp_Pnt& p1, const gp_Pnt& p2, const gp_Pnt& p3) -> bool{
            gp_Mat(p1.XYZ(), p2.XYZ(), p3.XYZ());

            return std::abs(gp_Mat().Determinant()) < Precision::Confusion();
        };

    std::list<utils::LineBoundary> bound_tmp;
    const size_t N = (this->elem_info->get_shape_type() == Element::Shape::TRI) ? 3 : 4;
    const size_t total_nodes_num = this->boundary_elem_info->get_nodes_per_element();

    for(const auto& e : this->boundary_mesh){
        for(size_t i = 0; i < N; ++i){
            const size_t j = (i+1) % N;
            const auto& n1 = e->nodes[i];
            const auto& n2 = e->nodes[j];
            // TODO: detect inner boundaries
            gp_Vec v(n1->point, n2->point);
            gp_Dir n(v.Y(), -v.X(), 0);
            std::vector<const Node*> other_nodes;
            for(size_t k = N; k < total_nodes_num; ++k){
                const auto& n3 = e->nodes[k];
                if(colinear(n1->point, n2->point, n3->point)){
                    other_nodes.push_back(n3);
                }
            }
            utils::LineBoundary b{{n1, n2}, false, e.get(), n, other_nodes};
            auto it = std::find(bound_tmp.begin(), bound_tmp.end(), b);
            if(it == bound_tmp.end()){
                bound_tmp.push_back(b);
            } else {
                bound_tmp.erase(it);
            }
        }
    }

    this->line_bounds.resize(bound_tmp.size());
    std::move(bound_tmp.begin(), bound_tmp.end(), this->line_bounds.begin());
    

    std::set<const Node*> nodes_tmp;
    for(const auto& l:this->line_bounds){
        nodes_tmp.insert(l.edges.begin(), l.edges.end());
        nodes_tmp.insert(l.other_nodes.begin(), l.other_nodes.end());
    }
    this->line_nodes.resize(nodes_tmp.size());
    std::copy(nodes_tmp.begin(), nodes_tmp.end(), this->line_nodes.begin());
}
