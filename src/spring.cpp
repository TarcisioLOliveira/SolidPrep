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

#include "logger.hpp"
#include "spring.hpp"
#include "utils.hpp"
#include "meshing.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <memory>
#include <set>

Spring::Spring(CrossSection cross_section, double thickness, gp_Dir normal, gp_Dir v, gp_Dir w, Material* mat, std::array<double, 3> L, std::array<double, 3> F, std::array<double, 3> M, MeshElementFactory* elem, BoundaryMeshElementFactory* bound_elem, utils::ProblemType type):
    S(std::move(cross_section)), 
    A(S.get_area()),
    thickness(thickness),
    rot2D{{normal.X(), v.X()},
          {normal.Y(), v.Y()}},
    rot3D{{normal.X(), v.X(), w.X()},
         {normal.Y(), v.Y(), w.Y()},
         {normal.Z(), v.Z(), w.Z()}},
    F(F), M(M),
    mat(mat),
    normal(normal),
    elem_info(elem),
    boundary_elem_info(bound_elem),
    v(v), w(w), L(L), 
    type(type){

    this->curvature = std::make_unique<Curvature>(mat, normal, v, w, rot2D, rot3D, this->boundary_elem_info, F[1], F[2], M[0], M[1], M[2]);
}

void Spring::apply_load_2D(std::vector<double>& load_vector) const{
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
        const auto E = this->mat->beam_E_2D(e->parent, c, this->normal);
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
void Spring::apply_load_3D(std::vector<double>& load_vector) const{
    const size_t bound_nodes_per_elem = this->elem_info->get_boundary_nodes_per_element();
    const size_t nodes_per_elem = this->elem_info->get_nodes_per_element();
    const size_t dof = this->elem_info->get_dof_per_node();

    const size_t N = 3;
    Eigen::Matrix<double, N, N> S{{0, 0, 0},
                                  {0, 0, 0},
                                  {0, 0, 0}};
    std::vector<double> Sn(N*N);

    double t_uv = 0, t_uw = 0;

    double FF = 0;
    for(size_t j = 0; j < submesh.size(); ++j){
        const auto& e = submesh[j];
        const auto& b = this->boundary_mesh[j];
        std::vector<gp_Pnt> points(bound_nodes_per_elem);
        for(size_t i = 0; i < bound_nodes_per_elem; ++i){
            points[i] = e->nodes[i]->point;
        }

        const gp_Pnt c = e->get_centroid(bound_nodes_per_elem);
        const auto E = this->mat->beam_E_3D(e->parent, c, this->normal);
        S(0,1) = E*this->curv[0];
        S(0,2) = E*this->curv[1];
        this->curvature->get_shear_in_3D(b.get(), t_uv, t_uw);
        Eigen::Vector<double, N> F0{this->F[0]/A, t_uv, t_uw};
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
        for(auto& ff:Rf){
            FF += ff;
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

std::vector<double> Spring::get_K(const MeshElement* const e, const gp_Pnt& p) const{
    if(type == utils::PROBLEM_TYPE_2D){

        auto EG = mat->beam_EG_2D(e, p, this->normal);

        Eigen::Matrix<double, 2, 2> Korig{{EG[0]/L[0], 0},
                                          {0, EG[1]/L[1]}};

        Eigen::Matrix<double, 2, 2> Ktmp = rot2D*Korig*rot2D.transpose();

        std::vector<double> K(4);
        std::copy(Ktmp.data(), Ktmp.data()+4, K.begin());

        return K;
    } else if(type == utils::PROBLEM_TYPE_3D){

        auto EG = mat->beam_EG_3D(e, p, this->normal);

        Eigen::Matrix<double, 3, 3> Korig{{EG[0]/L[0], 0, 0},
                                          {0, EG[1]/L[1], 0},
                                          {0, 0, EG[2]/L[2]}};

        Eigen::Matrix<double, 3, 3> Ktmp = rot3D*Korig*rot3D.transpose();

        std::vector<double> K(9);
        std::copy(Ktmp.data(), Ktmp.data()+9, K.begin());

        return K;
    }

    return std::vector<double>();
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

void Spring::generate_mesh(std::vector<BoundaryElement>& boundary_elements){
    const size_t bound_nodes_per_elem = this->elem_info->get_boundary_nodes_per_element();

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

    this->generate_boundary();

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

    this->generate_boundary();

    long npos = 0;
    long id = 0;
    for(size_t i = 0; i < this->boundary_nodes.size(); ++i){
        const auto it = std::find(this->line_nodes.begin(), this->line_nodes.end(), this->boundary_nodes[i].get());
        if(it != this->line_nodes.end()){
            this->boundary_nodes[i]->u_pos[0] = -1;
        } else {
            this->boundary_nodes[i]->u_pos[0] = npos;
            ++npos;
        }
        this->boundary_nodes[i]->id = id;
        ++id;
    }
    nodes.clear();
    this->phi_size = npos;
}

void Spring::generate_boundary() {
    std::list<utils::LineBoundary> bound_tmp;
    const size_t N = (this->elem_info->get_shape_type() == Element::Shape::TRI) ? 3 : 4;

    for(const auto& e : this->boundary_mesh){
        for(size_t i = 0; i < N; ++i){
            const size_t j = (i+1) % N;
            const auto& n1 = e->nodes[i];
            const auto& n2 = e->nodes[j];
            // TODO: detect inner boundaries
            utils::LineBoundary b{{n1, n2}, false, e.get()};
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
    }
    this->line_nodes.resize(nodes_tmp.size());
    std::copy(nodes_tmp.begin(), nodes_tmp.end(), this->line_nodes.begin());
}
