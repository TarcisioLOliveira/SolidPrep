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
#include <memory>
#include <set>

Spring::Spring(CrossSection cross_section, double thickness, gp_Dir normal, gp_Dir v, gp_Dir w, Material* mat, std::array<double, 3> L, MeshElementFactory* elem, utils::ProblemType type):
    S(std::move(cross_section)), 
    A(S.get_area()),
    thickness(thickness),
    rot2D{{normal.X(), v.X()},
          {normal.Y(), v.Y()}},
    rot3D{{normal.X(), v.X(), w.X()},
         {normal.Y(), v.Y(), w.Y()},
         {normal.Z(), v.Z(), w.Z()}},
    mat(mat),
    normal(normal),
    elem_info(elem),
    v(v), w(w), L(L), 
    type(type){

}

std::vector<double> Spring::get_K(const MeshElement* const e, const gp_Pnt& p) const{
    if(type == utils::PROBLEM_TYPE_2D){

        auto EG = mat->beam_EG_2D(e, p, this->rot2D);

        Eigen::Matrix<double, 2, 2> Korig{{EG[0]/L[0], 0},
                                          {0, EG[1]/L[1]}};

        Eigen::Matrix<double, 2, 2> Ktmp = rot2D*Korig*rot2D.transpose();

        std::vector<double> K(4);
        std::copy(Ktmp.data(), Ktmp.data()+4, K.begin());

        return K;
    } else if(type == utils::PROBLEM_TYPE_3D){

        auto EG = mat->beam_EG_3D(e, p, this->rot3D);

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

    nodes.clear();

    // Generate the proper mesh using the associated boundary elements
    this->submesh.resize(num_elems);
    size_t cur_elem = 0;
    ElementShape sh;
    sh.nodes.resize(bound_nodes_per_elem);
    for(size_t i = 0; i < boundary_elements.size(); ++i){
        if(apply_spring[i]){
            this->submesh[cur_elem] = &boundary_elements[i];
            ++cur_elem;
        }
    }
}
