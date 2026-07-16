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
#include "math/matrix.hpp"
#include "utils.hpp"
#include "meshing.hpp"
#include "utils/basis_tensor.hpp"
#include <memory>
#include <set>

Spring::Spring(CrossSection cross_section, double thickness, gp_Dir normal, gp_Dir v, gp_Dir w, utils::DelayedPointerView<Material> mat, std::array<double, 3> L, MeshElementFactory* elem, utils::ProblemType type):
    S(std::move(cross_section)), 
    A(S.get_area()),
    thickness(thickness),
    rot2D({normal.X(), v.X(),
           normal.Y(), v.Y()}, 2, 2),
    rot3D({normal.X(), v.X(), w.X(),
           normal.Y(), v.Y(), w.Y(),
           normal.Z(), v.Z(), w.Z()}, 3, 3),
    rot2D_S(utils::basis_tensor_2D_inv_T(rot2D)),
    rot3D_S(utils::basis_tensor_3D_inv_T(rot3D)),
    mat(mat),
    normal(normal),
    elem_info(elem),
    v(v), w(w), L(L), 
    type(type),
    mmat({mat}, type, false)
{

}

math::Matrix Spring::get_K(const BoundaryElement* const b) const{
    const size_t bnode_num = this->elem_info->get_boundary_nodes_per_element();
    std::vector<gp_Pnt> points(bnode_num);

    for(size_t i = 0; i < bnode_num; ++i){
        points[i] = b->nodes[i]->point;
    }
    if(type == utils::PROBLEM_TYPE_2D){
        return b->parent->get_R_spring(this->mmat, rot2D_S, rot2D, L, thickness, points);
    } else {
        return b->parent->get_R_spring(this->mmat, rot3D_S, rot3D, L, thickness, points);
    }

}

math::Matrix Spring::get_dK(const math::Matrix& S, const math::Matrix& dS) const{
    if(type == utils::PROBLEM_TYPE_2D){

        std::array<double, 2> EG = {
            -dS(0,0)/(S(0,0)*S(0,0)), 
            -dS(2,2)/(S(2,2)*S(2,2))
        };

        math::Matrix Korig({EG[0]/L[0], 0,
                            0, EG[1]/L[1]}, 2, 2);

        return rot2D*Korig*rot2D.T();
    } else if(type == utils::PROBLEM_TYPE_3D){

        std::array<double, 3> EG = {
            -dS(0,0)/(S(0,0)*S(0,0)), 
            -dS(3,3)/(S(3,3)*S(3,3)), 
            -dS(4,4)/(S(4,4)*S(4,4))
        };

        math::Matrix Korig({EG[0]/L[0], 0, 0,
                            0, EG[1]/L[1], 0,
                            0, 0, EG[2]/L[2]}, 3, 3);

        return rot3D*Korig*rot3D.T();
    }

    return math::Matrix();
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
}
