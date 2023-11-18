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
#ifndef SPRING_HPP
#define SPRING_HPP

#include <algorithm>
#include <array>
#include <Eigen/Core>
#include "element_factory.hpp"
#include "material.hpp"
#include "utils.hpp"
#include "cross_section.hpp"
#include "curvature.hpp"
#include "utils/boundary_nullifier.hpp"

class BoundaryElement;

class Spring{
    public:

    // ASSUMES CROSS SECTION IS PLANE
    // ASSUMES this->normal == this->S.normal

    Spring(CrossSection cross_section, double thickness, gp_Dir normal, gp_Dir v, gp_Dir w, Material* mat, std::array<double, 3> L, std::array<double, 3> F, std::array<double, 3> curv, MeshElementFactory* elem, BoundaryMeshElementFactory* bound_elem, utils::ProblemType type);
    Spring(Spring&&) = default;

    inline void clear_curvature_data(){
        this->boundary_mesh.clear();
        this->boundary_nodes.clear();
        this->curvature.reset();
    }

    std::vector<double> get_K(const gp_Pnt& p) const;

    void calculate_curvature(std::vector<BoundaryElement>& boundary_elements){
        this->generate_mesh(boundary_elements);
        this->curvature->generate_curvature_3D(this->boundary_mesh, this->phi_size, this->boundary_nodes.size());

        this->curv = this->curvature->get_curvatures();
        this->center = this->curvature->get_center();
    }

    void apply_load_2D(std::vector<double>& load_vector) const;
    void apply_load_3D(std::vector<double>& load_vector) const;

    const CrossSection S;
    const double A;
    const double thickness;
    const Eigen::Matrix<double, 2, 2> rot2D;
    const Eigen::Matrix<double, 3, 3> rot3D;
    const std::array<double, 3> F;
    const std::array<double, 3> M;
    const Material* mat;
    const gp_Dir normal;
    const MeshElementFactory* elem_info;
    const BoundaryMeshElementFactory* boundary_elem_info;
    std::vector<const BoundaryElement*> submesh;
    std::unique_ptr<Curvature> curvature;

    private:
    const gp_Dir v;
    const gp_Dir w;
    const std::array<double, 3> L;
    const utils::ProblemType type;

    std::array<double, 2> curv;
    gp_Pnt center;
    size_t phi_size;

    std::vector<std::unique_ptr<MeshNode>> boundary_nodes;
    std::vector<std::unique_ptr<BoundaryMeshElement>> boundary_mesh;

    std::vector<utils::LineBoundary> line_bounds;
    std::vector<const Node*> line_nodes;

    void generate_mesh(std::vector<BoundaryElement>& boundary_elements);

    void generate_boundary();

    // Base it off Meshing::apply_springs
    void apply_load(const std::vector<double>& load_vector) const;
};

#endif
