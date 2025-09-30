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

#ifndef INTERNAL_LOADS_HPP
#define INTERNAL_LOADS_HPP

#include <array>
#include <memory>
#include "element_factory.hpp"
#include "logger.hpp"
#include "material.hpp"
#include "utils.hpp"
#include "cross_section.hpp"
#include "curvature.hpp"
#include "utils/delayed_pointer.hpp"

class BoundaryElement;

class InternalLoads{
    public:

    // ASSUMES CROSS SECTION IS PLANE
    // ASSUMES this->normal == this->S.normal

    InternalLoads(CrossSection cross_section, double thickness, gp_Dir normal, gp_Dir v, gp_Dir w, utils::DelayedPointerView<Material> mat, std::array<double, 3> F, std::array<double, 3> M, MeshElementFactory* elem, utils::ProblemType type);
    InternalLoads(InternalLoads&&) = default;

    inline void clear_curvature_data(){
        this->boundary_mesh.clear();
        this->boundary_nodes.clear();
        this->curvature.reset();
    }
    inline void set_calculate_adjoint(){
        this->calculate_adjoint = true;
    }
    inline bool get_calculate_adjoint() const{
        return this->calculate_adjoint;
    }

    void calculate_curvature(std::vector<BoundaryElement>& boundary_elements);
    void curvature_adjoint(const std::vector<double>& u){
        this->curvature->solve_adjoint_problem(this->boundary_mesh, u);
    }

    inline double get_adjoint_gradient(const std::vector<double>& fea_u, const math::Matrix& S, const math::Matrix& D, const math::Matrix& dS, const math::Matrix& dD, const BoundaryMeshElement* const e) const{
        return this->curvature->get_adjoint_gradient(fea_u, S, D, dS, dD, e);
    }
    void apply_load_2D(const std::vector<long>& node_positions, std::vector<double>& load_vector) const;
    void apply_load_3D(const std::vector<long>& node_positions, std::vector<double>& load_vector) const;

    const CrossSection S;
    const double A;
    const double thickness;
    const math::Matrix rot2D;
    const math::Matrix Lek_basis;
    const math::Matrix rot3D;
    const std::array<double, 3> F;
    const std::array<double, 3> M;
    utils::DelayedPointerView<Material> mat;
    const gp_Dir normal;
    const MeshElementFactory* elem_info;
    std::vector<const BoundaryElement*> submesh;
    std::unique_ptr<Curvature> curvature;
    std::vector<std::unique_ptr<MeshNode>> boundary_nodes;
    std::vector<std::unique_ptr<BoundaryMeshElement>> boundary_mesh;

    private:
    bool calculate_adjoint = false;
    const gp_Dir v;
    const gp_Dir w;
    const utils::ProblemType type;

    std::array<double, 2> curv;
    gp_Pnt center;
    size_t phi_size;


    std::vector<const Node*> line_nodes;

    void generate_mesh(const std::vector<BoundaryElement>& boundary_elements);

    // Base it off Meshing::apply_springs
    void apply_load(const std::vector<double>& load_vector) const;
};

#endif
