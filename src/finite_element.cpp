/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
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

#include "finite_element.hpp"
#include "logger.hpp"
#include "project_data.hpp"

std::vector<double> FiniteElement::calculate_forces(const Meshing* mesh, const std::unique_ptr<Geometry>& geometry, const double t, const std::vector<double>& displacements, const std::unique_ptr<MeshElementFactory>& elem_maker) const{
    logger::quick_log("Calculating forces...");
    size_t dof = elem_maker->get_dof_per_node();
    std::vector<double> results(mesh->node_list.size()*dof, 0);
   
    const auto D = geometry->get_D(); 
    for(auto& e:mesh->element_list){
        auto f = e->get_internal_loads(D, t, displacements);
        for(size_t n = 0; n < e->nodes.size(); ++n){
            auto& node = e->nodes[n];
            for(size_t i = 0; i < dof; ++i){
                results[node->id*dof + i] += f[n*dof+i];
            }
        }
    }
    logger::quick_log("Done.");

    return results;
}
