/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   Foobar is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Foobar is distributed in the hope that it will be useful,
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

void FiniteElement::calculate_stresses(Meshing* mesh, const std::vector<float>& displacements, const std::vector<double>& density) const{
    logger::quick_log("Calculating stresses...");
    for(auto& e:mesh->element_list){
        for(size_t n = 0; n < e->nodes.size(); ++n){
            MeshNode* node = e->get_node(n);
            for(size_t k = 0; k < node->get_result_size(); ++k){
                node->results[k] = 0;
            }
        }
    }
    for(size_t i = 0; i < mesh->element_list.size(); ++i){
        auto& e = mesh->element_list[i];
        for(size_t n = 0; n < e->nodes.size(); ++n){
            bool calculated = false;
            MeshNode* node = e->get_node(n);
            for(size_t k = 0; k < node->get_result_size(); ++k){
                if(node->results[k] != 0){
                    calculated = true;
                    break;
                }
            }
            double rho = 1;
            if(density.size() > 0){
                rho = density[i];
            }
            if(!calculated){
                e->get_stresses(n, displacements, rho);
            }
        }
    }
    logger::quick_log("Done");
}
void FiniteElement::calculate_forces(Meshing* mesh, const std::vector<float>& displacements) const{
    logger::quick_log("Calculating forces...");
    for(auto& e:mesh->element_list){
        for(size_t n = 0; n < e->nodes.size(); ++n){
            MeshNode* node = e->get_node(n);
            for(size_t k = 0; k < node->get_result_size(); ++k){
                node->results[k] = 0;
            }
        }
    }
    for(auto& e:mesh->element_list){
        for(size_t n = 0; n < e->nodes.size(); ++n){
            bool calculated = false;
            MeshNode* node = e->get_node(n);
            for(size_t k = 0; k < node->get_result_size(); ++k){
                if(node->results[k] != 0){
                    calculated = true;
                    break;
                }
            }
            if(!calculated){
                e->get_internal_loads(n, displacements);
            }
        }
    }
    logger::quick_log("Done.");
}
