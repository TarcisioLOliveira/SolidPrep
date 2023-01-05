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

#include "geometry.hpp"
#include "logger.hpp"

Geometry::Geometry(const std::string& path, double scale, utils::ProblemType type,
        MeshElementFactory* elem_type, bool do_topopt, bool with_void, Material* material, std::vector<Material*> alt_materials ):
    shape(utils::load_shape(path, scale)),
    material(material), alternative_materials(alt_materials.begin(), alt_materials.end()),
    element_type(elem_type), do_topopt(do_topopt), with_void(with_void),
    constitutive_matrices(this->init_constitutive_matrices(type)), mesh(), 
    type(type){}

Geometry::Geometry(TopoDS_Shape shape, utils::ProblemType type,
        MeshElementFactory* elem_type, bool do_topopt, bool with_void, Material* material, std::vector<Material*> alt_materials ):
    shape(std::move(shape)),
    material(material), alternative_materials(alt_materials.begin(), alt_materials.end()),
    element_type(elem_type), do_topopt(do_topopt), with_void(with_void),
    constitutive_matrices(this->init_constitutive_matrices(type)), mesh(), 
    type(type){}

void Geometry::get_stresses(const std::vector<double>& u, std::vector<double>::iterator& stress_it) const{
    const size_t num_den = this->number_of_densities_needed();
    const size_t num_mat = this->number_of_materials();
    if(this->do_topopt){
        if(num_den == 1){
            if(num_mat == 1){
                const auto D = this->get_D(0);
                for(const auto& e:this->mesh){
                    const double S = e->get_stress_at(D, e->get_centroid(), u);
                    *stress_it = S;

                    ++stress_it;
                }
            } else {
                std::vector<std::vector<double>> D(num_mat);
                for(size_t i = 0; i < num_mat; ++i){
                    D[i] = this->get_D(i);
                }
                for(const auto& e:this->mesh){
                    for(size_t i = 0; i < num_mat; ++i){
                        const double S = e->get_stress_at(D[i], e->get_centroid(), u);
                        *stress_it = S;

                        ++stress_it;
                    }
                }
            }
        } else {
            logger::log_assert(num_den == 1, logger::ERROR, "problems that require more than 1 design variables are currently not supported.");
        }
    } else {
        const auto D = this->get_D(0);
        for(const auto& e:this->mesh){
            double S = e->get_stress_at(D, e->get_centroid(), u);
            *stress_it = S;
            ++stress_it;
        }
    }
}
