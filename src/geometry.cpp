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
        MeshElementFactory* elem_type, bool do_topopt, bool with_void, std::vector<Material*> materials):
    shape(utils::load_shape(path, scale)),
    materials(materials, type),
    element_type(elem_type), do_topopt(do_topopt), with_void(with_void), mesh(), 
    type(type){}

Geometry::Geometry(TopoDS_Shape shape, utils::ProblemType type,
        MeshElementFactory* elem_type, bool do_topopt, bool with_void, std::vector<Material*> materials):
    shape(std::move(shape)),
    materials(materials, type),
    element_type(elem_type), do_topopt(do_topopt), with_void(with_void), mesh(), 
    type(type){}

void Geometry::get_stresses(const std::vector<double>& u, const double pc, const double p_min, const double psi, std::vector<double>::const_iterator& rho_it, std::vector<double>::iterator& stress_it) const{
    const size_t num_mat = this->number_of_materials();
    if(this->do_topopt){
        if(num_mat == 1){
            const auto D = this->materials.get_D();
            for(const auto& e:this->mesh){
                const double S = e->get_stress_at(D, e->get_centroid(), u);
                *stress_it = S;

                ++stress_it;
            }
        } else {
            std::vector<double> D = this->materials.get_D();
            for(const auto& e:this->mesh){
                this->materials.get_D(rho_it, this->with_void, pc, p_min, psi, D);
                const double S = e->get_stress_at(D, e->get_centroid(), u);
                *stress_it = S;

                ++stress_it;
            }
        }
    } else {
        const auto D = this->materials.get_D();
        for(const auto& e:this->mesh){
            double S = e->get_stress_at(D, e->get_centroid(), u);
            *stress_it = S;
            ++stress_it;
        }
    }
}
