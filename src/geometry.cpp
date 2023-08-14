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
#include "element_factory.hpp"

Geometry::Geometry(const std::string& path, double scale, utils::ProblemType type,
        MeshElementFactory* elem_type, bool do_topopt, bool with_void, std::vector<Material*> materials):
    shape(utils::load_shape(path, scale)),
    materials(materials, type, with_void),
    element_type(elem_type), do_topopt(do_topopt), with_void(with_void), mesh(), 
    type(type){}

Geometry::Geometry(TopoDS_Shape shape, utils::ProblemType type,
        MeshElementFactory* elem_type, bool do_topopt, bool with_void, std::vector<Material*> materials):
    shape(std::move(shape)),
    materials(materials, type, with_void),
    element_type(elem_type), do_topopt(do_topopt), with_void(with_void), mesh(), 
    type(type){}

void Geometry::get_stresses(const std::vector<double>& u, const double pc, const double psi, std::vector<double>::const_iterator& rho_it, std::vector<double>::iterator& stress_it) const{
    const size_t s_size = this->element_type->get_D_dimension();
    if(this->do_topopt){
        auto D = std::vector<double>(s_size*s_size, 0);
        const size_t num_den = this->number_of_densities_needed();
        if(this->with_void){
            for(const auto& e:this->mesh){
                const gp_Pnt c = e->get_centroid();
                this->materials.get_D(rho_it, psi, c, D);
                const double S = std::pow(*rho_it, pc)*e->get_stress_at(D, c, u);
                *stress_it = S;

                rho_it += num_den;
                ++stress_it;
            }
        } else {
            for(const auto& e:this->mesh){
                const gp_Pnt c = e->get_centroid();
                this->materials.get_D(rho_it, psi, c, D);
                const double S = e->get_stress_at(D, c, u);
                *stress_it = S;

                rho_it += num_den;
                ++stress_it;
            }
        }
    } else {
        for(const auto& e:this->mesh){
            const gp_Pnt c = e->get_centroid();
            const auto D = this->materials.get_D(c);
            double S = e->get_stress_at(D, c, u);
            *stress_it = S;
            ++stress_it;
        }
    }
}
