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

#include "optimizer.hpp"

void Optimizer::initialize_optimizer(const Meshing* const mesh){
    this->number_of_elements = this->get_number_of_elements(mesh->geometries);
    this->volumes.resize(this->number_of_elements);
    this->stresses.resize(this->number_of_elements);

    this->get_volumes(mesh->geometries, mesh->thickness, this->volumes);
}

TopoDS_Shape Optimizer::make_shape(const std::vector<double>& x, const std::vector<Geometry*>& geometries, const double result_threshold) const{
    logger::quick_log(" "); 
    logger::quick_log("Saving resulting geometries...");
    // TODO: make it work with multiple geometries
    std::cout << "\r" << 0 << "%         ";
    TopoDS_Shape result = BRepBuilderAPI_Copy(geometries[0]->shape);
    for(size_t i = 0; i < x.size(); ++i){
        if(x[i] >= result_threshold){
            result = utils::cut_shape(result, geometries[0]->mesh[i]->get_shape());
        }
        double pc = i/(double)(x.size()-1);
        std::cout << "\r" << pc*100 << "%         ";
    }
    result = utils::cut_shape(geometries[0]->shape, result);
    std::cout << "\r" << 100 << "%         ";
    logger::quick_log(" "); 
    return result;
}

void Optimizer::get_stresses(const std::vector<Geometry*> geometries, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& stresses) const{
    (void)x;
    auto stress_it = stresses.begin();
    for(const auto& g:geometries){
        g->get_stresses(u, stress_it);
    }
}

void Optimizer::get_volumes(const std::vector<Geometry*> geometries, const double thickness, std::vector<double>& volumes) const{
    auto V_it = volumes.begin();
    for(const auto& g:geometries){
        for(const auto& e:g->mesh){
            *V_it = e->get_volume(thickness);
            ++V_it;
        }
    }
}

size_t Optimizer::get_number_of_elements(const std::vector<Geometry*> geometries) const{
    size_t elem_num = 0;
    for(const auto& g:geometries){
        elem_num += g->mesh.size();
    }

    return elem_num;
}

void Optimizer::apply_densities(const std::vector<Geometry*> geometries, const std::vector<double>& x, std::vector<double>& vals, const double pc) const{
    if(pc == 0){
        return;
    } else if(pc == 1){
        auto x_it = x.cbegin();
        auto v_it = vals.begin();
        for(const auto& g:geometries){
            if(g->do_topopt){
                for(auto xi = x_it; xi < x_it + g->mesh.size(); ++xi, ++v_it){
                    *v_it *= *xi;
                }
                x_it += g->mesh.size();
            } else {
                v_it += g->mesh.size();
            }
        }
    } else {
        auto x_it = x.cbegin();
        auto v_it = vals.begin();
        for(const auto& g:geometries){
            if(g->do_topopt){
                for(auto xi = x_it; xi < x_it + g->mesh.size(); ++xi, ++v_it){
                    *v_it *= std::pow(*xi, pc);
                }
                x_it += g->mesh.size();
            } else {
                v_it += g->mesh.size();
            }
        }
    }
}
