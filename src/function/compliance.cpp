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

#include "function/compliance.hpp"
#include <cblas.h>

namespace function{

Compliance::Compliance(const Meshing* const mesh, double pc):
    pc(pc), mesh(mesh){}

double Compliance::calculate(const std::vector<double>& u, const std::vector<double>& x){
    (void)x;
    const double c = cblas_ddot(u.size(), this->mesh->load_vector.data(), 1, u.data(), 1);

    return c;
}
double Compliance::calculate_with_gradient(const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
    size_t i = 0;
    for(const auto& g:this->mesh->geometries){
        const size_t num_den = g->number_of_densities_needed();
        const size_t num_mat = g->number_of_materials();
        if(g->do_topopt){
            if(num_den == 1){
                const auto D = g->get_D(0);
                size_t j = i;
                for(const auto& e:g->mesh){
                    double uKu = e->get_compliance(D, this->mesh->thickness, u);
                    grad[i] = -pc*std::pow(x[i], pc-1)*uKu;

                    ++i;
                }
                if(num_mat == 2){
                    const auto D = g->get_D(1);
                    for(const auto& e:g->mesh){
                        double uKu = e->get_compliance(D, this->mesh->thickness, u);
                        grad[j] += pc*std::pow(1 - x[j], pc-1)*uKu;

                        ++j;
                    }
                }
            } else {
                logger::log_assert(num_den == 1, logger::ERROR, "minimal compliance problems that require more than 1 design variables are currently not supported.");
            }
        }
    }
    const double c = cblas_ddot(u.size(), this->mesh->load_vector.data(), 1, u.data(), 1);

    return c;
}


}
