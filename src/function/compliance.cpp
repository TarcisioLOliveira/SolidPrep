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
#include <mpich-x86_64/mpi.h>

namespace function{

Compliance::Compliance(const Meshing* const mesh, double pc, double psi):
    pc(pc), psi(psi), mesh(mesh){}

double Compliance::calculate(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x){
    (void)x;
    (void)op;
    const double c = cblas_ddot(u.size(), this->mesh->load_vector.data(), 1, u.data(), 1);

    return c;
}
double Compliance::calculate_with_gradient(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
    (void)op;
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return 0;
    }

    auto rho_it = x.cbegin();
    auto grad_it = grad.begin();
    for(const auto& g:this->mesh->geometries){
        const size_t num_den = g->number_of_densities_needed();
        const size_t num_mat = g->number_of_materials();
        if(g->do_topopt){
            if(num_mat == 1){
                const auto D = g->materials.get_D();
                for(const auto& e:g->mesh){
                    const double uKu = e->get_compliance(D, this->mesh->thickness, u);
                    *grad_it = -(1-K_MIN)*pc*std::pow(*rho_it, pc-1)*uKu;

                    ++grad_it;
                    ++rho_it;
                }
            } else {
                std::vector<std::vector<double>> gradD(num_den, g->materials.get_D());
                for(const auto& e:g->mesh){
                    g->materials.get_gradD(rho_it, g->with_void, pc, K_MIN, psi, gradD);
                    for(auto& D:gradD){
                        const double uKu = e->get_compliance(D, this->mesh->thickness, u);
                        *grad_it = -uKu;

                        ++grad_it;
                    }
                }
            }
        }
    }
    const double c = cblas_ddot(u.size(), this->mesh->load_vector.data(), 1, u.data(), 1);

    return c;
}


}
