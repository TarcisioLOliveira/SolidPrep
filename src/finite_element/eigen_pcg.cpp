/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
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

#include <Eigen/IterativeLinearSolvers>
#include <mpich-x86_64/mpi.h>
#include "finite_element/eigen_pcg.hpp"
#include "math/matrix.hpp"
#include "project_data.hpp"

namespace finite_element{

EigenPCG::EigenPCG(const projspec::DataMap& data):
    FiniteElement(data.proj->contact_data),
    gsm(data.proj->contact_data.EPS_DISPL_SIMPLE)
{
    this->set_global_matrix(&this->gsm);    
}

void EigenPCG::generate_matrix_base(const Meshing* const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const ContactType type){

    this->l_num = l_num;
    this->u_size = u_size;
    this->gsm.generate(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, lambda, type);
    auto& K = this->gsm.get_K();
    this->cg.compute(K);
}

void EigenPCG::solve(std::vector<double>& load){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return;
    }

    Eigen::VectorXd f(load.size());
    std::copy(load.begin(), load.end(), f.begin());
    Eigen::VectorXd u = cg.solve(f);

    std::copy(u.cbegin(), u.cend(), load.begin());
}

void EigenPCG::reset_hessian(){
    this->gsm.reset_hessian();
}

using namespace projspec;
const bool EigenPCG::reg = Factory<FiniteElement>::add(
    [](const DataMap& data){
        return std::make_unique<EigenPCG>(data);
    },
    ObjectRequirements{
        "eigen_pcg",
        {}
    }
);

}
