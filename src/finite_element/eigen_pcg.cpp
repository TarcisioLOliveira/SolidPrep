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
#include "logger.hpp"

namespace finite_element{

EigenPCG::EigenPCG(NonlinearSolver* nl):FiniteElement(nl), gsm(){}

void EigenPCG::generate_matrix_base(const Meshing* const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache, const MatrixType type){

    this->l_num = l_num;
    this->u_size = u_size;
    this->gsm.generate(mesh, u_size, l_num, node_positions, topopt, D_cache, type);
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
bool EigenPCG::generate_hessian(std::vector<double>& lambda, const std::vector<double>& Ku){
    bool mod = this->gsm.generate_hessian(lambda, Ku);
    auto& K = this->gsm.get_K();
    this->cg.compute(K);
    return mod;
}
void EigenPCG::dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const{
    this->gsm.dot_vector(v, v_out);
}
double EigenPCG::get_newton_step(const std::vector<double>& delta, const std::vector<double>& lambda, const std::vector<double>& Ku){
    return this->gsm.get_newton_step(delta, lambda, Ku);
}

}
