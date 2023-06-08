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

#ifndef EIGEN_PCG_HPP
#define EIGEN_PCG_HPP

#include "finite_element.hpp"
#include "global_stiffness_matrix/eigen_sparse_symmetric.hpp"
#include "global_stiffness_matrix/eigen_sparse_asymmetric.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/IterativeLinearSolvers/BasicPreconditioners.h>
#include <Eigen/src/IterativeLinearSolvers/BiCGSTAB.h>
#include <Eigen/src/IterativeLinearSolvers/ConjugateGradient.h>
#include <Eigen/src/IterativeLinearSolvers/IncompleteLUT.h>

namespace finite_element{

class EigenPCG : public FiniteElement{
    public:
    EigenPCG();

    virtual std::vector<double> calculate_displacements(const Meshing* const mesh, std::vector<double> load, const std::vector<double>& density = std::vector<double>(), double pc = 3, double psi = 0.1) override;

    inline virtual void set_steps(size_t s) override{
        this->steps = s;
        this->u.resize(steps);
    }

    private:
    global_stiffness_matrix::EigenSparseAsymmetric gsm;
    Eigen::ConjugateGradient<global_stiffness_matrix::EigenSparseAsymmetric::Mat, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg;
    //Eigen::BiCGSTAB<global_stiffness_matrix::EigenSparseAsymmetric::Mat, Eigen::DiagonalPreconditioner<double>> cg;
    std::vector<Eigen::VectorXd> u;
};

}

#endif
