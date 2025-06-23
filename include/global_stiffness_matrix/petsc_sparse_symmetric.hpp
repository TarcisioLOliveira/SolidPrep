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

#ifndef PETSC_SPARSE_SYMMETRIC_HPP
#define PETSC_SPARSE_SYMMETRIC_HPP

#include <vector>
#include <petsc.h>
#include "logger.hpp"
#include "math/matrix.hpp"
#include "utils/coo.hpp"
#include "meshing.hpp"
#include "global_stiffness_matrix.hpp"

namespace global_stiffness_matrix{

class PETScSparseSymmetric : public GlobalStiffnessMatrix{
    public:
    enum class Backend{
        CPU,
        CUDA
    };
    PETScSparseSymmetric(double EPS_DISPL):GlobalStiffnessMatrix(EPS_DISPL){}
    virtual ~PETScSparseSymmetric();

    virtual void generate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type) override;

    Mat get_K() const{
        return this->K;
    }

    protected:
    Mat K = 0;
    bool first_time = true;
    size_t u_size, l_num;

    virtual void preallocate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type, const size_t mpi_id) = 0;
    virtual void assemble_matrix(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type, const size_t mpi_id) = 0;
    virtual void zero() = 0;
};

class PETScSparseSymmetricCPU : public PETScSparseSymmetric {
    public:
    PETScSparseSymmetricCPU(double EPS_DISPL):PETScSparseSymmetric(EPS_DISPL){}
    virtual ~PETScSparseSymmetricCPU() = default;

    virtual void dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const override;

    inline virtual void reset_hessian() override{
        MatCopy(this->H, this->K, SAME_NONZERO_PATTERN);
    };

    protected:
    Mat H = 0;
    virtual void preallocate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type, const size_t mpi_id) override;
    virtual void assemble_matrix(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type, const size_t mpi_id) override;
    inline virtual void zero() override{
        MatZeroEntries(this->K);
    }
    virtual void final_flush_matrix() override{
        MatAssemblyBegin(this->K, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(this->K, MAT_FINAL_ASSEMBLY);
    }

    inline virtual void insert_block_symmetric(const math::Matrix& k, const std::vector<long>& posi, const std::vector<long>& posj) override{
        // Requires 64-bit indices
        MatSetValues(this->K, posi.size(), posi.data(), posj.size(), posj.data(), k.data(), ADD_VALUES);
        // I think this is a hack, but it works
        MatSetOption(this->K, MAT_ROW_ORIENTED, PETSC_FALSE);
        MatSetValues(this->K, posj.size(), posj.data(), posi.size(), posi.data(), k.data(), ADD_VALUES);
        MatSetOption(this->K, MAT_ROW_ORIENTED, PETSC_TRUE);
    }

    inline virtual void insert_element_matrix(const math::Matrix& k, const std::vector<long>& pos) override{
        // Requires 64-bit indices
        MatSetValues(this->K, pos.size(), pos.data(), pos.size(), pos.data(), k.data(), ADD_VALUES);
    }

    inline virtual void add_to_matrix(size_t i, size_t j, double val) override{
        // Requires 64-bit indices
        MatSetValue(this->K, i, j, val, ADD_VALUES);
    }
};

class PETScSparseSymmetricCUDA : public PETScSparseSymmetric {
    public:
    PETScSparseSymmetricCUDA(double EPS_DISPL):PETScSparseSymmetric(EPS_DISPL){}
    virtual ~PETScSparseSymmetricCUDA() = default;

    inline virtual void dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const override{
        this->K_coo.dot_vector(v, v_out, true);
    }

    inline virtual void reset_hessian() override;

    protected:
    utils::COO<PetscInt> K_coo;

    virtual void preallocate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type, const size_t mpi_id) override;
    virtual void assemble_matrix(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type, const size_t mpi_id) override;
    inline virtual void zero() override{
        this->K_coo.zero();
    }

    inline virtual void insert_block_symmetric(const math::Matrix& k, const std::vector<long>& posi, const std::vector<long>& posj) override{
        // Requires 64-bit indices
        this->K_coo.insert_block(k, posi, posj, false);
        this->K_coo.insert_block(k, posi, posj, true);
    }

    inline virtual void insert_element_matrix(const math::Matrix& k, const std::vector<long>& pos) override{
        // Requires 64-bit indices
        this->K_coo.insert_matrix_general(k, pos);
    }

    virtual void final_flush_matrix() override{
        MatSetValuesCOO(this->K, this->K_coo.vals.data(), INSERT_VALUES);
    }

    inline virtual void add_to_matrix(size_t i, size_t j, double val) override{
        // Requires 64-bit indices
        this->K_coo.add(i, j, val);
    }
};

}

#endif
