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

#include "global_stiffness_matrix/eigen_sparse_asymmetric.hpp"
#include "logger.hpp"
#include "utils/sparse_matrix.hpp"

namespace global_stiffness_matrix{

namespace internal{

class EigenSparseAsymmetricTriplets : public GlobalStiffnessMatrix{
    public:
    typedef Eigen::Triplet<double, std::ptrdiff_t> T;

    virtual ~EigenSparseAsymmetricTriplets() = default;

    virtual void generate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache, const std::vector<double>& u_ext, const FiniteElement::ContactType type) override;

    inline virtual void dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const override{
        (void)v;
        (void)v_out;
        // Unused
    }
    inline virtual void reset_hessian() override{
       // Unused 
    };
    inline virtual bool generate_hessian(std::vector<double>& lambda, const std::vector<double>& Ku) override{
        (void)lambda;
        (void)Ku;
        // Unused
    };
    inline virtual double get_newton_step(const std::vector<double>& delta, const std::vector<double>& lambda, const std::vector<double>& Ku) override{
        (void) delta;
        (void) lambda;
        (void) Ku;
        // Unused
    }

    std::vector<T> triplets;

    protected:
    bool first_time = true;
    utils::SparseMatrix K;

    inline virtual void insert_block_symmetric(const std::vector<double>& k, const std::vector<long>& posi, const std::vector<long>& posj) override{
        this->K.insert_block(k, posi, posj, false);
        this->K.insert_block(k, posi, posj, true);
    }
    virtual inline void insert_element_matrix(const std::vector<double>& k, const std::vector<long>& pos) override{
        this->K.insert_matrix_general_mumps(k, pos);
    }
    inline virtual void add_to_matrix(size_t i, size_t j, double val) override{
        this->K.add(i, j, val);
    }
};

void EigenSparseAsymmetricTriplets::generate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache, const std::vector<double>& u_ext, const FiniteElement::ContactType type){
    (void) u_size;
    (void) l_num;
    this->generate_base(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, type);
    this->triplets = K.get_eigen_triplets();
    this->K.clear();
}

}

void EigenSparseAsymmetric::generate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache, const std::vector<double>& u_ext, const FiniteElement::ContactType type){
    logger::quick_log("Generating stiffness matrix...");
    this->u_size = u_size;
    this->l_num = l_num;
    if(this->first_time){
        size_t matrix_width = u_size;
        if(type != FiniteElement::ContactType::RIGID){
            matrix_width += l_num;
        }
        this->K = Mat(matrix_width, matrix_width);
        internal::EigenSparseAsymmetricTriplets trigen;
        trigen.generate(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, type);
        this->K.setFromTriplets(trigen.triplets.begin(), trigen.triplets.end());
        this->first_time = false;
    } else {
        std::fill(this->K.valuePtr(), this->K.valuePtr() + this->K.nonZeros(), 0);
        this->generate_base(mesh, u_size, l_num, node_positions, topopt, D_cache, u_ext, type);
    }
    if(type == FiniteElement::ContactType::FRICTIONLESS_DISPL){
        this->K_bkp = K;
    }
    logger::quick_log("Done.");
}

bool EigenSparseAsymmetric::generate_hessian(std::vector<double>& lambda, const std::vector<double>& Ku){
    const size_t hoffset = u_size + 2*l_num;
    for(size_t i = 0; i < l_num; ++i){
        const double l = lambda[2*l_num + i];
    }
    for(size_t i = 0; i < l_num; ++i){
        this->K.coeffRef(i + hoffset, i + hoffset) += 2*Ku[i + hoffset];
    }

    // TODO
    return false;
}

double EigenSparseAsymmetric::get_newton_step(const std::vector<double>& delta, const std::vector<double>& lambda, const std::vector<double>& Ku){
    const size_t hoffset = this->u_size + 2*this->l_num;
    double M = 1.0;
    for(size_t i = 0; i < l_num; ++i){
        const size_t ui = hoffset + i;
        const size_t li = 2*l_num + i;
        if(Ku[ui] > 0 || std::abs(delta[ui]) < 1e-14){
            continue;
        }
        double M_test = (std::sqrt((this->K_MIN - 2*Ku[ui])/this->K.coeff(ui,ui)) - lambda[li])/delta[ui];
        if(M_test > 0 && M_test < M){
            M = M_test;
        }
    }
    return M;
}

}
