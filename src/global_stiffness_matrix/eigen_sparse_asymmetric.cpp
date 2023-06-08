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

    virtual void generate(const Meshing * const mesh, const std::vector<double>& density, const double pc, const double psi) override;

    std::vector<T> triplets;

    protected:
    bool first_time = true;
    utils::SparseMatrix K;

    virtual inline void insert_element_matrix(const std::vector<double>& k, const std::vector<long>& pos) override{
        this->K.insert_matrix_general_mumps(k, pos);
    }
};

void EigenSparseAsymmetricTriplets::generate(const Meshing * const mesh, const std::vector<double>& density, const double pc, const double psi){
    this->generate_base(mesh, density, pc, psi);
    this->triplets = K.get_eigen_triplets();
    this->K.clear();
}

}

void EigenSparseAsymmetric::generate(const Meshing* const mesh, const std::vector<double>& density, const double pc, const double psi){
    logger::quick_log("Generating stiffness matrix...");
    if(this->first_time){
        this->K = Mat(mesh->load_vector.size(), mesh->load_vector.size());
        internal::EigenSparseAsymmetricTriplets trigen;
        trigen.generate(mesh, density, pc, psi);
        this->K.setFromTriplets(trigen.triplets.begin(), trigen.triplets.end());
        this->first_time = false;
    } else {
        std::fill(this->K.valuePtr(), this->K.valuePtr() + this->K.nonZeros(), 0);
        this->generate_base(mesh, density, pc, psi);
    }
    logger::quick_log("Done.");
}

}
