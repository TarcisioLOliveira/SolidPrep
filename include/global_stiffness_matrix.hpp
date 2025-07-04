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

#ifndef GLOBAL_STIFFNESS_MATRIX_HPP
#define GLOBAL_STIFFNESS_MATRIX_HPP

#include "math/matrix.hpp"
#include "meshing.hpp"
#include "finite_element.hpp"

class GlobalStiffnessMatrix{
    public:
    virtual ~GlobalStiffnessMatrix() = default;
    const double K_MIN = 1e-14;
    GlobalStiffnessMatrix(double EPS_DISPL_SIMPLE);

    virtual void generate(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type) = 0;

    virtual void dot_vector(const std::vector<double>& v, std::vector<double>& v_out) const = 0;

    virtual void reset_hessian() = 0;

    virtual void add_frictionless_part2(const Meshing * const mesh, const std::vector<long>& node_positions, const std::vector<double>& u_ext, const std::vector<double>& lambda, const std::vector<math::Matrix>& D_cache, bool topopt, bool stub = false);
    virtual void add_frictionless_simple(const Meshing * const mesh, const std::vector<long>& node_positions, const std::vector<double>& u_ext, const std::vector<double>& lambda, bool stub = false);

    void append_Ku_frictionless(const Meshing* const mesh, const std::vector<double>& u, std::vector<double>& Ku, const std::vector<math::Matrix>& D_cache, bool topopt) const;
    void append_dKu_frictionless(const Meshing* const mesh, const std::vector<double>& u, const std::vector<double>& du, const double eta, std::vector<double>& Ku, const std::vector<math::Matrix>& D_cache, bool topopt) const;

    void append_Ku_frictionless_simple(const Meshing* const mesh, const std::vector<double>& u, std::vector<double>& Ku) const;
    void append_dKu_frictionless_simple(const Meshing* const mesh, const std::vector<double>& u, const std::vector<double>& du, const double eta, std::vector<double>& Ku) const;

    inline void set_lag_displ_simple(double L){
        this->LAG_DISPL_SIMPLE = L;
    }
    inline double get_lag_displ_simple() const{
        return this->LAG_DISPL_SIMPLE;
    }

    protected:
    size_t W, N;
    const double EPS_PENALTY = 5e7;
    const double EPS_DISPL;
    double LAG_DISPL_SIMPLE = 1e4;
    bool first_time = true;

    virtual void final_flush_matrix(){}

    virtual void generate_base(const Meshing * const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const std::vector<double>& lambda, const FiniteElement::ContactType type);

    virtual void calculate_dimensions(const Meshing * const mesh, const std::vector<long>& node_positions, const size_t matrix_width);

    virtual void add_geometry(const Meshing * const mesh, const std::vector<long>& node_positions, const Geometry * const g);

    virtual void add_geometry(const Meshing * const mesh, const std::vector<long>& node_positions, const Geometry * const g, const size_t D_offset, const std::vector<math::Matrix>& D_cache);

    virtual void add_springs(const Meshing * const mesh, const std::vector<long>& node_positions);

    virtual void add_contacts(const Meshing * const mesh, const std::vector<long>& node_positions, const std::vector<double>& u_ext);

    virtual void insert_block_symmetric(const math::Matrix& k, const std::vector<long>& posi, const std::vector<long>& posj) = 0;

    virtual void insert_element_matrix(const math::Matrix& k, const std::vector<long>& pos) = 0;

    virtual void add_to_matrix(size_t i, size_t j, double val) = 0;
};

#endif
