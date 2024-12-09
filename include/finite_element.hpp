/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
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

#ifndef FINITE_ELEMENT_HPP
#define FINITE_ELEMENT_HPP

#include <vector>
#include "math/matrix.hpp"
#include "meshing.hpp"
#include "geometry.hpp"

class GlobalStiffnessMatrix;

class FiniteElement{
    public:
    const double K_MIN = 1e-6;

    enum ContactType{
        RIGID,
        FRICTIONLESS_PENALTY,
        FRICTIONLESS_DISPL
    };

    FiniteElement(ContactType contact_type, double rtol_abs, GlobalStiffnessMatrix* m);

    virtual ~FiniteElement() = default;

    void generate_matrix(const Meshing* const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext);

    void calculate_displacements(const Meshing* const mesh, std::vector<double>& load, const std::vector<double>& u0, std::vector<double>& lambda, const bool topopt, const std::vector<math::Matrix>& D_cache);

    virtual std::vector<double> calculate_forces(const Meshing* const mesh, const std::vector<double>& displacements) const;

    protected:
    const ContactType contact_type;
    const double rtol_abs;
    GlobalStiffnessMatrix* matrix = nullptr;
    size_t u_size, l_num;

    virtual void generate_matrix_base(const Meshing* const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u_ext, const ContactType type) = 0;

    virtual void solve(std::vector<double>& load) = 0;
    virtual void reset_hessian() = 0;

    void solve_rigid(std::vector<double>& load);
    void solve_frictionless_displ(const Meshing* const mesh, std::vector<double>& load, std::vector<double>& lambda);

    void solve_frictionless_penalty(const Meshing* const mesh, std::vector<double>& load, const bool topopt, const std::vector<math::Matrix>& D_cache, const std::vector<double>& u0);

    private:
};

#endif
