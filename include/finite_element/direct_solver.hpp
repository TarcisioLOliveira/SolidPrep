/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   Foobar is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Foobar is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with SolidPrep.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef DIRECT_SOLVER_HPP
#define DIRECT_SOLVER_HPP

#include "finite_element.hpp"

namespace finite_element{

class DirectSolver : public FiniteElement{
    public:
    virtual std::vector<double> calculate_displacements(ProjectData* data, Meshing* mesh, const std::vector<double>& density = std::vector<double>(), double pc = 3, const std::vector<double>& virtual_load = std::vector<double>()) const override;

    private:
    std::vector<double> calculate_displacements_simple(ProjectData* data, Meshing* mesh) const;

    void insert_element_matrix(std::vector<double>& K, const std::vector<double>& k, const std::vector<long>& pos, size_t w, size_t n) const;
    std::vector<double> expand_U(const std::vector<long>& new_pos, std::vector<double>& U, Meshing* mesh, size_t dof, size_t size) const;
};

}

#endif
