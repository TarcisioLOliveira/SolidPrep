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
    virtual std::vector<double> calculate_displacements(ProjectData* data, Meshing* mesh, const std::vector<double>& density = std::vector<double>(), double pc = 3, bool use_stored_matrix = false, const std::vector<double>& virtual_load = std::vector<double>()) override;

    private:
    std::vector<double> K = std::vector<double>();
    size_t W = 0;
    size_t N = 0;
    std::vector<double> calculate_displacements_simple(ProjectData* data, Meshing* mesh, bool use_stored_matrix);

    void insert_element_matrix(std::vector<double>& K, const std::vector<double>& k, const std::vector<long>& pos, size_t w, size_t n) const;
};

}

#endif
