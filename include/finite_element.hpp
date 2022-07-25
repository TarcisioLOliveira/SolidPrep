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

#include <memory>
#include <vector>
#include "element.hpp"
#include "meshing.hpp"
#include "element_factory.hpp"
#include "geometry.hpp"

class ProjectData;

class FiniteElement{
    public:
    const double K_MIN = 1e-6;

    virtual std::vector<double> calculate_displacements(Meshing* mesh, std::vector<double> load, const std::vector<double>& density = std::vector<double>(), double pc = 3) = 0;

    inline virtual void set_steps(size_t s){
        this->steps = s;
    }

    virtual std::vector<double> calculate_forces(const Meshing* const mesh, const std::vector<double>& displacements) const;

    protected:
    size_t steps = 1;
    std::vector<double> K = std::vector<double>();
    size_t W = 0;
    size_t N = 0;
    size_t current_step = 0;
    bool recalculated_dimensions = false;

    virtual void calculate_dimensions(const MeshElementFactory* const element, Meshing* mesh, const std::vector<double>& load);

    virtual void generate_K(Meshing* mesh, const std::vector<double>& density, const double pc);

    virtual void add_geometry_to_K(Meshing* mesh, Geometry* g);

    virtual void add_geometry_to_K(Meshing* mesh, Geometry* g, std::vector<double>::const_iterator rho, const double pc);

    virtual void insert_element_matrix(std::vector<double>& K, const std::vector<double>& k, const std::vector<long>& pos, size_t n) const;
};

#endif
