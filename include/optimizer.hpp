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

#ifndef OPTIMIZER_HPP
#define OPTIMIZER_HPP

#include <BRepBuilderAPI_Copy.hxx>
#include <TopoDS_Shape.hxx>
#include "finite_element.hpp"
#include "visualization.hpp"

class DensityBasedFunction;

class Constraint{
    public:
    enum class Type{
        LESS_THAN,
        GREATER_THAN,
        EQUAL
    };

    Constraint(std::unique_ptr<DensityBasedFunction> fun, std::vector<Type> types, std::vector<double> bounds);

    std::unique_ptr<DensityBasedFunction> fun;
    std::vector<Type> types;
    std::vector<double> bounds;
};

class Optimizer{
    public:
    const double K_MIN = 1e-6;

    virtual ~Optimizer() = default;

    virtual void initialize_views(Visualization* viz) = 0;
    virtual TopoDS_Shape optimize(FiniteElement* fem, Meshing* mesh) = 0;

    inline size_t get_number_of_elements() const{
        return this->number_of_elements;
    }
    inline const std::vector<double>& get_volumes() const{
        return this->volumes;
    }
    inline const std::vector<double>& get_stresses() const{
        return this->stresses;
    }
    inline const std::vector<double>& get_filtered_densities() const{
        return this->filtered_densities;
    }

    protected:
    size_t number_of_elements;
    std::vector<double> volumes;
    std::vector<double> stresses;
    std::vector<double> filtered_densities;

    void initialize_optimizer(const Meshing* const mesh);

    TopoDS_Shape make_shape(const std::vector<double>& x, const std::vector<Geometry*>& geometries, const double result_threshold) const;
    void get_stresses(const std::vector<Geometry*> geometries, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& stresses, double pc, double psi) const;
    void get_volumes(const std::vector<Geometry*> geometries, const double thickness, std::vector<double>& volumes) const;
    size_t get_number_of_elements(const std::vector<Geometry*> geometries) const;

    void apply_densities(const std::vector<Geometry*> geometries, const std::vector<double>& x, std::vector<double>& vals, const double pc = 1.0) const;
};


#endif
