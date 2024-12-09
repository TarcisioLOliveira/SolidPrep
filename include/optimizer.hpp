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
#include "math/matrix.hpp"
#include "shape_handler.hpp"
#include "utils.hpp"
#include "visualization.hpp"
#include "solver_manager.hpp"
#include "function.hpp"

class DensityBasedConstraint;

class Constraint{
    public:
    enum class Type{
        LESS_THAN,
        GREATER_THAN,
        EQUAL
    };

    Constraint() = default;
    Constraint(const Constraint&) = default;
    Constraint(Constraint&&) = default;
    Constraint& operator=(const Constraint&) = default;
    Constraint& operator=(Constraint&&) = default;

    Constraint(std::vector<Type> types, std::vector<double> bounds);
    virtual ~Constraint() = default;

    std::vector<Type> types;
    std::vector<double> bounds;
};

class DensityBasedConstraint : public Constraint{
    public:
    DensityBasedConstraint() = default;
    DensityBasedConstraint(const DensityBasedConstraint&) = delete;
    DensityBasedConstraint(DensityBasedConstraint&&) = default;
    DensityBasedConstraint& operator=(const DensityBasedConstraint&) = delete;
    DensityBasedConstraint& operator=(DensityBasedConstraint&&) = default;

    DensityBasedConstraint(std::unique_ptr<DensityBasedFunction> fun, std::vector<Type> types, std::vector<double> bounds);
    virtual ~DensityBasedConstraint() = default;

    std::unique_ptr<DensityBasedFunction> fun;
};

class NodeShapeBasedConstraint : public Constraint{
    public:
    NodeShapeBasedConstraint() = default;
    NodeShapeBasedConstraint(const NodeShapeBasedConstraint&) = delete;
    NodeShapeBasedConstraint(NodeShapeBasedConstraint&&) = default;
    NodeShapeBasedConstraint& operator=(const NodeShapeBasedConstraint&) = delete;
    NodeShapeBasedConstraint& operator=(NodeShapeBasedConstraint&&) = default;


    NodeShapeBasedConstraint(std::unique_ptr<NodeShapeBasedFunction> fun, std::vector<Type> types, std::vector<double> bounds);
    virtual ~NodeShapeBasedConstraint() = default;

    std::unique_ptr<NodeShapeBasedFunction> fun;
};

class Optimizer{
    public:

    virtual ~Optimizer() = default;

    virtual void initialize_views(Visualization* viz) = 0;
    virtual TopoDS_Shape optimize(SolverManager* fem, Meshing* mesh) = 0;

    inline size_t get_number_of_elements() const{
        return this->number_of_elements;
    }
    inline const std::vector<double>& get_volumes() const{
        return this->volumes;
    }
    inline const std::vector<double>& get_stresses() const{
        return this->stresses;
    }

    protected:
    size_t number_of_elements;
    std::vector<double> volumes;
    std::vector<double> stresses;

    void initialize_optimizer(const Meshing* const mesh);

    void get_stresses(const std::vector<Geometry*> geometries, const bool is_topopt, const std::vector<double>& u, const std::vector<math::Matrix>& D_cache, std::vector<double>& stresses) const;
    void get_volumes(const std::vector<Geometry*> geometries, const double thickness, std::vector<double>& volumes) const;
    size_t get_number_of_elements(const std::vector<Geometry*> geometries) const;

    // Workaround for make_shape with 3D shapes because for some reason cutting
    // does not work when using 3D shapes except if I turn it into STEP and
    // back into TopoDS_Shape.
    TopoDS_Shape STEP_workaround(const TopoDS_Shape& s) const;

};

class DensityBasedOptimizer : public Optimizer{
    public:
    virtual ~DensityBasedOptimizer() = default;

    inline const std::vector<double>& get_filtered_densities() const{
        return this->filtered_densities;
    }

    protected:
    std::vector<double> filtered_densities;

    TopoDS_Shape make_shape(const std::vector<double>& x, const std::vector<Geometry*>& geometries, const double result_threshold, const utils::ProblemType type) const;

    void apply_densities(const std::vector<Geometry*> geometries, const std::vector<double>& x, std::vector<double>& vals, const double pc = 1.0) const;
};

class NodeShapeBasedOptimizer : public Optimizer{
    public:
    NodeShapeBasedOptimizer(ShapeHandler sh):
        shape_handler(std::move(sh)){}

    ShapeHandler shape_handler;

    virtual ~NodeShapeBasedOptimizer() = default;

    protected:

    TopoDS_Shape make_shape(const std::vector<Geometry*>& geometries, const utils::ProblemType type) const;
};


#endif
