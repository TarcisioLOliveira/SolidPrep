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

#ifndef GT9_HPP
#define GT9_HPP

#include "element.hpp"
#include "material.hpp"
#include "utils.hpp"

class ProjectData;

namespace element{

class GT9 : public MeshElement{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 2;
    static const size_t NODE_DOF       = 3;
    static const size_t NODES_PER_ELEM = 3;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;

    static const utils::ProblemType PROBLEM_TYPE = utils::PROBLEM_TYPE_2D;

    GT9(ElementShape s, ProjectData* data);

    virtual std::vector<double> get_k() const override;
    virtual std::vector<double> get_internal_loads(size_t node, const std::vector<double>& u) const override;
    virtual double get_stress_at(gp_Pnt p, const std::vector<double>& u) const override;
    virtual std::vector<double> get_stress_tensor(gp_Pnt p, const std::vector<double>& u) const override;
    virtual std::vector<double> get_loads_at(gp_Pnt p, const std::vector<double>& u) const override;
    virtual std::vector<gp_Pnt> get_intersection_points(const TopoDS_Shape& crosssection) const override;
    virtual double get_compliance(const std::vector<double>& u, const std::vector<double>& l = std::vector<double>()) const override;
    virtual double get_volume() const override;
    virtual void get_virtual_load(double mult, gp_Pnt point, const std::vector<double>& u, std::vector<double>& l) const override;
    virtual std::vector<double> get_f(gp_Dir dir, double norm, std::vector<gp_Pnt> points) const override;
    virtual TopoDS_Shape get_shape(std::vector<gp_Vec> disp = std::vector<gp_Vec>()) const override;
    virtual gp_Pnt get_centroid() const override;

    private:
    Material const * const mat;
    double t;

    std::vector<double> get_DB(gp_Pnt point) const;
};

}

#endif
