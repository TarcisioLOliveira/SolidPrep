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

#ifndef TRI3_HPP
#define TRI3_HPP

#include "element.hpp"
#include "material.hpp"
#include <vector>
#include "utils.hpp"
#include "element_factory.hpp"

class ProjectData;

namespace element{

class TRI3 : public MeshElement{
    public:
    static const size_t ORDER          = 1;
    static const size_t GMSH_TYPE      = 2;
    static const size_t NODE_DOF       = 2;
    static const size_t NODES_PER_ELEM = 3;
    static const size_t K_DIM          = NODE_DOF*NODES_PER_ELEM;

    static const utils::ProblemType PROBLEM_TYPE = utils::PROBLEM_TYPE_2D;

    TRI3(ElementShape s, ProjectData* data);

    virtual std::vector<double> get_k() const override;
    virtual std::vector<double> get_internal_loads(const std::vector<double>& u) const override;
    virtual std::vector<gp_Pnt> get_intersection_points(const TopoDS_Shape& crosssection) const override;
    virtual double get_stress_at(gp_Pnt p, const std::vector<double>& u) const override;
    virtual std::vector<double> get_stress_tensor(const gp_Pnt& p, const std::vector<double>& u) const override;
    virtual double get_compliance(const std::vector<double>& u, const std::vector<double>& l = std::vector<double>()) const override;
    virtual void get_virtual_load(double mult, const gp_Pnt& point, const std::vector<double>& u, std::vector<double>& l) const override;
    virtual std::vector<double> get_f(const gp_Dir& dir, double norm, const std::vector<gp_Pnt>& points) const override;
    virtual double get_volume() const override;
    virtual TopoDS_Shape get_shape(const std::vector<gp_Vec>& disp = std::vector<gp_Vec>()) const override;
    virtual gp_Pnt get_centroid() const override;
    virtual inline std::unique_ptr<MeshElementFactory> get_element_info() const override{
        return std::unique_ptr<MeshElementFactory>(new MeshElementFactoryImpl<TRI3>());
    }

    private:
    Material const * const mat;
    double t;

    std::vector<double> get_DB(const gp_Pnt& point) const;
};

}

#endif
