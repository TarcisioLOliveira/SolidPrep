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

#ifndef GT9_HPP
#define GT9_HPP

#include "element.hpp"
#include "material.hpp"

class ProjectData;

namespace element{

class GT9 : public MeshElement{
    public:
    GT9(ElementShape s, ProjectData* data);

    virtual std::vector<float> get_k() const override;
    virtual MeshNode* get_stresses(size_t node, const std::vector<float>& u, double density = 1) const override;
    virtual MeshNode* get_internal_loads(size_t node, const std::vector<float>& u) const override;
    virtual double get_stress_at(gp_Pnt p, const std::vector<float>& u) const override;
    virtual size_t get_gmsh_element_type() const override{ return 2;};
    virtual double get_compliance(const std::vector<float>& u, const std::vector<float>& l = std::vector<float>()) const override;
    virtual double get_volume() const override;
    virtual void get_virtual_load(double P, gp_Pnt point, std::vector<float>& u, std::vector<float>& l) const override;
    virtual TopoDS_Shape get_shape() const override;
    virtual gp_Pnt get_centroid() const override;

    private:
    Material const * const mat;
    float t;

    std::vector<float> get_DB(gp_Pnt point) const;
};

}

#endif
