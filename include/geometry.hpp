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

#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <string>
#include <gp_Pnt.hxx>
#include <TopoDS_Shape.hxx>
#include "element.hpp"
#include "utils.hpp"
#include "material.hpp"
#include <BRepClass3d_SolidClassifier.hxx>

class Geometry{
    public:
    Geometry(const std::string& path, double scale, utils::ProblemType type, MeshElementFactory* elem_type, bool do_topopt, Material* material, std::vector<Material*> alt_materials = std::vector<Material*>());

    inline bool is_inside(const gp_Pnt& p) const{
        if(this->type == utils::PROBLEM_TYPE_2D){
            return this->is_inside_2D(p);
        } else if(this->type == utils::PROBLEM_TYPE_3D){
            return this->is_inside_3D(p);
        }
        return false;
    }

    inline std::vector<double> get_D() const{
        if(this->type == utils::PROBLEM_TYPE_2D){
            return this->material->stiffness_2D();
        } else if(this->type == utils::PROBLEM_TYPE_3D){
            return this->material->stiffness_3D();
        }
        return std::vector<double>();
    }

    const TopoDS_Shape shape;
    const double scale;
    const Material* const material;
    const std::vector<const Material*> alternate_materials;
    const MeshElementFactory* const element_type;
    const bool do_topopt;

    std::vector<MeshElement*> mesh;
    private:
    utils::ProblemType type;

    TopoDS_Shape load_shape(const std::string& path, double scale) const;

    inline bool is_inside_2D(const gp_Pnt& p) const{
        BRepClass3d_SolidClassifier insider(this->shape, p, 0.01);
        return insider.State() == TopAbs_ON;
    }

    inline bool is_inside_3D(const gp_Pnt& p) const{
        BRepClass3d_SolidClassifier insider(this->shape, p, 0.01);
        return insider.State() == TopAbs_IN;
    }
};

#endif
