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
#include "multimaterial.hpp"
#include "utils.hpp"
#include "material.hpp"
#include <BRepClass3d_SolidClassifier.hxx>

class Meshing;
class BoundaryElement;

class Geometry{
    public:
    Geometry(const std::string& path, double scale, utils::ProblemType type, MeshElementFactory* elem_type, bool do_topopt, bool with_void, std::vector<Material*> materials, size_t id);

    Geometry(TopoDS_Shape shape, utils::ProblemType type, MeshElementFactory* elem_type, bool do_topopt, bool with_void, std::vector<Material*> materials, size_t id);

    void get_stresses(const std::vector<double>& u, const double p, const double psi, std::vector<double>::const_iterator& rho_it, std::vector<double>::iterator& stress_it) const;

    /**
     * Check if point is inside the geometry.
     * Very slow, avoid using it a lot.
     *
     * @param p Point
     *
     * @return Whether the point is inside the geometry or not.
     */
    inline bool is_inside(const gp_Pnt& p) const{
        if(this->type == utils::PROBLEM_TYPE_2D){
            return this->is_inside_2D(p);
        } else if(this->type == utils::PROBLEM_TYPE_3D){
            return this->is_inside_3D(p);
        }
        return false;
    }

    /**
     * Returns the number of available materials for use in topology
     * optimization.
     *
     * @return Number of available materials.
     */
    inline size_t number_of_materials() const{
        return this->materials.number_of_materials();
    }

    inline size_t number_of_densities_needed() const{
        if(this->materials.number_of_materials() == 1){
            return 1;
        } else {
            if(with_void){
                return this->materials.number_of_materials();
            } else {
                return this->materials.number_of_materials() - 1;
            }
        }
    }

    TopoDS_Shape shape;
    const MultiMaterial materials;
    const MeshElementFactory* const element_type;
    const bool do_topopt;
    const bool with_void;
    const size_t id;

    std::vector<std::unique_ptr<MeshElement>> mesh;
    std::vector<BoundaryElement*> boundary_mesh;
    // Only used for non-linear FEA
    std::vector<std::unique_ptr<MeshNode>> node_list; 
    private:
    const utils::ProblemType type;

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
