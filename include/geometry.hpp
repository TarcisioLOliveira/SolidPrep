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

class Meshing;

class Geometry{
    public:
    Geometry(const std::string& path, double scale, utils::ProblemType type, MeshElementFactory* elem_type, bool do_topopt, bool with_void, Material* material, std::vector<Material*> alt_materials = std::vector<Material*>());

    Geometry(TopoDS_Shape shape, utils::ProblemType type, MeshElementFactory* elem_type, bool do_topopt, bool with_void, Material* material, std::vector<Material*> alt_materials = std::vector<Material*>());

    void get_stresses(const std::vector<double>& u, std::vector<double>::iterator& stress_it) const;

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
        return this->constitutive_matrices.size();
    }

    inline size_t number_of_densities_needed() const{
        if(this->alternative_materials.size() == 0){
            return 1;
        } else {
            if(with_void){
                return this->alternative_materials.size() + 1;
            } else {
                return this->alternative_materials.size();
            }
        }
    }

    /**
     * Returns the constitutive matrix for the specified material.
     * i = 0: default material
     * i > 0: alternative materials (for topology optimization)
     *
     * @param i Material
     *
     * @return Selected constitutive matrix.
     */
    inline std::vector<double> get_D(const size_t i) const{
        return this->constitutive_matrices[i];
    }


    /**
     * Single density.
     */
    inline std::vector<double> get_D_topopt(const double p, const double pc) const{
        if(this->constitutive_matrices.size() == 1){
            auto D = this->constitutive_matrices[0];
            const double p_pen = std::pow(p, pc);
            for(auto& d:D){
                d *= p_pen;
            }
            return D;
        } else {
            return this->D_topopt_mult({p}, pc, 1);
        }
    }

    /**
     * Single density, using minimum density.
     */
    inline std::vector<double> get_D_topopt(const double p, const double pc, const double p_min) const{
        if(this->constitutive_matrices.size() == 1){
            auto D = this->constitutive_matrices[0];
            const double p_norm = p_min + (1-p_min)*std::pow(p, pc);
            for(auto& d:D){
                d *= p_norm;
            }
            return D;
        } else {
            return this->D_topopt_mult({p}, pc, 1);
        }
    }

    /**
     * Multiple densities.
     */
    inline std::vector<double> get_D_topopt(const std::vector<double>& p, const double pc) const{
        if(!this->with_void){
            return this->D_topopt_mult(p, pc, p.size());
        } else {
            const size_t p_size = p.size();
            auto D = this->D_topopt_mult(p, pc, p_size-1);
            const size_t i = p_size - 1;
            const double p_penv = std::pow(p[i], pc);
            for(auto& d:D){
                d *= p_penv;
            }
            return D;
        }
    }

    /**
     * Multiple densities with minimum density.
     */
    inline std::vector<double> get_D_topopt(const std::vector<double>& p, const double pc, const double p_min) const{
        if(!this->with_void){
            return this->D_topopt_mult(p, pc, p.size());
        } else {
            const size_t p_size = p.size();
            auto D = this->D_topopt_mult(p, pc, p_size-1);
            const size_t i = p_size - 1;
            const double p_penv = p_min + (1-p_min)*std::pow(p[i], pc);
            for(auto& d:D){
                d *= p_penv;
            }
            return D;
        }
    }

    TopoDS_Shape shape;
    const Material* const material;
    const std::vector<const Material*> alternative_materials;
    const MeshElementFactory* const element_type;
    const bool do_topopt;
    const bool with_void;
    const std::vector<std::vector<double>> constitutive_matrices;

    std::vector<std::unique_ptr<MeshElement>> mesh;
    // Only used for non-linear FEA
    std::vector<std::unique_ptr<MeshNode>> node_list; 
    private:
    const utils::ProblemType type;

    inline std::vector<std::vector<double>> init_constitutive_matrices(const utils::ProblemType type) const{
        std::vector<std::vector<double>> D;
        D.reserve(1 + this->alternative_materials.size());
        
        if(type == utils::PROBLEM_TYPE_2D){
            D.push_back(material->stiffness_2D());
            for(const auto& m:alternative_materials){
                D.push_back(m->stiffness_2D());
            }
        } else if(type == utils::PROBLEM_TYPE_3D){
            D.push_back(material->stiffness_3D());
            for(const auto& m:alternative_materials){
                D.push_back(m->stiffness_3D());
            }
        }
        return D;
    }

    inline std::vector<double> D_topopt_mult(const std::vector<double>& p, const double pc, size_t p_size) const{
        // TODO: make this actually work for 3 or more materials.
        size_t i = 0;
        std::vector<double> D(this->constitutive_matrices[i].size(),0);
        for(i = 0; i < p_size; ++i){
            const auto& Di = this->constitutive_matrices[i];
            const double p_peni = std::pow(p[i], pc);
            for(size_t j = 0; j < Di.size(); ++j){
                D[j] += Di[j]*p_peni;
            }
        }
        i = p_size - 1;
        const auto& Di = this->constitutive_matrices[i+1];
        const double p_peni = std::pow(1 - p[i], pc);
        for(size_t j = 0; j < Di.size(); ++j){
            D[j] += Di[j]*p_peni;
        }
        return D;
    }

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
