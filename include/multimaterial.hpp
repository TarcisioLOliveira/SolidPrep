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

#ifndef MULTIMATERIAL_HPP
#define MULTIMATERIAL_HPP

#include "material.hpp"
#include "utils.hpp"

class MultiMaterial{
    public:
    MultiMaterial() = default;
    MultiMaterial(std::vector<Material*> materials, utils::ProblemType type, bool has_void);

    inline std::vector<double> get_D(const MeshElement* e, const gp_Pnt& p) const{
        if(this->problem_type == utils::PROBLEM_TYPE_2D){
            return this->materials[0]->stiffness_2D(e, p);
        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
            return this->materials[0]->stiffness_3D(e, p);
        }
        return std::vector<double>();
    }

    void get_D(std::vector<double>::const_iterator rho, const double mix, const MeshElement* e, const gp_Pnt& p, std::vector<double>& D) const;

    void get_gradD(std::vector<double>::const_iterator rho, const double mix, const MeshElement* e, const gp_Pnt& p, std::vector<std::vector<double>>& gradD) const;

    double get_density(std::vector<double>::const_iterator rho, const MeshElement* e, const gp_Pnt& p) const;

    double get_density_deriv(std::vector<double>::const_iterator rho, const MeshElement* e, const gp_Pnt& p, std::vector<double>::iterator& grad) const;

    inline size_t number_of_materials() const{
        return this->materials.size();
    }

    inline const std::vector<Material*>& get_materials() const{
        return this->materials;
    }

    private:
    std::vector<Material*> materials;
    utils::ProblemType problem_type;
    bool has_void;

    void get_D_internal(std::vector<double>::const_iterator& rho, const double mix, const MeshElement* e, const gp_Pnt& p, std::vector<double>& D) const;
    void get_gradD_internal(std::vector<double>::const_iterator& rho, const double mix, const MeshElement* e, const gp_Pnt& p, std::vector<std::vector<double>>& gradD) const;
};

#endif
