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
    MultiMaterial(std::vector<Material*> materials, utils::ProblemType type);

    inline std::vector<double> get_D() const{
        if(this->problem_type == utils::PROBLEM_TYPE_2D){
            return this->materials[0]->stiffness_2D();
        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
            return this->materials[0]->stiffness_3D();
        }
        return std::vector<double>();
    }

    void get_D(std::vector<double>::const_iterator& rho, const bool has_void, const double p, const double p_min, const double mix, std::vector<double>& D) const;

    void get_gradD(std::vector<double>::const_iterator& rho, const bool has_void, const double p, const double p_min, const double mix, std::vector<std::vector<double>>& gradD) const;

    double get_density(std::vector<double>::const_iterator& rho, const bool has_void) const;

    double get_density_deriv(std::vector<double>::const_iterator& rho, const bool has_void, std::vector<double>::iterator& grad) const;

    inline size_t number_of_materials() const{
        return this->materials.size();
    }

    inline const std::vector<Material*>& get_materials() const{
        return this->materials;
    }

    private:
    const std::vector<Material*> materials;
    const utils::ProblemType problem_type;

    void get_D_internal(std::vector<double>::const_iterator& rho, const double* pos, const size_t posN, const double mix, std::vector<double>& D) const;
    void get_gradD_internal(std::vector<double>::const_iterator& rho, const bool has_void, const double* pos, const size_t posN, const double mix, std::vector<std::vector<double>>& gradD) const;

    inline void apply_void(const double rho, const double* pos, const size_t posN, std::vector<double>& D) const{
        for(size_t i = 0; i < posN; ++i){
            D[pos[i]] *= rho;
        }
    }

    std::vector<double> invert_2D(const std::vector<double>& d) const;
    std::vector<double> invert_3D(const std::vector<double>& d) const;
    std::vector<double> square_2D(const std::vector<double>& d) const;
    std::vector<double> square_3D(const std::vector<double>& d) const;
    std::vector<double> mult_2D(const std::vector<double>& d, const std::vector<double>& e) const;
    std::vector<double> mult_3D(const std::vector<double>& d, const std::vector<double>& e) const;
    std::vector<double> square_root_2D(const std::vector<double>& d) const;
    std::vector<double> square_root_3D(const std::vector<double>& d) const;
};

#endif
