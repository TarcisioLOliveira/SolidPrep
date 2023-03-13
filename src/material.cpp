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

#include "material.hpp"
#include <cmath>

Material::Material(const std::string& name, const double density, std::vector<double> Smax, std::vector<double> Tmax):
    name(name), density(density){
    if(Smax.size() == 1){
        this->Smax = std::vector<double>(6, Smax[0]);
    } else if(Smax.size() == 3){
        this->Smax = std::move(Smax);
        this->Smax.insert(Smax.end(), Smax.begin(), Smax.end());
    } else {
        this->Smax = std::move(Smax);
    }
    if(Tmax.size() == 1){
        this->Tmax = std::vector<double>(3, Tmax[0]);
    } else {
        this->Tmax = std::move(Tmax);
    }
}


double Material::get_max_Von_Mises_2D() const{
    double S11 = std::min(this->Smax[0], this->Smax[3]);
    double S22 = std::min(this->Smax[1], this->Smax[4]);
    double T12 = this->Tmax[0];
    return std::sqrt(0.5*(std::pow(S11-S22, 2)
                          + std::pow(S11, 2) + std::pow(S22,2)
                          + 6*std::pow(T12, 2)));
}
double Material::get_max_Von_Mises_3D() const{
    double S11 = std::min(this->Smax[0], this->Smax[3]);
    double S22 = std::min(this->Smax[1], this->Smax[4]);
    double S33 = std::min(this->Smax[2], this->Smax[5]);
    double T12 = this->Tmax[0];
    double T13 = this->Tmax[1];
    double T23 = this->Tmax[2];
    return std::sqrt(0.5*(std::pow(S11-S22, 2)
                          + std::pow(S22-S33, 2)
                          + std::pow(S11-S33, 2)
                          + 6*(std::pow(T12, 2)
                               + std::pow(T13, 2)
                               + std::pow(T23, 2))));
}
