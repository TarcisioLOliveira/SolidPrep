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

#include "material.hpp"

Material::Material(std::vector<double> Smax, std::vector<double> Tmax){
    if(Smax.size() == 1){
        this->Smax = std::vector<double>(3, Smax[0]);
    } else {
        this->Smax = std::move(Smax);
    }
    if(Tmax.size() == 1){
        this->Tmax = std::vector<double>(3, Tmax[0]);
    } else {
        this->Tmax = std::move(Tmax);
    }
}
