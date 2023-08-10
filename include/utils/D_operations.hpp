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

#ifndef UTILS_D_OPERATIONS_HPP
#define UTILS_D_OPERATIONS_HPP

#include <vector>

namespace utils::D_op{

std::vector<double> invert_2D(const std::vector<double>& d);
std::vector<double> invert_3D(const std::vector<double>& d);
std::vector<double> square_2D(const std::vector<double>& d);
std::vector<double> square_3D(const std::vector<double>& d);
std::vector<double> mult_2D(const std::vector<double>& d, const std::vector<double>& e);
std::vector<double> mult_3D(const std::vector<double>& d, const std::vector<double>& e);
std::vector<double> square_root_2D(const std::vector<double>& d);
std::vector<double> square_root_3D(const std::vector<double>& d);

}

#endif
