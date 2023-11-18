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

#ifndef UTILS_LIN3B_HPP
#define UTILS_LIN3B_HPP

#include <cstddef>
#include <Eigen/Core>

namespace utils{

class LIN3B{
    public:
    LIN3B(std::array<double, 3> Y);

    Eigen::Matrix<double, 3, 3> absorption() const;
    Eigen::Vector<double, 3> source(std::function<double(double)> fn) const;

    double N(double y, size_t i) const{
        return a[i] + b[i]*y + c[i]*y*y;
    }
    double dN(double y, size_t i) const{
        return b[i] + 2*c[i]*y;
    }
    std::array<double, 3> Y;

    private:
    double a[3], b[3], c[3];
    double len;

    Eigen::Vector<double, 3> NN(double y) const{
        return{N(y,0), N(y,1), N(y,2)};
    }
    Eigen::Vector<double, 3> dNN(double y) const{
        return{dN(y,0), dN(y,1), dN(y,2)};
    }
};

}

#endif
