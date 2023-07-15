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

#ifndef GAUSS_LEGENDRE_HPP
#define GAUSS_LEGENDRE_HPP

#include <array>
#include <cstddef>
#include <gsl/gsl_integration.h>

namespace utils{

struct GLPoint{
    double x, w;
};

template<size_t N>
class GaussLegendre{
    public:
    GaussLegendre(const GaussLegendre& g) = delete;
    GaussLegendre(GaussLegendre&& g) = delete;
    void operator=(const GaussLegendre& g) = delete;
    void operator=(GaussLegendre&& g) = delete;

    static const GaussLegendre& get(){
        static GaussLegendre instance;

        return instance;
    }

    auto begin() const{
        return p.cbegin();
    }

    auto end() const{
        return p.cend();
    }

    private:
    GaussLegendre(){
        auto table = gsl_integration_glfixed_table_alloc(N);

        for(size_t i = 0; i < N; ++i){
            gsl_integration_glfixed_point(-1, 1, i, &p[i].x, &p[i].w, table);
        }

        gsl_integration_glfixed_table_free(table);
    }

    std::array<GLPoint, N> p;
};

}

#endif
