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
#include "logger.hpp"

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


struct GLPointTri{
    double a, b, c, w;
};

namespace internal{

template<size_t N>
class GaussLegendreTriSuper{
    public:
    GaussLegendreTriSuper(const GaussLegendreTriSuper& g) = delete;
    GaussLegendreTriSuper(GaussLegendreTriSuper&& g) = delete;
    void operator=(const GaussLegendreTriSuper& g) = delete;
    void operator=(GaussLegendreTriSuper&& g) = delete;

    //static const GaussLegendreTriSuper& get(){
    //    static GaussLegendreTriSuper instance;

    //    return instance;
    //}

    auto begin() const{
        return p.cbegin();
    }

    auto end() const{
        return p.cend();
    }

    protected:
    GaussLegendreTriSuper(std::array<GLPointTri, N> p): 
        p(std::move(p))
    {
    }

    std::array<GLPointTri, N> p;
};

}

constexpr size_t P_to_N(size_t P){
    switch(P){
        case 1: return 1;
        case 2: return 3;
        case 3: return 4;
        case 4: return 6;
        case 5: return 7;
        case 6: return 12;
        case 7: return 13;
    }
    return 0;
}

template<size_t P>
class GaussLegendreTri : public internal::GaussLegendreTriSuper<P_to_N(P)>{
    public:
    GaussLegendreTri(const GaussLegendreTri& g) = delete;
    GaussLegendreTri(GaussLegendreTri&& g) = delete;
    void operator=(const GaussLegendreTri& g) = delete;
    void operator=(GaussLegendreTri&& g) = delete;

    static const GaussLegendreTri& get(){
        logger::log_assert(false, logger::ERROR, "Gauss Legendre for triangles not implemented for polynomials of order {}.", P);
        static GaussLegendreTri instance;

        return instance;
    }

    private:
    GaussLegendreTri()
    : internal::GaussLegendreTriSuper<P_to_N(P)>(std::array<GLPointTri, P_to_N(P)>())
    {
    }
};


// (DUNAVANT, 1985)
template<>
class GaussLegendreTri<1> : public internal::GaussLegendreTriSuper<1>{
    public:
    GaussLegendreTri(const GaussLegendreTri& g) = delete;
    GaussLegendreTri(GaussLegendreTri&& g) = delete;
    void operator=(const GaussLegendreTri& g) = delete;
    void operator=(GaussLegendreTri&& g) = delete;

    static const GaussLegendreTri& get(){
        static GaussLegendreTri instance;

        return instance;
    }

    private:
    GaussLegendreTri()
    : internal::GaussLegendreTriSuper<1>(std::array<GLPointTri, 1>
            {{{1.0/3.0, 1.0/3.0, 1.0/3.0, 1}}}
        )
    {
    }
};
template<>
class GaussLegendreTri<2> : public internal::GaussLegendreTriSuper<3>{
    public:
    GaussLegendreTri(const GaussLegendreTri& g) = delete;
    GaussLegendreTri(GaussLegendreTri&& g) = delete;
    void operator=(const GaussLegendreTri& g) = delete;
    void operator=(GaussLegendreTri&& g) = delete;

    static const GaussLegendreTri& get(){
        static GaussLegendreTri instance;

        return instance;
    }

    private:
    GaussLegendreTri()
    : internal::GaussLegendreTriSuper<3>(std::array<GLPointTri, 3>
            {{
                {0.666666666666667, 0.166666666666667, 0.166666666666667, 0.333333333333333},
                {0.166666666666667, 0.666666666666667, 0.166666666666667, 0.333333333333333},
                {0.166666666666667, 0.166666666666667, 0.666666666666667, 0.333333333333333}
            }}
        )
    {
    }
};
template<>
class GaussLegendreTri<3> : public internal::GaussLegendreTriSuper<4>{
    public:
    GaussLegendreTri(const GaussLegendreTri& g) = delete;
    GaussLegendreTri(GaussLegendreTri&& g) = delete;
    void operator=(const GaussLegendreTri& g) = delete;
    void operator=(GaussLegendreTri&& g) = delete;

    static const GaussLegendreTri& get(){
        static GaussLegendreTri instance;

        return instance;
    }

    private:
    GaussLegendreTri()
    : internal::GaussLegendreTriSuper<4>(std::array<GLPointTri, 4>
            {{
                {0.333333333333333, 0.333333333333333, 0.333333333333333, -0.562500000000000},
                {0.600000000000000, 0.200000000000000, 0.200000000000000,  0.520833333333333},
                {0.200000000000000, 0.600000000000000, 0.200000000000000,  0.520833333333333},
                {0.200000000000000, 0.200000000000000, 0.600000000000000,  0.520833333333333}
            }}
        )
    {
    }
};
template<>
class GaussLegendreTri<4> : public internal::GaussLegendreTriSuper<6>{
    public:
    GaussLegendreTri(const GaussLegendreTri& g) = delete;
    GaussLegendreTri(GaussLegendreTri&& g) = delete;
    void operator=(const GaussLegendreTri& g) = delete;
    void operator=(GaussLegendreTri&& g) = delete;

    static const GaussLegendreTri& get(){
        static GaussLegendreTri instance;

        return instance;
    }

    private:
    GaussLegendreTri()
    : internal::GaussLegendreTriSuper<6>(std::array<GLPointTri, 6>
            {{
                {0.108103018168070, 0.445948490915965, 0.445948490915965, 0.223381589678011}, 
                {0.816847572980459, 0.091576213509771, 0.091576213509771, 0.109951743655322}, 
                {0.445948490915965, 0.108103018168070, 0.445948490915965, 0.223381589678011}, 
                {0.091576213509771, 0.816847572980459, 0.091576213509771, 0.109951743655322}, 
                {0.445948490915965, 0.445948490915965, 0.108103018168070, 0.223381589678011}, 
                {0.091576213509771, 0.091576213509771, 0.816847572980459, 0.109951743655322} 
            }}
        )
    {
    }
};
template<>
class GaussLegendreTri<5> : public internal::GaussLegendreTriSuper<7>{
    public:
    GaussLegendreTri(const GaussLegendreTri& g) = delete;
    GaussLegendreTri(GaussLegendreTri&& g) = delete;
    void operator=(const GaussLegendreTri& g) = delete;
    void operator=(GaussLegendreTri&& g) = delete;

    static const GaussLegendreTri& get(){
        static GaussLegendreTri instance;

        return instance;
    }
    
    private:
    GaussLegendreTri()
    : internal::GaussLegendreTriSuper<7>(std::array<GLPointTri, 7>
            {{
                {0.333333333333333, 0.333333333333333, 0.333333333333333, 0.225000000000000},
                {0.059715871789770, 0.470142064105115, 0.470142064105115, 0.132394152788506},
                {0.797426985353087, 0.101286507323456, 0.101286507323456, 0.125939180544827},
                {0.470142064105115, 0.059715871789770, 0.470142064105115, 0.132394152788506},
                {0.101286507323456, 0.797426985353087, 0.101286507323456, 0.125939180544827},
                {0.470142064105115, 0.470142064105115, 0.059715871789770, 0.132394152788506},
                {0.101286507323456, 0.101286507323456, 0.797426985353087, 0.125939180544827}
            }}
        )
    {
    }
};
template<>
class GaussLegendreTri<6> : public internal::GaussLegendreTriSuper<12>{
    public:
    GaussLegendreTri(const GaussLegendreTri& g) = delete;
    GaussLegendreTri(GaussLegendreTri&& g) = delete;
    void operator=(const GaussLegendreTri& g) = delete;
    void operator=(GaussLegendreTri&& g) = delete;

    static const GaussLegendreTri& get(){
        static GaussLegendreTri instance;

        return instance;
    }

    private:
    GaussLegendreTri()
    : internal::GaussLegendreTriSuper<12>(std::array<GLPointTri, 12>
            {{
                {0.501426509658179, 0.249286745170910, 0.249286745170910, 0.116786275726379},
                {0.873821971016996, 0.063089014491502, 0.063089014491502, 0.050844906370207},
                {0.249286745170910, 0.501426509658179, 0.249286745170910, 0.116786275726379},
                {0.063089014491502, 0.873821971016996, 0.063089014491502, 0.050844906370207},
                {0.249286745170910, 0.249286745170910, 0.501426509658179, 0.116786275726379},
                {0.063089014491502, 0.063089014491502, 0.873821971016996, 0.050844906370207},
                {0.053145049844817, 0.310352451033784, 0.636502499121399, 0.082851075618374},
                {0.310352451033784, 0.053145049844817, 0.636502499121399, 0.082851075618374},
                {0.310352451033784, 0.636502499121399, 0.053145049844817, 0.082851075618374},
                {0.053145049844817, 0.636502499121399, 0.310352451033784, 0.082851075618374},
                {0.636502499121399, 0.053145049844817, 0.310352451033784, 0.082851075618374},
                {0.636502499121399, 0.310352451033784, 0.053145049844817, 0.082851075618374}
            }}
        )
    {
    }
};

}

#endif
