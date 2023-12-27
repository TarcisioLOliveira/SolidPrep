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

/// TET ELEMENTS

struct GLPointTet{
    double a, b, c, d, w;
};

namespace internal{

template<size_t N>
class GaussLegendreTetSuper{
    public:
    GaussLegendreTetSuper(const GaussLegendreTetSuper& g) = delete;
    GaussLegendreTetSuper(GaussLegendreTetSuper&& g) = delete;
    void operator=(const GaussLegendreTetSuper& g) = delete;
    void operator=(GaussLegendreTetSuper&& g) = delete;

    auto begin() const{
        return p.cbegin();
    }

    auto end() const{
        return p.cend();
    }

    protected:
    GaussLegendreTetSuper(std::array<GLPointTet, N> p): 
        p(std::move(p))
    {
    }

    std::array<GLPointTet, N> p;
};

}

constexpr size_t P_to_N_Tet(size_t P){
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
class GaussLegendreTet : public internal::GaussLegendreTetSuper<P_to_N(P)>{
    public:
    GaussLegendreTet(const GaussLegendreTet& g) = delete;
    GaussLegendreTet(GaussLegendreTet&& g) = delete;
    void operator=(const GaussLegendreTet& g) = delete;
    void operator=(GaussLegendreTet&& g) = delete;

    static const GaussLegendreTet& get(){
        logger::log_assert(false, logger::ERROR, "Gauss Legendre for triangles not implemented for polynomials of order {}.", P);
        static GaussLegendreTet instance;

        return instance;
    }

    private:
    GaussLegendreTet()
    : internal::GaussLegendreTetSuper<P_to_N(P)>(std::array<GLPointTet, P_to_N(P)>())
    {
    }
};


// (Shunn & Ham, 2012)
template<>
class GaussLegendreTet<1> : public internal::GaussLegendreTetSuper<1>{
    public:
    GaussLegendreTet(const GaussLegendreTet& g) = delete;
    GaussLegendreTet(GaussLegendreTet&& g) = delete;
    void operator=(const GaussLegendreTet& g) = delete;
    void operator=(GaussLegendreTet&& g) = delete;

    static const GaussLegendreTet& get(){
        static GaussLegendreTet instance;

        return instance;
    }

    private:
    GaussLegendreTet()
    : internal::GaussLegendreTetSuper<1>(std::array<GLPointTet, 1>
            {{{1.0/4.0, 1.0/4.0, 1.0/4.0, 1.0/4.0, 1}}}
        )
    {
    }
};

template<>
class GaussLegendreTet<2> : public internal::GaussLegendreTetSuper<4>{
    public:
    GaussLegendreTet(const GaussLegendreTet& g) = delete;
    GaussLegendreTet(GaussLegendreTet&& g) = delete;
    void operator=(const GaussLegendreTet& g) = delete;
    void operator=(GaussLegendreTet&& g) = delete;

    static const GaussLegendreTet& get(){
        static GaussLegendreTet instance;

        return instance;
    }

    private:
    GaussLegendreTet()
    : internal::GaussLegendreTetSuper<4>(std::array<GLPointTet, 4>
            {{
                {0.5854101966249680, 0.1381966011250110, 0.1381966011250110, 0.1381966011250110, 0.2500000000000000},
                {0.1381966011250110, 0.5854101966249680, 0.1381966011250110, 0.1381966011250110, 0.2500000000000000},
                {0.1381966011250110, 0.1381966011250110, 0.5854101966249680, 0.1381966011250110, 0.2500000000000000},
                {0.1381966011250110, 0.1381966011250110, 0.1381966011250110, 0.5854101966249680, 0.2500000000000000}
            }}
        )
    {
    }
};

template<>
class GaussLegendreTet<3> : public internal::GaussLegendreTetSuper<10>{
    public:
    GaussLegendreTet(const GaussLegendreTet& g) = delete;
    GaussLegendreTet(GaussLegendreTet&& g) = delete;
    void operator=(const GaussLegendreTet& g) = delete;
    void operator=(GaussLegendreTet&& g) = delete;

    static const GaussLegendreTet& get(){
        static GaussLegendreTet instance;

        return instance;
    }

    private:
    GaussLegendreTet()
    : internal::GaussLegendreTetSuper<10>(std::array<GLPointTet, 10>
            {{
                {0.7784952948213300, 0.0738349017262234, 0.0738349017262234, 0.0738349017262234, 0.0476331348432089},
                {0.0738349017262234, 0.7784952948213300, 0.0738349017262234, 0.0738349017262234, 0.0476331348432089},
                {0.0738349017262234, 0.0738349017262234, 0.7784952948213300, 0.0738349017262234, 0.0476331348432089},
                {0.0738349017262234, 0.0738349017262234, 0.0738349017262234, 0.7784952948213300, 0.0476331348432089},
                {0.4062443438840510, 0.4062443438840510, 0.0937556561159491, 0.0937556561159491, 0.1349112434378610},
                {0.4062443438840510, 0.0937556561159491, 0.4062443438840510, 0.0937556561159491, 0.1349112434378610},
                {0.4062443438840510, 0.0937556561159491, 0.0937556561159491, 0.4062443438840510, 0.1349112434378610},
                {0.0937556561159491, 0.4062443438840510, 0.4062443438840510, 0.0937556561159491, 0.1349112434378610},
                {0.0937556561159491, 0.4062443438840510, 0.0937556561159491, 0.4062443438840510, 0.1349112434378610},
                {0.0937556561159491, 0.0937556561159491, 0.4062443438840510, 0.4062443438840510, 0.1349112434378610}
            }}
        )
    {
    }
};

template<>
class GaussLegendreTet<4> : public internal::GaussLegendreTetSuper<20>{
    public:
    GaussLegendreTet(const GaussLegendreTet& g) = delete;
    GaussLegendreTet(GaussLegendreTet&& g) = delete;
    void operator=(const GaussLegendreTet& g) = delete;
    void operator=(GaussLegendreTet&& g) = delete;

    static const GaussLegendreTet& get(){
        static GaussLegendreTet instance;

        return instance;
    }

    private:
    GaussLegendreTet()
    : internal::GaussLegendreTetSuper<20>(std::array<GLPointTet, 20>
            {{
                {0.9029422158182680, 0.0323525947272439, 0.0323525947272439, 0.0323525947272439, 0.0070670747944695},
                {0.0323525947272439, 0.9029422158182680, 0.0323525947272439, 0.0323525947272439, 0.0070670747944695},
                {0.0323525947272439, 0.0323525947272439, 0.9029422158182680, 0.0323525947272439, 0.0070670747944695},
                {0.0323525947272439, 0.0323525947272439, 0.0323525947272439, 0.9029422158182680, 0.0070670747944695},
                {0.2626825838877790, 0.6165965330619370, 0.0603604415251421, 0.0603604415251421, 0.0469986689718877},
                {0.6165965330619370, 0.2626825838877790, 0.0603604415251421, 0.0603604415251421, 0.0469986689718877},
                {0.2626825838877790, 0.0603604415251421, 0.6165965330619370, 0.0603604415251421, 0.0469986689718877},
                {0.6165965330619370, 0.0603604415251421, 0.2626825838877790, 0.0603604415251421, 0.0469986689718877},
                {0.2626825838877790, 0.0603604415251421, 0.0603604415251421, 0.6165965330619370, 0.0469986689718877},
                {0.6165965330619370, 0.0603604415251421, 0.0603604415251421, 0.2626825838877790, 0.0469986689718877},
                {0.0603604415251421, 0.2626825838877790, 0.6165965330619370, 0.0603604415251421, 0.0469986689718877},
                {0.0603604415251421, 0.6165965330619370, 0.2626825838877790, 0.0603604415251421, 0.0469986689718877},
                {0.0603604415251421, 0.2626825838877790, 0.0603604415251421, 0.6165965330619370, 0.0469986689718877},
                {0.0603604415251421, 0.6165965330619370, 0.0603604415251421, 0.2626825838877790, 0.0469986689718877},
                {0.0603604415251421, 0.0603604415251421, 0.2626825838877790, 0.6165965330619370, 0.0469986689718877},
                {0.0603604415251421, 0.0603604415251421, 0.6165965330619370, 0.2626825838877790, 0.0469986689718877},
                {0.3097693042728620, 0.3097693042728620, 0.3097693042728620, 0.0706920871814129, 0.1019369182898680},
                {0.3097693042728620, 0.3097693042728620, 0.0706920871814129, 0.3097693042728620, 0.1019369182898680},
                {0.3097693042728620, 0.0706920871814129, 0.3097693042728620, 0.3097693042728620, 0.1019369182898680},
                {0.0706920871814129, 0.3097693042728620, 0.3097693042728620, 0.3097693042728620, 0.1019369182898680}
            }}
        )
    {
    }
};

template<>
class GaussLegendreTet<5> : public internal::GaussLegendreTetSuper<35>{
    public:
    GaussLegendreTet(const GaussLegendreTet& g) = delete;
    GaussLegendreTet(GaussLegendreTet&& g) = delete;
    void operator=(const GaussLegendreTet& g) = delete;
    void operator=(GaussLegendreTet&& g) = delete;

    static const GaussLegendreTet& get(){
        static GaussLegendreTet instance;

        return instance;
    }

    private:
    GaussLegendreTet()
    : internal::GaussLegendreTetSuper<35>(std::array<GLPointTet, 35>
            {{
                {0.9197896733368800, 0.0267367755543735, 0.0267367755543735, 0.0267367755543735, 0.0021900463965388},
                {0.0267367755543735, 0.9197896733368800, 0.0267367755543735, 0.0267367755543735, 0.0021900463965388},
                {0.0267367755543735, 0.0267367755543735, 0.9197896733368800, 0.0267367755543735, 0.0021900463965388},
                {0.0267367755543735, 0.0267367755543735, 0.0267367755543735, 0.9197896733368800, 0.0021900463965388},
                {0.1740356302468940, 0.7477598884818090, 0.0391022406356488, 0.0391022406356488, 0.0143395670177665},
                {0.7477598884818090, 0.1740356302468940, 0.0391022406356488, 0.0391022406356488, 0.0143395670177665},
                {0.1740356302468940, 0.0391022406356488, 0.7477598884818090, 0.0391022406356488, 0.0143395670177665},
                {0.7477598884818090, 0.0391022406356488, 0.1740356302468940, 0.0391022406356488, 0.0143395670177665},
                {0.1740356302468940, 0.0391022406356488, 0.0391022406356488, 0.7477598884818090, 0.0143395670177665},
                {0.7477598884818090, 0.0391022406356488, 0.0391022406356488, 0.1740356302468940, 0.0143395670177665},
                {0.0391022406356488, 0.1740356302468940, 0.7477598884818090, 0.0391022406356488, 0.0143395670177665},
                {0.0391022406356488, 0.7477598884818090, 0.1740356302468940, 0.0391022406356488, 0.0143395670177665},
                {0.0391022406356488, 0.1740356302468940, 0.0391022406356488, 0.7477598884818090, 0.0143395670177665},
                {0.0391022406356488, 0.7477598884818090, 0.0391022406356488, 0.1740356302468940, 0.0143395670177665},
                {0.0391022406356488, 0.0391022406356488, 0.1740356302468940, 0.7477598884818090, 0.0143395670177665},
                {0.0391022406356488, 0.0391022406356488, 0.7477598884818090, 0.1740356302468940, 0.0143395670177665},
                {0.4547545999844830, 0.4547545999844830, 0.0452454000155172, 0.0452454000155172, 0.0250305395686746},
                {0.4547545999844830, 0.0452454000155172, 0.4547545999844830, 0.0452454000155172, 0.0250305395686746},
                {0.4547545999844830, 0.0452454000155172, 0.0452454000155172, 0.4547545999844830, 0.0250305395686746},
                {0.0452454000155172, 0.4547545999844830, 0.4547545999844830, 0.0452454000155172, 0.0250305395686746},
                {0.0452454000155172, 0.4547545999844830, 0.0452454000155172, 0.4547545999844830, 0.0250305395686746},
                {0.0452454000155172, 0.0452454000155172, 0.4547545999844830, 0.4547545999844830, 0.0250305395686746},
                {0.5031186450145980, 0.2232010379623150, 0.2232010379623150, 0.0504792790607720, 0.0479839333057554},
                {0.2232010379623150, 0.5031186450145980, 0.2232010379623150, 0.0504792790607720, 0.0479839333057554},
                {0.2232010379623150, 0.2232010379623150, 0.5031186450145980, 0.0504792790607720, 0.0479839333057554},
                {0.5031186450145980, 0.2232010379623150, 0.0504792790607720, 0.2232010379623150, 0.0479839333057554},
                {0.2232010379623150, 0.5031186450145980, 0.0504792790607720, 0.2232010379623150, 0.0479839333057554},
                {0.2232010379623150, 0.2232010379623150, 0.0504792790607720, 0.5031186450145980, 0.0479839333057554},
                {0.5031186450145980, 0.0504792790607720, 0.2232010379623150, 0.2232010379623150, 0.0479839333057554},
                {0.2232010379623150, 0.0504792790607720, 0.5031186450145980, 0.2232010379623150, 0.0479839333057554},
                {0.2232010379623150, 0.0504792790607720, 0.2232010379623150, 0.5031186450145980, 0.0479839333057554},
                {0.0504792790607720, 0.5031186450145980, 0.2232010379623150, 0.2232010379623150, 0.0479839333057554},
                {0.0504792790607720, 0.2232010379623150, 0.5031186450145980, 0.2232010379623150, 0.0479839333057554},
                {0.0504792790607720, 0.2232010379623150, 0.2232010379623150, 0.5031186450145980, 0.0479839333057554},
                {0.2500000000000000, 0.2500000000000000, 0.2500000000000000, 0.2500000000000000, 0.0931745731195340}
            }}
        )
    {
    }
};

}

#endif
