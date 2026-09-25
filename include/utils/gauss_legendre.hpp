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
        case 8: return 16;
        case 9: return 19;
        case 10: return 25;
        case 11: return 27;
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

template<>
class GaussLegendreTri<7> : public internal::GaussLegendreTriSuper<13>{
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
    : internal::GaussLegendreTriSuper<13>(std::array<GLPointTri, 13>
            {{
            {0.333333333333333, 0.333333333333333, 0.333333333333333, -0.149570044467682},
            {0.479308067841920, 0.260345966079040, 0.260345966079040, 0.175615257433208},
            {0.260345966079040, 0.260345966079040, 0.479308067841920, 0.175615257433208},
            {0.260345966079040, 0.479308067841920, 0.260345966079040, 0.175615257433208},
            {0.869739794195568, 0.065130102902216, 0.065130102902216, 0.053347235608838},
            {0.065130102902216, 0.065130102902216, 0.869739794195568, 0.053347235608838},
            {0.065130102902216, 0.869739794195568, 0.065130102902216, 0.053347235608838},
            {0.048690315425316, 0.312865496004874, 0.638444188569810, 0.077113760890257},
            {0.312865496004874, 0.638444188569810, 0.048690315425316, 0.077113760890257},
            {0.638444188569810, 0.048690315425316, 0.312865496004874, 0.077113760890257},
            {0.048690315425316, 0.638444188569810, 0.312865496004874, 0.077113760890257},
            {0.312865496004874, 0.048690315425316, 0.638444188569810, 0.077113760890257},
            {0.638444188569810, 0.312865496004874, 0.048690315425316, 0.077113760890257},
            }}
        )
    {
    }
};

template<>
class GaussLegendreTri<8> : public internal::GaussLegendreTriSuper<16>{
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
    : internal::GaussLegendreTriSuper<16>(std::array<GLPointTri, 16>
            {{
            {0.333333333333333, 0.333333333333333, 0.333333333333333, 0.144315607677787},
            {0.081414823414554, 0.459292588292723, 0.459292588292723, 0.095091634267285},
            {0.459292588292723, 0.459292588292723, 0.081414823414554, 0.095091634267285},
            {0.459292588292723, 0.081414823414554, 0.459292588292723, 0.095091634267285},
            {0.658861384496480, 0.170569307751760, 0.170569307751760, 0.103217370534718},
            {0.170569307751760, 0.170569307751760, 0.658861384496480, 0.103217370534718},
            {0.170569307751760, 0.658861384496480, 0.170569307751760, 0.103217370534718},
            {0.898905543365938, 0.050547228317031, 0.050547228317031, 0.032458497623198},
            {0.050547228317031, 0.050547228317031, 0.898905543365938, 0.032458497623198},
            {0.050547228317031, 0.898905543365938, 0.050547228317031, 0.032458497623198},
            {0.008394777409958, 0.263112829634638, 0.728492392955404, 0.027230314174435},
            {0.263112829634638, 0.728492392955404, 0.008394777409958, 0.027230314174435},
            {0.728492392955404, 0.008394777409958, 0.263112829634638, 0.027230314174435},
            {0.008394777409958, 0.728492392955404, 0.263112829634638, 0.027230314174435},
            {0.263112829634638, 0.008394777409958, 0.728492392955404, 0.027230314174435},
            {0.728492392955404, 0.263112829634638, 0.008394777409958, 0.027230314174435},
            }}
        )
    {
    }
};

template<>
class GaussLegendreTri<9> : public internal::GaussLegendreTriSuper<19>{
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
    : internal::GaussLegendreTriSuper<19>(std::array<GLPointTri, 19>
            {{
            {0.333333333333333, 0.333333333333333, 0.333333333333333, 0.097135796282799},
            {0.020634961602525, 0.489682519198738, 0.489682519198738, 0.031334700227139},
            {0.489682519198738, 0.489682519198738, 0.020634961602525, 0.031334700227139},
            {0.489682519198738, 0.020634961602525, 0.489682519198738, 0.031334700227139},
            {0.125820817014127, 0.437089591492937, 0.437089591492937, 0.077827541004774},
            {0.437089591492937, 0.437089591492937, 0.125820817014127, 0.077827541004774},
            {0.437089591492937, 0.125820817014127, 0.437089591492937, 0.077827541004774},
            {0.623592928761935, 0.188203535619033, 0.188203535619033, 0.079647738927210},
            {0.188203535619033, 0.188203535619033, 0.623592928761935, 0.079647738927210},
            {0.188203535619033, 0.623592928761935, 0.188203535619033, 0.079647738927210},
            {0.910540973211095, 0.044729513394453, 0.044729513394453, 0.025577675658698},
            {0.044729513394453, 0.044729513394453, 0.910540973211095, 0.025577675658698},
            {0.044729513394453, 0.910540973211095, 0.044729513394453, 0.025577675658698},
            {0.036838412054736, 0.221962989160766, 0.741198598784498, 0.043283539377289},
            {0.221962989160766, 0.741198598784498, 0.036838412054736, 0.043283539377289},
            {0.741198598784498, 0.036838412054736, 0.221962989160766, 0.043283539377289},
            {0.036838412054736, 0.741198598784498, 0.221962989160766, 0.043283539377289},
            {0.221962989160766, 0.036838412054736, 0.741198598784498, 0.043283539377289},
            {0.741198598784498, 0.221962989160766, 0.036838412054736, 0.043283539377289},
            }}
        )
    {
    }
};

template<>
class GaussLegendreTri<10> : public internal::GaussLegendreTriSuper<25>{
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
    : internal::GaussLegendreTriSuper<25>(std::array<GLPointTri, 25>
            {{
            {0.333333333333333, 0.333333333333333, 0.333333333333333, 0.090817990382754},
            {0.028844733232685, 0.485577633383657, 0.485577633383657, 0.036725957756467},
            {0.485577633383657, 0.485577633383657, 0.028844733232685, 0.036725957756467},
            {0.485577633383657, 0.028844733232685, 0.485577633383657, 0.036725957756467},
            {0.781036849029926, 0.109481575485037, 0.109481575485037, 0.045321059435528},
            {0.109481575485037, 0.109481575485037, 0.781036849029926, 0.045321059435528},
            {0.109481575485037, 0.781036849029926, 0.109481575485037, 0.045321059435528},
            {0.141707219414880, 0.307939838764121, 0.550352941820999, 0.072757916845420},
            {0.307939838764121, 0.550352941820999, 0.141707219414880, 0.072757916845420},
            {0.550352941820999, 0.141707219414880, 0.307939838764121, 0.072757916845420},
            {0.141707219414880, 0.550352941820999, 0.307939838764121, 0.072757916845420},
            {0.307939838764121, 0.141707219414880, 0.550352941820999, 0.072757916845420},
            {0.550352941820999, 0.307939838764121, 0.141707219414880, 0.072757916845420},
            {0.025003534762686, 0.246672560639903, 0.728323904597411, 0.028327242531057},
            {0.246672560639903, 0.728323904597411, 0.025003534762686, 0.028327242531057},
            {0.728323904597411, 0.025003534762686, 0.246672560639903, 0.028327242531057},
            {0.025003534762686, 0.728323904597411, 0.246672560639903, 0.028327242531057},
            {0.246672560639903, 0.025003534762686, 0.728323904597411, 0.028327242531057},
            {0.728323904597411, 0.246672560639903, 0.025003534762686, 0.028327242531057},
            {0.009540815400299, 0.066803251012200, 0.923655933587500, 0.009421666963733},
            {0.066803251012200, 0.923655933587500, 0.009540815400299, 0.009421666963733},
            {0.923655933587500, 0.009540815400299, 0.066803251012200, 0.009421666963733},
            {0.009540815400299, 0.923655933587500, 0.066803251012200, 0.009421666963733},
            {0.066803251012200, 0.009540815400299, 0.923655933587500, 0.009421666963733},
            {0.923655933587500, 0.066803251012200, 0.009540815400299, 0.009421666963733},
            }}
        )
    {
    }
};

template<>
class GaussLegendreTri<11> : public internal::GaussLegendreTriSuper<27>{
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
    : internal::GaussLegendreTriSuper<27>(std::array<GLPointTri, 27>
            {{
            {-0.069222096541517, 0.534611048270758, 0.534611048270758, 0.000927006328961},
            {0.534611048270758, 0.534611048270758, -0.069222096541517, 0.000927006328961},
            {0.534611048270758, -0.069222096541517, 0.534611048270758, 0.000927006328961},
            {0.202061394068290, 0.398969302965855, 0.398969302965855, 0.077149534914813},
            {0.398969302965855, 0.398969302965855, 0.202061394068290, 0.077149534914813},
            {0.398969302965855, 0.202061394068290, 0.398969302965855, 0.077149534914813},
            {0.593380199137435, 0.203309900431282, 0.203309900431282, 0.059322977380774},
            {0.203309900431282, 0.203309900431282, 0.593380199137435, 0.059322977380774},
            {0.203309900431282, 0.593380199137435, 0.203309900431282, 0.059322977380774},
            {0.761298175434837, 0.119350912282581, 0.119350912282581, 0.036184540503418},
            {0.119350912282581, 0.119350912282581, 0.761298175434837, 0.036184540503418},
            {0.119350912282581, 0.761298175434837, 0.119350912282581, 0.036184540503418},
            {0.935270103777448, 0.032364948111276, 0.032364948111276, 0.013659731002678},
            {0.032364948111276, 0.032364948111276, 0.935270103777448, 0.013659731002678},
            {0.032364948111276, 0.935270103777448, 0.032364948111276, 0.013659731002678},
            {0.050178138310495, 0.356620648261293, 0.593201213428213, 0.052337111962204},
            {0.356620648261293, 0.593201213428213, 0.050178138310495, 0.052337111962204},
            {0.593201213428213, 0.050178138310495, 0.356620648261293, 0.052337111962204},
            {0.050178138310495, 0.593201213428213, 0.356620648261293, 0.052337111962204},
            {0.356620648261293, 0.050178138310495, 0.593201213428213, 0.052337111962204},
            {0.593201213428213, 0.356620648261293, 0.050178138310495, 0.052337111962204},
            {0.021022016536166, 0.171488980304042, 0.807489003159792, 0.020707659639141},
            {0.171488980304042, 0.807489003159792, 0.021022016536166, 0.020707659639141},
            {0.807489003159792, 0.021022016536166, 0.171488980304042, 0.020707659639141},
            {0.021022016536166, 0.807489003159792, 0.171488980304042, 0.020707659639141},
            {0.171488980304042, 0.021022016536166, 0.807489003159792, 0.020707659639141},
            {0.807489003159792, 0.171488980304042, 0.021022016536166, 0.020707659639141},
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
        case 2: return 4;
        case 3: return 8;
        case 4: return 14;
        case 5: return 14;
        case 6: return 24;
        case 7: return 35;
        case 8: return 46;
        case 9: return 59;
        case 10: return 81;
        case 11: return 110;
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


// (Jaśkowiec & Sukumar, 2020)
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
            {0.138196601125010515, 0.138196601125010515, 0.138196601125010515, 0.585410196624968454, 0.25},
            {0.585410196624968454, 0.138196601125010515, 0.138196601125010515, 0.138196601125010515, 0.25},
            {0.138196601125010515, 0.585410196624968454, 0.138196601125010515, 0.138196601125010515, 0.25},
            {0.138196601125010515, 0.138196601125010515, 0.585410196624968454, 0.138196601125010515, 0.25},
            }}
        )
    {
    }
};

template<>
class GaussLegendreTet<3> : public internal::GaussLegendreTetSuper<8>{
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
    : internal::GaussLegendreTetSuper<8>(std::array<GLPointTet, 8>
            {{
            {0.328681146665348992, 0.328681146665348992, 0.328681146665348992, 0.013956560003953021, 0.127491311557506424},
            {0.013956560003953021, 0.328681146665348992, 0.328681146665348992, 0.328681146665348992, 0.127491311557506424},
            {0.328681146665348992, 0.013956560003953021, 0.328681146665348992, 0.328681146665348992, 0.127491311557506424},
            {0.328681146665348992, 0.328681146665348992, 0.013956560003953021, 0.328681146665348992, 0.127491311557506424},
            {0.111920727509291548, 0.111920727509291548, 0.111920727509291548, 0.664237817472125354, 0.122508688442493575},
            {0.664237817472125354, 0.111920727509291548, 0.111920727509291548, 0.111920727509291548, 0.122508688442493575},
            {0.111920727509291548, 0.664237817472125354, 0.111920727509291548, 0.111920727509291548, 0.122508688442493575},
            {0.111920727509291548, 0.111920727509291548, 0.664237817472125354, 0.111920727509291548, 0.122508688442493575},
            }}
        )
    {
    }
};

template<>
class GaussLegendreTet<4> : public internal::GaussLegendreTetSuper<14>{
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
    : internal::GaussLegendreTetSuper<14>(std::array<GLPointTet, 14>
            {{
            {0.309582165317951906, 0.309582165317951906, 0.309582165317951906, 0.071253504046144279, 0.065409163634452181},
            {0.071253504046144279, 0.309582165317951906, 0.309582165317951906, 0.309582165317951906, 0.065409163634452181},
            {0.309582165317951906, 0.071253504046144279, 0.309582165317951906, 0.309582165317951906, 0.065409163634452181},
            {0.309582165317951906, 0.309582165317951906, 0.071253504046144279, 0.309582165317951906, 0.065409163634452181},
            {0.083442772303977572, 0.083442772303977572, 0.083442772303977572, 0.749671683088067283, 0.059355264231004755},
            {0.749671683088067283, 0.083442772303977572, 0.083442772303977572, 0.083442772303977572, 0.059355264231004755},
            {0.083442772303977572, 0.749671683088067283, 0.083442772303977572, 0.083442772303977572, 0.059355264231004755},
            {0.083442772303977572, 0.083442772303977572, 0.749671683088067283, 0.083442772303977572, 0.059355264231004755},
            {0.077220954070058665, 0.422779045929941334, 0.422779045929941334, 0.077220954070058665, 0.083490381423028708},
            {0.422779045929941334, 0.077220954070058665, 0.422779045929941334, 0.077220954070058665, 0.083490381423028708},
            {0.077220954070058665, 0.077220954070058665, 0.422779045929941334, 0.422779045929941334, 0.083490381423028708},
            {0.422779045929941334, 0.422779045929941334, 0.077220954070058665, 0.077220954070058665, 0.083490381423028708},
            {0.077220954070058665, 0.422779045929941334, 0.077220954070058665, 0.422779045929941334, 0.083490381423028708},
            {0.422779045929941334, 0.077220954070058665, 0.077220954070058665, 0.422779045929941334, 0.083490381423028708},
            }}
        )
    {
    }
};

template<>
class GaussLegendreTet<5> : public internal::GaussLegendreTetSuper<14>{
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
    : internal::GaussLegendreTetSuper<14>(std::array<GLPointTet, 14>
            {{
            {0.310885919263300609, 0.310885919263300609, 0.310885919263300609, 0.067342242210098170, 0.112687925718015850},
            {0.067342242210098170, 0.310885919263300609, 0.310885919263300609, 0.310885919263300609, 0.112687925718015850},
            {0.310885919263300609, 0.067342242210098170, 0.310885919263300609, 0.310885919263300609, 0.112687925718015850},
            {0.310885919263300609, 0.310885919263300609, 0.067342242210098170, 0.310885919263300609, 0.112687925718015850},
            {0.092735250310891226, 0.092735250310891226, 0.092735250310891226, 0.721794249067326320, 0.073493043116361949},
            {0.721794249067326320, 0.092735250310891226, 0.092735250310891226, 0.092735250310891226, 0.073493043116361949},
            {0.092735250310891226, 0.721794249067326320, 0.092735250310891226, 0.092735250310891226, 0.073493043116361949},
            {0.092735250310891226, 0.092735250310891226, 0.721794249067326320, 0.092735250310891226, 0.073493043116361949},
            {0.045503704125649649, 0.454496295874350350, 0.454496295874350350, 0.045503704125649649, 0.042546020777081466},
            {0.454496295874350350, 0.045503704125649649, 0.454496295874350350, 0.045503704125649649, 0.042546020777081466},
            {0.045503704125649649, 0.045503704125649649, 0.454496295874350350, 0.454496295874350350, 0.042546020777081466},
            {0.454496295874350350, 0.454496295874350350, 0.045503704125649649, 0.045503704125649649, 0.042546020777081466},
            {0.045503704125649649, 0.454496295874350350, 0.045503704125649649, 0.454496295874350350, 0.042546020777081466},
            {0.454496295874350350, 0.045503704125649649, 0.045503704125649649, 0.454496295874350350, 0.042546020777081466},
            }}
        )
    {
    }
};

template<>
class GaussLegendreTet<6> : public internal::GaussLegendreTetSuper<24>{
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
    : internal::GaussLegendreTetSuper<24>(std::array<GLPointTet, 24>
            {{
            {0.040673958534611353, 0.040673958534611353, 0.040673958534611353, 0.877978124396165940, 0.010077211055320642},
            {0.877978124396165940, 0.040673958534611353, 0.040673958534611353, 0.040673958534611353, 0.010077211055320642},
            {0.040673958534611353, 0.877978124396165940, 0.040673958534611353, 0.040673958534611353, 0.010077211055320642},
            {0.040673958534611353, 0.040673958534611353, 0.877978124396165940, 0.040673958534611353, 0.010077211055320642},
            {0.322337890142275510, 0.322337890142275510, 0.322337890142275510, 0.032986329573173468, 0.055357181543654722},
            {0.032986329573173468, 0.322337890142275510, 0.322337890142275510, 0.322337890142275510, 0.055357181543654722},
            {0.322337890142275510, 0.032986329573173468, 0.322337890142275510, 0.322337890142275510, 0.055357181543654722},
            {0.322337890142275510, 0.322337890142275510, 0.032986329573173468, 0.322337890142275510, 0.055357181543654722},
            {0.214602871259152029, 0.214602871259152029, 0.214602871259152029, 0.356191386222543912, 0.039922750258167492},
            {0.356191386222543912, 0.214602871259152029, 0.214602871259152029, 0.214602871259152029, 0.039922750258167492},
            {0.214602871259152029, 0.356191386222543912, 0.214602871259152029, 0.214602871259152029, 0.039922750258167492},
            {0.214602871259152029, 0.214602871259152029, 0.356191386222543912, 0.214602871259152029, 0.039922750258167492},
            {0.603005664791649141, 0.063661001875017525, 0.063661001875017525, 0.269672331458315808, 0.048214285714285714},
            {0.063661001875017525, 0.603005664791649141, 0.063661001875017525, 0.269672331458315808, 0.048214285714285714},
            {0.063661001875017525, 0.063661001875017525, 0.603005664791649141, 0.269672331458315808, 0.048214285714285714},
            {0.269672331458315808, 0.063661001875017525, 0.063661001875017525, 0.603005664791649141, 0.048214285714285714},
            {0.269672331458315808, 0.603005664791649141, 0.063661001875017525, 0.063661001875017525, 0.048214285714285714},
            {0.269672331458315808, 0.063661001875017525, 0.603005664791649141, 0.063661001875017525, 0.048214285714285714},
            {0.063661001875017525, 0.269672331458315808, 0.063661001875017525, 0.603005664791649141, 0.048214285714285714},
            {0.603005664791649141, 0.269672331458315808, 0.063661001875017525, 0.063661001875017525, 0.048214285714285714},
            {0.063661001875017525, 0.269672331458315808, 0.603005664791649141, 0.063661001875017525, 0.048214285714285714},
            {0.063661001875017525, 0.063661001875017525, 0.269672331458315808, 0.603005664791649141, 0.048214285714285714},
            {0.603005664791649141, 0.063661001875017525, 0.269672331458315808, 0.063661001875017525, 0.048214285714285714},
            {0.063661001875017525, 0.603005664791649141, 0.269672331458315808, 0.063661001875017525, 0.048214285714285714},
            }}
        )
    {
    }
};

template<>
class GaussLegendreTet<7> : public internal::GaussLegendreTetSuper<35>{
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
            {0.25, 0.25, 0.25, 0.25, 0.095485289464130848},
            {0.315701149778202799, 0.315701149778202799, 0.315701149778202799, 0.052896550665391601, 0.042329581209967029},
            {0.052896550665391601, 0.315701149778202799, 0.315701149778202799, 0.315701149778202799, 0.042329581209967029},
            {0.315701149778202799, 0.052896550665391601, 0.315701149778202799, 0.315701149778202799, 0.042329581209967029},
            {0.315701149778202799, 0.315701149778202799, 0.052896550665391601, 0.315701149778202799, 0.042329581209967029},
            {0.050489822598396368, 0.449510177401603631, 0.449510177401603631, 0.050489822598396368, 0.031896927832857579},
            {0.449510177401603631, 0.050489822598396368, 0.449510177401603631, 0.050489822598396368, 0.031896927832857579},
            {0.050489822598396368, 0.050489822598396368, 0.449510177401603631, 0.449510177401603631, 0.031896927832857579},
            {0.449510177401603631, 0.449510177401603631, 0.050489822598396368, 0.050489822598396368, 0.031896927832857579},
            {0.050489822598396368, 0.449510177401603631, 0.050489822598396368, 0.449510177401603631, 0.031896927832857579},
            {0.449510177401603631, 0.050489822598396368, 0.050489822598396368, 0.449510177401603631, 0.031896927832857579},
            {0.810830241098548561, 0.021265472541483245, 0.021265472541483245, 0.146638813818484946, 0.008110770829903341},
            {0.021265472541483245, 0.810830241098548561, 0.021265472541483245, 0.146638813818484946, 0.008110770829903341},
            {0.021265472541483245, 0.021265472541483245, 0.810830241098548561, 0.146638813818484946, 0.008110770829903341},
            {0.146638813818484946, 0.021265472541483245, 0.021265472541483245, 0.810830241098548561, 0.008110770829903341},
            {0.146638813818484946, 0.810830241098548561, 0.021265472541483245, 0.021265472541483245, 0.008110770829903341},
            {0.146638813818484946, 0.021265472541483245, 0.810830241098548561, 0.021265472541483245, 0.008110770829903341},
            {0.021265472541483245, 0.146638813818484946, 0.021265472541483245, 0.810830241098548561, 0.008110770829903341},
            {0.810830241098548561, 0.146638813818484946, 0.021265472541483245, 0.021265472541483245, 0.008110770829903341},
            {0.021265472541483245, 0.146638813818484946, 0.810830241098548561, 0.021265472541483245, 0.008110770829903341},
            {0.021265472541483245, 0.021265472541483245, 0.146638813818484946, 0.810830241098548561, 0.008110770829903341},
            {0.810830241098548561, 0.021265472541483245, 0.146638813818484946, 0.021265472541483245, 0.008110770829903341},
            {0.021265472541483245, 0.810830241098548561, 0.146638813818484946, 0.021265472541483245, 0.008110770829903341},
            {0.575171637587000023, 0.188833831026001047, 0.188833831026001047, 0.047160700360997881, 0.037207130728334621},
            {0.188833831026001047, 0.575171637587000023, 0.188833831026001047, 0.047160700360997881, 0.037207130728334621},
            {0.188833831026001047, 0.188833831026001047, 0.575171637587000023, 0.047160700360997881, 0.037207130728334621},
            {0.047160700360997881, 0.188833831026001047, 0.188833831026001047, 0.575171637587000023, 0.037207130728334621},
            {0.047160700360997881, 0.575171637587000023, 0.188833831026001047, 0.188833831026001047, 0.037207130728334621},
            {0.047160700360997881, 0.188833831026001047, 0.575171637587000023, 0.188833831026001047, 0.037207130728334621},
            {0.188833831026001047, 0.047160700360997881, 0.188833831026001047, 0.575171637587000023, 0.037207130728334621},
            {0.575171637587000023, 0.047160700360997881, 0.188833831026001047, 0.188833831026001047, 0.037207130728334621},
            {0.188833831026001047, 0.047160700360997881, 0.575171637587000023, 0.188833831026001047, 0.037207130728334621},
            {0.188833831026001047, 0.188833831026001047, 0.047160700360997881, 0.575171637587000023, 0.037207130728334621},
            {0.575171637587000023, 0.188833831026001047, 0.047160700360997881, 0.188833831026001047, 0.037207130728334621},
            {0.188833831026001047, 0.575171637587000023, 0.047160700360997881, 0.188833831026001047, 0.037207130728334621},
            }}
        )
    {
    }
};

template<>
class GaussLegendreTet<8> : public internal::GaussLegendreTetSuper<46>{
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
    : internal::GaussLegendreTetSuper<46>(std::array<GLPointTet, 46>
            {{
            {0.186812758270720971, 0.186812758270720971, 0.186812758270720971, 0.439561725187837086, 0.047812988034306654},
            {0.439561725187837086, 0.186812758270720971, 0.186812758270720971, 0.186812758270720971, 0.047812988034306654},
            {0.186812758270720971, 0.439561725187837086, 0.186812758270720971, 0.186812758270720971, 0.047812988034306654},
            {0.186812758270720971, 0.186812758270720971, 0.439561725187837086, 0.186812758270720971, 0.047812988034306654},
            {0.114462406761230491, 0.114462406761230491, 0.114462406761230491, 0.656612779716308524, 0.029727854086794222},
            {0.656612779716308524, 0.114462406761230491, 0.114462406761230491, 0.114462406761230491, 0.029727854086794222},
            {0.114462406761230491, 0.656612779716308524, 0.114462406761230491, 0.114462406761230491, 0.029727854086794222},
            {0.114462406761230491, 0.114462406761230491, 0.656612779716308524, 0.114462406761230491, 0.029727854086794222},
            {0.313806592272401620, 0.313806592272401620, 0.313806592272401620, 0.058580223182795139, 0.043527324218975264},
            {0.058580223182795139, 0.313806592272401620, 0.313806592272401620, 0.313806592272401620, 0.043527324218975264},
            {0.313806592272401620, 0.058580223182795139, 0.313806592272401620, 0.313806592272401620, 0.043527324218975264},
            {0.313806592272401620, 0.313806592272401620, 0.058580223182795139, 0.313806592272401620, 0.043527324218975264},
            {0.044577379448846248, 0.044577379448846248, 0.044577379448846248, 0.866267861653461255, 0.008632736020814974},
            {0.866267861653461255, 0.044577379448846248, 0.044577379448846248, 0.044577379448846248, 0.008632736020814974},
            {0.044577379448846248, 0.866267861653461255, 0.044577379448846248, 0.044577379448846248, 0.008632736020814974},
            {0.044577379448846248, 0.044577379448846248, 0.866267861653461255, 0.044577379448846248, 0.008632736020814974},
            {0.434516509456157605, 0.065483490543842394, 0.065483490543842394, 0.434516509456157605, 0.036890541269866362},
            {0.065483490543842394, 0.434516509456157605, 0.065483490543842394, 0.434516509456157605, 0.036890541269866362},
            {0.434516509456157605, 0.434516509456157605, 0.065483490543842394, 0.065483490543842394, 0.036890541269866362},
            {0.065483490543842394, 0.065483490543842394, 0.434516509456157605, 0.434516509456157605, 0.036890541269866362},
            {0.434516509456157605, 0.065483490543842394, 0.434516509456157605, 0.065483490543842394, 0.036890541269866362},
            {0.065483490543842394, 0.434516509456157605, 0.434516509456157605, 0.065483490543842394, 0.036890541269866362},
            {0.005073320750421590, 0.203893174662110274, 0.203893174662110274, 0.587140329925357860, 0.014497026710859505},
            {0.203893174662110274, 0.005073320750421590, 0.203893174662110274, 0.587140329925357860, 0.014497026710859505},
            {0.203893174662110274, 0.203893174662110274, 0.005073320750421590, 0.587140329925357860, 0.014497026710859505},
            {0.587140329925357860, 0.203893174662110274, 0.203893174662110274, 0.005073320750421590, 0.014497026710859505},
            {0.587140329925357860, 0.005073320750421590, 0.203893174662110274, 0.203893174662110274, 0.014497026710859505},
            {0.587140329925357860, 0.203893174662110274, 0.005073320750421590, 0.203893174662110274, 0.014497026710859505},
            {0.203893174662110274, 0.587140329925357860, 0.203893174662110274, 0.005073320750421590, 0.014497026710859505},
            {0.005073320750421590, 0.587140329925357860, 0.203893174662110274, 0.203893174662110274, 0.014497026710859505},
            {0.203893174662110274, 0.587140329925357860, 0.005073320750421590, 0.203893174662110274, 0.014497026710859505},
            {0.203893174662110274, 0.203893174662110274, 0.587140329925357860, 0.005073320750421590, 0.014497026710859505},
            {0.005073320750421590, 0.203893174662110274, 0.587140329925357860, 0.203893174662110274, 0.014497026710859505},
            {0.203893174662110274, 0.005073320750421590, 0.587140329925357860, 0.203893174662110274, 0.014497026710859505},
            {0.714676926393304887, 0.021247779669571673, 0.021247779669571673, 0.242827514267551765, 0.007157401867243607},
            {0.021247779669571673, 0.714676926393304887, 0.021247779669571673, 0.242827514267551765, 0.007157401867243607},
            {0.021247779669571673, 0.021247779669571673, 0.714676926393304887, 0.242827514267551765, 0.007157401867243607},
            {0.242827514267551765, 0.021247779669571673, 0.021247779669571673, 0.714676926393304887, 0.007157401867243607},
            {0.242827514267551765, 0.714676926393304887, 0.021247779669571673, 0.021247779669571673, 0.007157401867243607},
            {0.242827514267551765, 0.021247779669571673, 0.714676926393304887, 0.021247779669571673, 0.007157401867243607},
            {0.021247779669571673, 0.242827514267551765, 0.021247779669571673, 0.714676926393304887, 0.007157401867243607},
            {0.714676926393304887, 0.242827514267551765, 0.021247779669571673, 0.021247779669571673, 0.007157401867243607},
            {0.021247779669571673, 0.242827514267551765, 0.714676926393304887, 0.021247779669571673, 0.007157401867243607},
            {0.021247779669571673, 0.021247779669571673, 0.242827514267551765, 0.714676926393304887, 0.007157401867243607},
            {0.714676926393304887, 0.021247779669571673, 0.242827514267551765, 0.021247779669571673, 0.007157401867243607},
            {0.021247779669571673, 0.714676926393304887, 0.242827514267551765, 0.021247779669571673, 0.007157401867243607},
            }}
        )
    {
    }
};

template<>
class GaussLegendreTet<9> : public internal::GaussLegendreTetSuper<59>{
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
    : internal::GaussLegendreTetSuper<59>(std::array<GLPointTet, 59>
            {{
            {0.25, 0.25, 0.25, 0.25, 0.056866624253551728},
            {0.167906605236742789, 0.167906605236742789, 0.167906605236742789, 0.496280184289771631, 0.025694276689523175},
            {0.496280184289771631, 0.167906605236742789, 0.167906605236742789, 0.167906605236742789, 0.025694276689523175},
            {0.167906605236742789, 0.496280184289771631, 0.167906605236742789, 0.167906605236742789, 0.025694276689523175},
            {0.167906605236742789, 0.167906605236742789, 0.496280184289771631, 0.167906605236742789, 0.025694276689523175},
            {0.091436270514076258, 0.091436270514076258, 0.091436270514076258, 0.725691188457771224, 0.002298923353028369},
            {0.725691188457771224, 0.091436270514076258, 0.091436270514076258, 0.091436270514076258, 0.002298923353028369},
            {0.091436270514076258, 0.725691188457771224, 0.091436270514076258, 0.091436270514076258, 0.002298923353028369},
            {0.091436270514076258, 0.091436270514076258, 0.725691188457771224, 0.091436270514076258, 0.002298923353028369},
            {0.321855664817653266, 0.321855664817653266, 0.321855664817653266, 0.034433005547040201, 0.030456038667591007},
            {0.034433005547040201, 0.321855664817653266, 0.321855664817653266, 0.321855664817653266, 0.030456038667591007},
            {0.321855664817653266, 0.034433005547040201, 0.321855664817653266, 0.321855664817653266, 0.030456038667591007},
            {0.321855664817653266, 0.321855664817653266, 0.034433005547040201, 0.321855664817653266, 0.030456038667591007},
            {0.041837695900365602, 0.041837695900365602, 0.041837695900365602, 0.874486912298903193, 0.007123722340238880},
            {0.874486912298903193, 0.041837695900365602, 0.041837695900365602, 0.041837695900365602, 0.007123722340238880},
            {0.041837695900365602, 0.874486912298903193, 0.041837695900365602, 0.041837695900365602, 0.007123722340238880},
            {0.041837695900365602, 0.041837695900365602, 0.874486912298903193, 0.041837695900365602, 0.007123722340238880},
            {0.393270561425153639, 0.106729438574846360, 0.106729438574846360, 0.393270561425153639, 0.036843655480943882},
            {0.106729438574846360, 0.393270561425153639, 0.106729438574846360, 0.393270561425153639, 0.036843655480943882},
            {0.393270561425153639, 0.393270561425153639, 0.106729438574846360, 0.106729438574846360, 0.036843655480943882},
            {0.106729438574846360, 0.106729438574846360, 0.393270561425153639, 0.393270561425153639, 0.036843655480943882},
            {0.393270561425153639, 0.106729438574846360, 0.393270561425153639, 0.106729438574846360, 0.036843655480943882},
            {0.106729438574846360, 0.393270561425153639, 0.393270561425153639, 0.106729438574846360, 0.036843655480943882},
            {0.718393093981424399, 0.033957168183088885, 0.033957168183088885, 0.213692569652397829, 0.010322056782209846},
            {0.033957168183088885, 0.718393093981424399, 0.033957168183088885, 0.213692569652397829, 0.010322056782209846},
            {0.033957168183088885, 0.033957168183088885, 0.718393093981424399, 0.213692569652397829, 0.010322056782209846},
            {0.213692569652397829, 0.033957168183088885, 0.033957168183088885, 0.718393093981424399, 0.010322056782209846},
            {0.213692569652397829, 0.718393093981424399, 0.033957168183088885, 0.033957168183088885, 0.010322056782209846},
            {0.213692569652397829, 0.033957168183088885, 0.718393093981424399, 0.033957168183088885, 0.010322056782209846},
            {0.033957168183088885, 0.213692569652397829, 0.033957168183088885, 0.718393093981424399, 0.010322056782209846},
            {0.718393093981424399, 0.213692569652397829, 0.033957168183088885, 0.033957168183088885, 0.010322056782209846},
            {0.033957168183088885, 0.213692569652397829, 0.718393093981424399, 0.033957168183088885, 0.010322056782209846},
            {0.033957168183088885, 0.033957168183088885, 0.213692569652397829, 0.718393093981424399, 0.010322056782209846},
            {0.718393093981424399, 0.033957168183088885, 0.213692569652397829, 0.033957168183088885, 0.010322056782209846},
            {0.033957168183088885, 0.718393093981424399, 0.213692569652397829, 0.033957168183088885, 0.010322056782209846},
            {0.078683634476687930, 0.460658181054777632, 0.460658181054777632, 0.000000003413756803, 0.007634887153980854},
            {0.460658181054777632, 0.078683634476687930, 0.460658181054777632, 0.000000003413756803, 0.007634887153980854},
            {0.460658181054777632, 0.460658181054777632, 0.078683634476687930, 0.000000003413756803, 0.007634887153980854},
            {0.000000003413756803, 0.460658181054777632, 0.460658181054777632, 0.078683634476687930, 0.007634887153980854},
            {0.000000003413756803, 0.078683634476687930, 0.460658181054777632, 0.460658181054777632, 0.007634887153980854},
            {0.000000003413756803, 0.460658181054777632, 0.078683634476687930, 0.460658181054777632, 0.007634887153980854},
            {0.460658181054777632, 0.000000003413756803, 0.460658181054777632, 0.078683634476687930, 0.007634887153980854},
            {0.078683634476687930, 0.000000003413756803, 0.460658181054777632, 0.460658181054777632, 0.007634887153980854},
            {0.460658181054777632, 0.000000003413756803, 0.078683634476687930, 0.460658181054777632, 0.007634887153980854},
            {0.460658181054777632, 0.460658181054777632, 0.000000003413756803, 0.078683634476687930, 0.007634887153980854},
            {0.078683634476687930, 0.460658181054777632, 0.000000003413756803, 0.460658181054777632, 0.007634887153980854},
            {0.460658181054777632, 0.078683634476687930, 0.000000003413756803, 0.460658181054777632, 0.007634887153980854},
            {0.597210722761849931, 0.184387943500015203, 0.184387943500015203, 0.034013390238119661, 0.020358022618747569},
            {0.184387943500015203, 0.597210722761849931, 0.184387943500015203, 0.034013390238119661, 0.020358022618747569},
            {0.184387943500015203, 0.184387943500015203, 0.597210722761849931, 0.034013390238119661, 0.020358022618747569},
            {0.034013390238119661, 0.184387943500015203, 0.184387943500015203, 0.597210722761849931, 0.020358022618747569},
            {0.034013390238119661, 0.597210722761849931, 0.184387943500015203, 0.184387943500015203, 0.020358022618747569},
            {0.034013390238119661, 0.184387943500015203, 0.597210722761849931, 0.184387943500015203, 0.020358022618747569},
            {0.184387943500015203, 0.034013390238119661, 0.184387943500015203, 0.597210722761849931, 0.020358022618747569},
            {0.597210722761849931, 0.034013390238119661, 0.184387943500015203, 0.184387943500015203, 0.020358022618747569},
            {0.184387943500015203, 0.034013390238119661, 0.597210722761849931, 0.184387943500015203, 0.020358022618747569},
            {0.184387943500015203, 0.184387943500015203, 0.034013390238119661, 0.597210722761849931, 0.020358022618747569},
            {0.597210722761849931, 0.184387943500015203, 0.034013390238119661, 0.184387943500015203, 0.020358022618747569},
            {0.184387943500015203, 0.597210722761849931, 0.034013390238119661, 0.184387943500015203, 0.020358022618747569},
            }}
        )
    {
    }
};

template<>
class GaussLegendreTet<10> : public internal::GaussLegendreTetSuper<81>{
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
    : internal::GaussLegendreTetSuper<81>(std::array<GLPointTet, 81>
            {{
            {0.25, 0.25, 0.25, 0.25, 0.051650892256097467},
            {0.072551757417437801, 0.072551757417437801, 0.072551757417437801, 0.782344727747686596, 0.013467972671660261},
            {0.782344727747686596, 0.072551757417437801, 0.072551757417437801, 0.072551757417437801, 0.013467972671660261},
            {0.072551757417437801, 0.782344727747686596, 0.072551757417437801, 0.072551757417437801, 0.013467972671660261},
            {0.072551757417437801, 0.072551757417437801, 0.782344727747686596, 0.072551757417437801, 0.013467972671660261},
            {0.307895947866922893, 0.307895947866922893, 0.307895947866922893, 0.076312156399231318, 0.026466919249689534},
            {0.076312156399231318, 0.307895947866922893, 0.307895947866922893, 0.307895947866922893, 0.026466919249689534},
            {0.307895947866922893, 0.076312156399231318, 0.307895947866922893, 0.307895947866922893, 0.026466919249689534},
            {0.307895947866922893, 0.307895947866922893, 0.076312156399231318, 0.307895947866922893, 0.026466919249689534},
            {0.476280860392779071, 0.023719139607220928, 0.023719139607220928, 0.476280860392779071, 0.005107047828143525},
            {0.023719139607220928, 0.476280860392779071, 0.023719139607220928, 0.476280860392779071, 0.005107047828143525},
            {0.476280860392779071, 0.476280860392779071, 0.023719139607220928, 0.023719139607220928, 0.005107047828143525},
            {0.023719139607220928, 0.023719139607220928, 0.476280860392779071, 0.476280860392779071, 0.005107047828143525},
            {0.476280860392779071, 0.023719139607220928, 0.476280860392779071, 0.023719139607220928, 0.005107047828143525},
            {0.023719139607220928, 0.476280860392779071, 0.476280860392779071, 0.023719139607220928, 0.005107047828143525},
            {0.096579773069107259, 0.403420226930892740, 0.403420226930892740, 0.096579773069107259, 0.023649996158245595},
            {0.403420226930892740, 0.096579773069107259, 0.403420226930892740, 0.096579773069107259, 0.023649996158245595},
            {0.096579773069107259, 0.096579773069107259, 0.403420226930892740, 0.403420226930892740, 0.023649996158245595},
            {0.403420226930892740, 0.403420226930892740, 0.096579773069107259, 0.096579773069107259, 0.023649996158245595},
            {0.096579773069107259, 0.403420226930892740, 0.096579773069107259, 0.403420226930892740, 0.023649996158245595},
            {0.403420226930892740, 0.096579773069107259, 0.096579773069107259, 0.403420226930892740, 0.023649996158245595},
            {0.073854181018486210, 0.201738606847679876, 0.201738606847679876, 0.522668605286154036, 0.025317366650968218},
            {0.201738606847679876, 0.073854181018486210, 0.201738606847679876, 0.522668605286154036, 0.025317366650968218},
            {0.201738606847679876, 0.201738606847679876, 0.073854181018486210, 0.522668605286154036, 0.025317366650968218},
            {0.522668605286154036, 0.201738606847679876, 0.201738606847679876, 0.073854181018486210, 0.025317366650968218},
            {0.522668605286154036, 0.073854181018486210, 0.201738606847679876, 0.201738606847679876, 0.025317366650968218},
            {0.522668605286154036, 0.201738606847679876, 0.073854181018486210, 0.201738606847679876, 0.025317366650968218},
            {0.201738606847679876, 0.522668605286154036, 0.201738606847679876, 0.073854181018486210, 0.025317366650968218},
            {0.073854181018486210, 0.522668605286154036, 0.201738606847679876, 0.201738606847679876, 0.025317366650968218},
            {0.201738606847679876, 0.522668605286154036, 0.073854181018486210, 0.201738606847679876, 0.025317366650968218},
            {0.201738606847679876, 0.201738606847679876, 0.522668605286154036, 0.073854181018486210, 0.025317366650968218},
            {0.073854181018486210, 0.201738606847679876, 0.522668605286154036, 0.201738606847679876, 0.025317366650968218},
            {0.201738606847679876, 0.073854181018486210, 0.522668605286154036, 0.201738606847679876, 0.025317366650968218},
            {0.687504855883039547, 0.154839700743861806, 0.154839700743861806, 0.002815742629236839, 0.005416263407447926},
            {0.154839700743861806, 0.687504855883039547, 0.154839700743861806, 0.002815742629236839, 0.005416263407447926},
            {0.154839700743861806, 0.154839700743861806, 0.687504855883039547, 0.002815742629236839, 0.005416263407447926},
            {0.002815742629236839, 0.154839700743861806, 0.154839700743861806, 0.687504855883039547, 0.005416263407447926},
            {0.002815742629236839, 0.687504855883039547, 0.154839700743861806, 0.154839700743861806, 0.005416263407447926},
            {0.002815742629236839, 0.154839700743861806, 0.687504855883039547, 0.154839700743861806, 0.005416263407447926},
            {0.154839700743861806, 0.002815742629236839, 0.154839700743861806, 0.687504855883039547, 0.005416263407447926},
            {0.687504855883039547, 0.002815742629236839, 0.154839700743861806, 0.154839700743861806, 0.005416263407447926},
            {0.154839700743861806, 0.002815742629236839, 0.687504855883039547, 0.154839700743861806, 0.005416263407447926},
            {0.154839700743861806, 0.154839700743861806, 0.002815742629236839, 0.687504855883039547, 0.005416263407447926},
            {0.687504855883039547, 0.154839700743861806, 0.002815742629236839, 0.154839700743861806, 0.005416263407447926},
            {0.154839700743861806, 0.687504855883039547, 0.002815742629236839, 0.154839700743861806, 0.005416263407447926},
            {0.174143716493932110, 0.412599150499363910, 0.412599150499363910, 0.000657982507340068, 0.008153659012054226},
            {0.412599150499363910, 0.174143716493932110, 0.412599150499363910, 0.000657982507340068, 0.008153659012054226},
            {0.412599150499363910, 0.412599150499363910, 0.174143716493932110, 0.000657982507340068, 0.008153659012054226},
            {0.000657982507340068, 0.412599150499363910, 0.412599150499363910, 0.174143716493932110, 0.008153659012054226},
            {0.000657982507340068, 0.174143716493932110, 0.412599150499363910, 0.412599150499363910, 0.008153659012054226},
            {0.000657982507340068, 0.412599150499363910, 0.174143716493932110, 0.412599150499363910, 0.008153659012054226},
            {0.412599150499363910, 0.000657982507340068, 0.412599150499363910, 0.174143716493932110, 0.008153659012054226},
            {0.174143716493932110, 0.000657982507340068, 0.412599150499363910, 0.412599150499363910, 0.008153659012054226},
            {0.412599150499363910, 0.000657982507340068, 0.174143716493932110, 0.412599150499363910, 0.008153659012054226},
            {0.412599150499363910, 0.412599150499363910, 0.000657982507340068, 0.174143716493932110, 0.008153659012054226},
            {0.174143716493932110, 0.412599150499363910, 0.000657982507340068, 0.412599150499363910, 0.008153659012054226},
            {0.412599150499363910, 0.174143716493932110, 0.000657982507340068, 0.412599150499363910, 0.008153659012054226},
            {0.265757506562719879, 0.038890727557806512, 0.038890727557806512, 0.656461038321667095, 0.011125971759612777},
            {0.038890727557806512, 0.265757506562719879, 0.038890727557806512, 0.656461038321667095, 0.011125971759612777},
            {0.038890727557806512, 0.038890727557806512, 0.265757506562719879, 0.656461038321667095, 0.011125971759612777},
            {0.656461038321667095, 0.038890727557806512, 0.038890727557806512, 0.265757506562719879, 0.011125971759612777},
            {0.656461038321667095, 0.265757506562719879, 0.038890727557806512, 0.038890727557806512, 0.011125971759612777},
            {0.656461038321667095, 0.038890727557806512, 0.265757506562719879, 0.038890727557806512, 0.011125971759612777},
            {0.038890727557806512, 0.656461038321667095, 0.038890727557806512, 0.265757506562719879, 0.011125971759612777},
            {0.265757506562719879, 0.656461038321667095, 0.038890727557806512, 0.038890727557806512, 0.011125971759612777},
            {0.038890727557806512, 0.656461038321667095, 0.265757506562719879, 0.038890727557806512, 0.011125971759612777},
            {0.038890727557806512, 0.038890727557806512, 0.656461038321667095, 0.265757506562719879, 0.011125971759612777},
            {0.265757506562719879, 0.038890727557806512, 0.656461038321667095, 0.038890727557806512, 0.011125971759612777},
            {0.038890727557806512, 0.265757506562719879, 0.656461038321667095, 0.038890727557806512, 0.011125971759612777},
            {0.093241349035258407, 0.005604847951021337, 0.005604847951021337, 0.895548955062698917, 0.001325678848264235},
            {0.005604847951021337, 0.093241349035258407, 0.005604847951021337, 0.895548955062698917, 0.001325678848264235},
            {0.005604847951021337, 0.005604847951021337, 0.093241349035258407, 0.895548955062698917, 0.001325678848264235},
            {0.895548955062698917, 0.005604847951021337, 0.005604847951021337, 0.093241349035258407, 0.001325678848264235},
            {0.895548955062698917, 0.093241349035258407, 0.005604847951021337, 0.005604847951021337, 0.001325678848264235},
            {0.895548955062698917, 0.005604847951021337, 0.093241349035258407, 0.005604847951021337, 0.001325678848264235},
            {0.005604847951021337, 0.895548955062698917, 0.005604847951021337, 0.093241349035258407, 0.001325678848264235},
            {0.093241349035258407, 0.895548955062698917, 0.005604847951021337, 0.005604847951021337, 0.001325678848264235},
            {0.005604847951021337, 0.895548955062698917, 0.093241349035258407, 0.005604847951021337, 0.001325678848264235},
            {0.005604847951021337, 0.005604847951021337, 0.895548955062698917, 0.093241349035258407, 0.001325678848264235},
            {0.093241349035258407, 0.005604847951021337, 0.895548955062698917, 0.005604847951021337, 0.001325678848264235},
            {0.005604847951021337, 0.093241349035258407, 0.895548955062698917, 0.005604847951021337, 0.001325678848264235},
            }}
        )
    {
    }
};

template<>
class GaussLegendreTet<11> : public internal::GaussLegendreTetSuper<110>{
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
    : internal::GaussLegendreTetSuper<110>(std::array<GLPointTet, 110>
            {{
            {0.299313092399342490, 0.299313092399342490, 0.299313092399342490, 0.102060722801972527, 0.014689343586130748},
            {0.102060722801972527, 0.299313092399342490, 0.299313092399342490, 0.299313092399342490, 0.014689343586130748},
            {0.299313092399342490, 0.102060722801972527, 0.299313092399342490, 0.299313092399342490, 0.014689343586130748},
            {0.299313092399342490, 0.299313092399342490, 0.102060722801972527, 0.299313092399342490, 0.014689343586130748},
            {0.032632843965181795, 0.032632843965181795, 0.032632843965181795, 0.902101468104454613, 0.003380869530621680},
            {0.902101468104454613, 0.032632843965181795, 0.032632843965181795, 0.032632843965181795, 0.003380869530621680},
            {0.032632843965181795, 0.902101468104454613, 0.032632843965181795, 0.032632843965181795, 0.003380869530621680},
            {0.032632843965181795, 0.032632843965181795, 0.902101468104454613, 0.032632843965181795, 0.003380869530621680},
            {0.325795584515310834, 0.174204415484689165, 0.174204415484689165, 0.325795584515310834, 0.016766768677225626},
            {0.174204415484689165, 0.325795584515310834, 0.174204415484689165, 0.325795584515310834, 0.016766768677225626},
            {0.325795584515310834, 0.325795584515310834, 0.174204415484689165, 0.174204415484689165, 0.016766768677225626},
            {0.174204415484689165, 0.174204415484689165, 0.325795584515310834, 0.325795584515310834, 0.016766768677225626},
            {0.325795584515310834, 0.174204415484689165, 0.325795584515310834, 0.174204415484689165, 0.016766768677225626},
            {0.174204415484689165, 0.325795584515310834, 0.325795584515310834, 0.174204415484689165, 0.016766768677225626},
            {0.098784624141625061, 0.401215375858374938, 0.401215375858374938, 0.098784624141625061, 0.019323640347758510},
            {0.401215375858374938, 0.098784624141625061, 0.401215375858374938, 0.098784624141625061, 0.019323640347758510},
            {0.098784624141625061, 0.098784624141625061, 0.401215375858374938, 0.401215375858374938, 0.019323640347758510},
            {0.401215375858374938, 0.401215375858374938, 0.098784624141625061, 0.098784624141625061, 0.019323640347758510},
            {0.098784624141625061, 0.401215375858374938, 0.098784624141625061, 0.401215375858374938, 0.019323640347758510},
            {0.401215375858374938, 0.098784624141625061, 0.098784624141625061, 0.401215375858374938, 0.019323640347758510},
            {0.005701905122700351, 0.494298094877299648, 0.494298094877299648, 0.005701905122700351, 0.001653253610375796},
            {0.494298094877299648, 0.005701905122700351, 0.494298094877299648, 0.005701905122700351, 0.001653253610375796},
            {0.005701905122700351, 0.005701905122700351, 0.494298094877299648, 0.494298094877299648, 0.001653253610375796},
            {0.494298094877299648, 0.494298094877299648, 0.005701905122700351, 0.005701905122700351, 0.001653253610375796},
            {0.005701905122700351, 0.494298094877299648, 0.005701905122700351, 0.494298094877299648, 0.001653253610375796},
            {0.494298094877299648, 0.005701905122700351, 0.005701905122700351, 0.494298094877299648, 0.001653253610375796},
            {0.727330346346383446, 0.123120654974722356, 0.123120654974722356, 0.026428343704171840, 0.009799971247032873},
            {0.123120654974722356, 0.727330346346383446, 0.123120654974722356, 0.026428343704171840, 0.009799971247032873},
            {0.123120654974722356, 0.123120654974722356, 0.727330346346383446, 0.026428343704171840, 0.009799971247032873},
            {0.026428343704171840, 0.123120654974722356, 0.123120654974722356, 0.727330346346383446, 0.009799971247032873},
            {0.026428343704171840, 0.727330346346383446, 0.123120654974722356, 0.123120654974722356, 0.009799971247032873},
            {0.026428343704171840, 0.123120654974722356, 0.727330346346383446, 0.123120654974722356, 0.009799971247032873},
            {0.123120654974722356, 0.026428343704171840, 0.123120654974722356, 0.727330346346383446, 0.009799971247032873},
            {0.727330346346383446, 0.026428343704171840, 0.123120654974722356, 0.123120654974722356, 0.009799971247032873},
            {0.123120654974722356, 0.026428343704171840, 0.727330346346383446, 0.123120654974722356, 0.009799971247032873},
            {0.123120654974722356, 0.123120654974722356, 0.026428343704171840, 0.727330346346383446, 0.009799971247032873},
            {0.727330346346383446, 0.123120654974722356, 0.026428343704171840, 0.123120654974722356, 0.009799971247032873},
            {0.123120654974722356, 0.727330346346383446, 0.026428343704171840, 0.123120654974722356, 0.009799971247032873},
            {0.170228073489264195, 0.013003870485590236, 0.013003870485590236, 0.803764185539555330, 0.002036278491304825},
            {0.013003870485590236, 0.170228073489264195, 0.013003870485590236, 0.803764185539555330, 0.002036278491304825},
            {0.013003870485590236, 0.013003870485590236, 0.170228073489264195, 0.803764185539555330, 0.002036278491304825},
            {0.803764185539555330, 0.013003870485590236, 0.013003870485590236, 0.170228073489264195, 0.002036278491304825},
            {0.803764185539555330, 0.170228073489264195, 0.013003870485590236, 0.013003870485590236, 0.002036278491304825},
            {0.803764185539555330, 0.013003870485590236, 0.170228073489264195, 0.013003870485590236, 0.002036278491304825},
            {0.013003870485590236, 0.803764185539555330, 0.013003870485590236, 0.170228073489264195, 0.002036278491304825},
            {0.170228073489264195, 0.803764185539555330, 0.013003870485590236, 0.013003870485590236, 0.002036278491304825},
            {0.013003870485590236, 0.803764185539555330, 0.170228073489264195, 0.013003870485590236, 0.002036278491304825},
            {0.013003870485590236, 0.013003870485590236, 0.803764185539555330, 0.170228073489264195, 0.002036278491304825},
            {0.170228073489264195, 0.013003870485590236, 0.803764185539555330, 0.013003870485590236, 0.002036278491304825},
            {0.013003870485590236, 0.170228073489264195, 0.803764185539555330, 0.013003870485590236, 0.002036278491304825},
            {0.618259809789115788, 0.042935401171925739, 0.042935401171925739, 0.295869387867032733, 0.008207613752916718},
            {0.042935401171925739, 0.618259809789115788, 0.042935401171925739, 0.295869387867032733, 0.008207613752916718},
            {0.042935401171925739, 0.042935401171925739, 0.618259809789115788, 0.295869387867032733, 0.008207613752916718},
            {0.295869387867032733, 0.042935401171925739, 0.042935401171925739, 0.618259809789115788, 0.008207613752916718},
            {0.295869387867032733, 0.618259809789115788, 0.042935401171925739, 0.042935401171925739, 0.008207613752916718},
            {0.295869387867032733, 0.042935401171925739, 0.618259809789115788, 0.042935401171925739, 0.008207613752916718},
            {0.042935401171925739, 0.295869387867032733, 0.042935401171925739, 0.618259809789115788, 0.008207613752916718},
            {0.618259809789115788, 0.295869387867032733, 0.042935401171925739, 0.042935401171925739, 0.008207613752916718},
            {0.042935401171925739, 0.295869387867032733, 0.618259809789115788, 0.042935401171925739, 0.008207613752916718},
            {0.042935401171925739, 0.042935401171925739, 0.295869387867032733, 0.618259809789115788, 0.008207613752916718},
            {0.618259809789115788, 0.042935401171925739, 0.295869387867032733, 0.042935401171925739, 0.008207613752916718},
            {0.042935401171925739, 0.618259809789115788, 0.295869387867032733, 0.042935401171925739, 0.008207613752916718},
            {0.539900116855701610, 0.184424032047077767, 0.184424032047077767, 0.091251819050142855, 0.015104987366663492},
            {0.184424032047077767, 0.539900116855701610, 0.184424032047077767, 0.091251819050142855, 0.015104987366663492},
            {0.184424032047077767, 0.184424032047077767, 0.539900116855701610, 0.091251819050142855, 0.015104987366663492},
            {0.091251819050142855, 0.184424032047077767, 0.184424032047077767, 0.539900116855701610, 0.015104987366663492},
            {0.091251819050142855, 0.539900116855701610, 0.184424032047077767, 0.184424032047077767, 0.015104987366663492},
            {0.091251819050142855, 0.184424032047077767, 0.539900116855701610, 0.184424032047077767, 0.015104987366663492},
            {0.184424032047077767, 0.091251819050142855, 0.184424032047077767, 0.539900116855701610, 0.015104987366663492},
            {0.539900116855701610, 0.091251819050142855, 0.184424032047077767, 0.184424032047077767, 0.015104987366663492},
            {0.184424032047077767, 0.091251819050142855, 0.539900116855701610, 0.184424032047077767, 0.015104987366663492},
            {0.184424032047077767, 0.184424032047077767, 0.091251819050142855, 0.539900116855701610, 0.015104987366663492},
            {0.539900116855701610, 0.184424032047077767, 0.091251819050142855, 0.184424032047077767, 0.015104987366663492},
            {0.184424032047077767, 0.539900116855701610, 0.091251819050142855, 0.184424032047077767, 0.015104987366663492},
            {0.455435319621959524, 0.259169407464282648, 0.259169407464282648, 0.026225865449475177, 0.012742392666106386},
            {0.259169407464282648, 0.455435319621959524, 0.259169407464282648, 0.026225865449475177, 0.012742392666106386},
            {0.259169407464282648, 0.259169407464282648, 0.455435319621959524, 0.026225865449475177, 0.012742392666106386},
            {0.026225865449475177, 0.259169407464282648, 0.259169407464282648, 0.455435319621959524, 0.012742392666106386},
            {0.026225865449475177, 0.455435319621959524, 0.259169407464282648, 0.259169407464282648, 0.012742392666106386},
            {0.026225865449475177, 0.259169407464282648, 0.455435319621959524, 0.259169407464282648, 0.012742392666106386},
            {0.259169407464282648, 0.026225865449475177, 0.259169407464282648, 0.455435319621959524, 0.012742392666106386},
            {0.455435319621959524, 0.026225865449475177, 0.259169407464282648, 0.259169407464282648, 0.012742392666106386},
            {0.259169407464282648, 0.026225865449475177, 0.455435319621959524, 0.259169407464282648, 0.012742392666106386},
            {0.259169407464282648, 0.259169407464282648, 0.026225865449475177, 0.455435319621959524, 0.012742392666106386},
            {0.455435319621959524, 0.259169407464282648, 0.026225865449475177, 0.259169407464282648, 0.012742392666106386},
            {0.259169407464282648, 0.455435319621959524, 0.026225865449475177, 0.259169407464282648, 0.012742392666106386},
            {0.010977572365960306, 0.535142000931459522, 0.112428560506761202, 0.341451866195818968, 0.005273427059689129},
            {0.010977572365960306, 0.112428560506761202, 0.535142000931459522, 0.341451866195818968, 0.005273427059689129},
            {0.535142000931459522, 0.010977572365960306, 0.112428560506761202, 0.341451866195818968, 0.005273427059689129},
            {0.112428560506761202, 0.010977572365960306, 0.535142000931459522, 0.341451866195818968, 0.005273427059689129},
            {0.535142000931459522, 0.112428560506761202, 0.010977572365960306, 0.341451866195818968, 0.005273427059689129},
            {0.112428560506761202, 0.535142000931459522, 0.010977572365960306, 0.341451866195818968, 0.005273427059689129},
            {0.341451866195818968, 0.535142000931459522, 0.112428560506761202, 0.010977572365960306, 0.005273427059689129},
            {0.341451866195818968, 0.112428560506761202, 0.535142000931459522, 0.010977572365960306, 0.005273427059689129},
            {0.341451866195818968, 0.010977572365960306, 0.112428560506761202, 0.535142000931459522, 0.005273427059689129},
            {0.341451866195818968, 0.010977572365960306, 0.535142000931459522, 0.112428560506761202, 0.005273427059689129},
            {0.341451866195818968, 0.112428560506761202, 0.010977572365960306, 0.535142000931459522, 0.005273427059689129},
            {0.341451866195818968, 0.535142000931459522, 0.010977572365960306, 0.112428560506761202, 0.005273427059689129},
            {0.535142000931459522, 0.341451866195818968, 0.112428560506761202, 0.010977572365960306, 0.005273427059689129},
            {0.112428560506761202, 0.341451866195818968, 0.535142000931459522, 0.010977572365960306, 0.005273427059689129},
            {0.010977572365960306, 0.341451866195818968, 0.112428560506761202, 0.535142000931459522, 0.005273427059689129},
            {0.010977572365960306, 0.341451866195818968, 0.535142000931459522, 0.112428560506761202, 0.005273427059689129},
            {0.112428560506761202, 0.341451866195818968, 0.010977572365960306, 0.535142000931459522, 0.005273427059689129},
            {0.535142000931459522, 0.341451866195818968, 0.010977572365960306, 0.112428560506761202, 0.005273427059689129},
            {0.535142000931459522, 0.112428560506761202, 0.341451866195818968, 0.010977572365960306, 0.005273427059689129},
            {0.112428560506761202, 0.535142000931459522, 0.341451866195818968, 0.010977572365960306, 0.005273427059689129},
            {0.010977572365960306, 0.112428560506761202, 0.341451866195818968, 0.535142000931459522, 0.005273427059689129},
            {0.010977572365960306, 0.535142000931459522, 0.341451866195818968, 0.112428560506761202, 0.005273427059689129},
            {0.112428560506761202, 0.010977572365960306, 0.341451866195818968, 0.535142000931459522, 0.005273427059689129},
            {0.535142000931459522, 0.010977572365960306, 0.341451866195818968, 0.112428560506761202, 0.005273427059689129},
            }}
        )
    {
    }
};

}

#endif
