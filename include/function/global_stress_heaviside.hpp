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

#ifndef FUNCTION_GLOBAL_STRESS_HEAVISIDE_HPP
#define FUNCTION_GLOBAL_STRESS_HEAVISIDE_HPP

#include "function.hpp"
#include "meshing.hpp"
#include "finite_element.hpp"

namespace function{

class GlobalStressHeaviside : public DensityBasedFunction{
    public:
    GlobalStressHeaviside(const Meshing* const mesh, FiniteElement* fem, double max_stress, double C, double pc, double pt, double psiK, double psiS);

    virtual ~GlobalStressHeaviside() = default;

    virtual double calculate(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x) override;
    virtual double calculate_with_gradient(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad) override;
    virtual size_t additional_steps() const override{
        return 1;
    }

    private:
    const Meshing* const mesh;
    FiniteElement* fem;
    const double max_stress;
    const double C;
    const double pc;
    const double pt;
    const double psiK;
    const double psiS;
    const double vm_eps = 1e-14;
    const double rho_eps = 0.2;
    size_t elem_number = 0;

    inline double heaviside(const double x) const{
        return std::atan(C*(x - max_stress))/M_PI + 0.5;
    }
    inline double heaviside_grad(const double x) const{
        const double X = C*(x - max_stress);
        const double datan = 1.0/(1.0 + X*X);
        return C*datan/M_PI;
    }

    inline double relaxed_rho(const double x) const{
        return x/(rho_eps*(1.0 - x) + x);
    }
    inline double relaxed_rho_grad(const double x) const{
        const double d = (rho_eps*(1.0 - x) + x);
        return rho_eps/(d*d);
    }
};

}

#endif
