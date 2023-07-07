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

#include <cmath>
#include "logger.hpp"
#include "projection/heaviside.hpp"

namespace projection{
Heaviside::Heaviside(Parameter beta, double eta):
    beta(std::move(beta)), eta(eta){}

void Heaviside::update(const size_t iteration){
    if(iteration > 0 && iteration % beta.iteration_step == 0 && beta.value < beta.final_value){
        beta.value = std::min(beta.final_value, beta.value*beta.value_step);
        logger::quick_log("beta: ", beta.value);
    }
}

void Heaviside::project_densities(std::vector<double>& new_x) const{
    const double b = this->beta.value;
    #pragma omp parallel for
    for(auto& x:new_x){
        x = 1.0/(1.0+std::exp(-2*b*(x - eta)));
    }
}

void Heaviside::project_gradient(std::vector<double>& new_df, const std::vector<double>& new_x) const{
    const double b = this->beta.value;
    #pragma omp parallel for
    for(size_t i = 0; i < new_df.size(); ++i){
        auto& df = new_df[i];
        const auto x = new_x[i];
        const double s = 1.0/(1.0+std::exp(-2*b*(x - eta)));
        const double ds = s*(1.0-s)*(2*b);
        df *= ds;
    }
}

}
