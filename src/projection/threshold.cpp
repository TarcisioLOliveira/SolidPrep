/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
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
#include "projection/threshold.hpp"
#include "project_data.hpp"

namespace projection{

Threshold::Threshold(Parameter beta, double eta):
    beta(std::move(beta)), eta(eta){}
Threshold::Threshold(const projspec::DataMap& data):
    beta(data.proj->projection_parameters),
    eta(data.get_double("eta")){

}

void Threshold::update(const size_t iteration){
    if(iteration > 0 && iteration % beta.iteration_step == 0 && beta.value < beta.final_value){
        beta.value = std::min(beta.final_value, beta.value*beta.value_step);
        logger::quick_log("beta: ", beta.value);
    }
}

void Threshold::project_densities(std::vector<double>& new_x) const{
    const double b = this->beta.value;
    #pragma omp parallel for
    for(auto& x:new_x){
        x = (std::tanh(b*eta) + std::tanh(b*(x-eta)))/(std::tanh(b*eta) + std::tanh(b*(1.0-eta)));
    }
}

void Threshold::project_gradient(std::vector<double>& new_df, const std::vector<double>& new_x) const{
    const double b = this->beta.value;
    #pragma omp parallel for
    for(size_t i = 0; i < new_df.size(); ++i){
        auto& df = new_df[i];
        const auto x = new_x[i];
        const double t = std::tanh(b*(x-eta));
        df *= ((1-t*t)*b)/(std::tanh(b*eta) + std::tanh(b*(1.0-eta)));
    }
}

using namespace projspec;
const bool Threshold::reg = Factory<Projection>::add(
    [](const DataMap& data){
        return std::make_unique<Threshold>(data);
    },
    ObjectRequirements{
        "threshold",
        {
            DataEntry{.name = "eta", .type = TYPE_DOUBLE, .required = true}
        }
    }
);

}
