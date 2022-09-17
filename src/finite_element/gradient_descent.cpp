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

#include <algorithm>
#include <cblas.h>
#include "finite_element/gradient_descent.hpp"
#include "logger.hpp"

namespace finite_element{

GradientDescent::GradientDescent(const double eps, const double step):
    eps(eps),step(step), displacement(1){}

std::vector<double> GradientDescent::calculate_displacements(const Meshing* const mesh, std::vector<double> load, const std::vector<double>& density, double pc){

    if(this->W == 0 || this->N == 0){
        this->calculate_dimensions(mesh, load);
        for(auto& u:this->displacement){
            u.resize(W,0);
        }
    }

    if(this->current_step == 0){
        this->generate_K(mesh, density, pc);
    }
    logger::quick_log("Done.");
    logger::quick_log("Calculating displacements...");
    logger::quick_log("W: ",W," N: ", N);

    auto& u = this->displacement[this->current_step];
    auto d = load;
    size_t it = 0;
    if(this->step > 0){
        auto comp_abs = [](const double a, const double b){
            return std::abs(a) < std::abs(b);
        };
        while(*std::max_element(d.begin(), d.end(),  comp_abs) > this->eps){
            cblas_dcopy(W, load.data(), 1, d.data(), 1);
            cblas_dsymm(CblasColMajor, CblasLeft, CblasLower, W, N-1, 1.0, this->K.data(), N, u.data(), W, -1.0, d.data(), W);
            cblas_daxpy(W, -step, d.data(), 1, u.data(), 1);
            ++it;
        }
    } else {
        double alpha = 1000;
        double top, bot;
        auto d2 = d;
        while(alpha > this->eps){
            cblas_dcopy(W, load.data(), 1, d.data(), 1);
            cblas_dsymm(CblasColMajor, CblasLeft, CblasLower, W, N-1, 1.0, this->K.data(), N, u.data(), W, -1.0, d.data(), W);
            
            top = cblas_ddot(W, d.data(), 1, d.data(), 1);
            cblas_dsymm(CblasColMajor, CblasLeft, CblasLower, W, N-1, 1.0, this->K.data(), N, d.data(), W, 0.0, d2.data(), W);
            bot = cblas_ddot(W, d.data(), 1, d2.data(), 1);
            alpha = top/bot;

            cblas_daxpy(W, -alpha, d.data(), 1, u.data(), 1);
            ++it;
        }
    }

    logger::quick_log("Required iterations: ", it);
    logger::quick_log("");
    logger::quick_log("Done.");

    this->current_step = (this->current_step + 1) % this->steps;
   
    return load; 
}

}
