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

#include <dmumps_c.h>
#include "finite_element/mumps_solver.hpp"
#include "logger.hpp"

namespace finite_element{

std::vector<double> MUMPSSolver::calculate_displacements(const Meshing* const mesh, std::vector<double> load, const std::vector<double>& density, double pc){
    if(this->current_step == 0){

        this->generate_K(mesh, density, pc);

        logger::quick_log("Decomposing...");

        logger::quick_log("Done.");
    }
    logger::quick_log("Calculating displacements...");


    logger::quick_log("Done.");

    this->current_step = (this->current_step + 1) % this->steps;
   
    return load; 
}

}
