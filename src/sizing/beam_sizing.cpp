/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   Foobar is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Foobar is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with SolidPrep.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "sizing/beam_sizing.hpp"
#include "project_data.hpp"
#include <lapacke.h>
#include <vector>
#include <cmath>
#include <cstring>
#include <logger.hpp>

namespace sizing{

TopoDS_Shape BeamSizing::run(){
    this->graph.run();
    if(this->data->type == ProjectData::TYPE_2D){
        size_t graph_size = this->graph.size();
        for(size_t i = 0; i < graph_size; ++i){

        }
    }

    // return TopoDS_Shape();
}

}
