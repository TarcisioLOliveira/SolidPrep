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

#ifndef VISUALIZATION_HPP
#define VISUALIZATION_HPP

#include <vector>
#include "element.hpp"
#include "meshing.hpp"
#include <gmsh.h>
#include "utils.hpp"

class Visualization{
    public:

    inline void start() const {gmsh::initialize();}

    void load_mesh(Meshing* mesh, utils::ProblemType type);
    void update_view();
    void show();
    inline void hide(){
        this->shown = false;
    }

    inline void end(){
        this->hide();
        gmsh::finalize();
    }

    private:
    const std::string MODEL_NAME = "loaded_model";
    const std::string STRESS_VIEW = "stress_view";
    bool shown = false;
    Meshing* mesh = nullptr;
    int tag = 0;
};

#endif
