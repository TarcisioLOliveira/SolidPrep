/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
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

#ifndef VISUALIZATION_HPP
#define VISUALIZATION_HPP

#include <vector>
#include "element.hpp"
#include "meshing.hpp"
#include <gmsh.h>
#include "utils.hpp"
#include "view_handler.hpp"

class Visualization{
    public:

    inline void start() const {gmsh::initialize();}

    void load_mesh(Meshing* mesh, utils::ProblemType type);
    ViewHandler* add_view(const std::string& view_name, ViewHandler::ViewType view_type, ViewHandler::DataType data_type);

    void show();
    inline void hide(){
        this->shown = false;
    }
    inline void redraw() const{
        if(this->shown){
            gmsh::graphics::draw();
        }
    }

    void wait();
    inline void end(){
        this->handler_list.clear();
        this->hide();
        gmsh::finalize();
    }

    private:
    const std::string MODEL_NAME = "loaded_model";
    bool shown = false;
    Meshing* mesh = nullptr;
    int mesh_tag = 0;
    int last_view_tag = 0;
    utils::ProblemType type = utils::PROBLEM_TYPE_2D;
    std::vector<std::unique_ptr<ViewHandler>> handler_list;
};

#endif
