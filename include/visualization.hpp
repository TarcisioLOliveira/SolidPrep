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
#include "utils.hpp"
#include "view_handler.hpp"
#include "spview.hpp"

class Visualization{
    public:
    Visualization();

    void load_mesh(Meshing* mesh, utils::ProblemType type);
    ViewHandler* add_view(const std::string& view_name, spview::defs::ViewType view_type, spview::defs::DataType data_type);

    inline void start(){
        this->server.start();
    }

    inline void wait(){
        this->server.wait();
    }

    inline void end(){
        this->handler_list.clear();
        this->server.close_client();
    }

    private:
    const std::string MODEL_NAME = "loaded_model";
    bool shown = false;
    Meshing* mesh = nullptr;
    int mesh_tag = 0;
    int last_view_tag = 0;
    utils::ProblemType type = utils::PROBLEM_TYPE_2D;
    spview::Server server;
    std::vector<std::unique_ptr<ViewHandler>> handler_list;

    inline std::string server_name() const{
        srand(std::time(nullptr));
        return "SolidPrep" + std::to_string(rand());
    }

    size_t get_number_of_material_colors() const;
};

#endif
