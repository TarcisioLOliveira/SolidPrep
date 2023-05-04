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

#ifndef VIEW_HANDLER_HPP
#define VIEW_HANDLER_HPP

#include "meshing.hpp"
#include "spview.hpp"
#include <string>
#include <vector>

class ViewHandler{
    public:

    ViewHandler(const Meshing* const mesh, spview::Server* server, const std::string& view_name, const spview::defs::ViewType view_type, const spview::defs::DataType data_type, utils::ProblemType problem_type, const size_t view_id);
    ~ViewHandler(){
        this->remove_view();
    }

    /**
     * The `geometries` parameter is currently unused for nodal views.
     */
    void update_view(const std::vector<double>& data, const std::vector<size_t>& geometries = std::vector<size_t>()) const;

    inline bool is_removed() const{
        return this->removed;
    }
    inline void remove_view() const{
        if(!this->removed){
            this->server->remove_view(this->view_id);
        }
    }

    const spview::defs::ViewType view_type;
    const spview::defs::DataType data_type;

    private:
    const size_t view_id;
    const Meshing* const mesh;
    const size_t elem_num;
    const size_t node_num;
    const size_t mat_color_num;
    utils::ProblemType problem_type;
    spview::Server* const server;

    bool removed = false;

    size_t get_number_of_elements() const;
    size_t get_number_of_nodes() const;
    size_t get_number_of_material_colors() const;

};


#endif
