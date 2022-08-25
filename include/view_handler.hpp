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
#include <gmsh.h>
#include <string>
#include <vector>

class ViewHandler{
    public:
    enum ViewType{
        ELEMENTAL,
        NODAL,
        VECTOR,
        TENSOR
    };

    enum DataType{
        STRESS,
        MATERIAL,
        DISPLACEMENT,
        OTHER
    };

    ViewHandler(const Meshing* const mesh, const std::string& model_name, const std::string& view_name, const ViewType view_type, const DataType data_type, const size_t view_id);

    /**
     * The `geometries` parameter is currently unused for nodal views.
     */
    void update_view(const std::vector<double>& data, const std::vector<size_t>& geometries = std::vector<size_t>()) const;

    const size_t view_id;
    const std::string model_name;
    const ViewType view_type;
    const DataType data_type;

    private:
    const Meshing* const mesh;
    const size_t elem_num;
    const size_t node_num;
    const size_t mat_color_num;

    size_t get_number_of_elements() const;
    size_t get_number_of_nodes() const;
    size_t get_number_of_material_colors() const;

    inline void update_elemental(const std::vector<double>& data, const std::vector<size_t>& tags) const{
        this->gmsh_update_view(data, tags, 1, "ElementData");
    }
    inline void update_vector(const std::vector<double>& data, const std::vector<size_t>& tags) const{
        this->update_nodal(data, tags, 3);
    }
    inline void update_tensor(const std::vector<double>& data, const std::vector<size_t>& tags) const{
        this->update_nodal(data, tags, 9);
    }
    inline void update_nodal(const std::vector<double>& data, const std::vector<size_t>& tags, const int num_components = 1) const{
        this->gmsh_update_view(data, tags, num_components, "NodeData");
    }
    inline void gmsh_update_view(const std::vector<double>& data, const std::vector<size_t>& tags, const int num_components, const std::string& data_type) const{
        gmsh::view::addHomogeneousModelData(this->view_id, 0, this->model_name, data_type, tags, data, 0, num_components);
    }
};


#endif
