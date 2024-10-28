/*
 *   Copyright (C) 2024 Tarcísio Ladeia de Oliveira.
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

#ifndef MESH_FILE_HPP
#define MESH_FILE_HPP

#include "meshing.hpp"

class ProjectData;

namespace meshing{

class MeshFile : public Meshing{
    public:
    MeshFile(const std::vector<std::unique_ptr<Geometry>>& geometries,
         const MeshElementFactory* const elem_type,
         const ProjectData* const proj_data,
         double thickness, const std::string& file_path,
         const std::string& elem_name, bool load);

    void save_mesh(const std::vector<std::unique_ptr<MeshNode>>& node_list,
                   const std::vector<MeshNode*>& boundary_node_list,
                   const std::vector<size_t>& geom_elem_mapping, 
                   const std::vector<size_t>& elem_node_tags, 
                   const std::vector<size_t>& bound_elem_node_tags,
                   std::unordered_map<size_t, MeshNode*>& id_map,
                   std::unordered_map<size_t, size_t>& duplicate_map) const;

    inline std::string get_element_name() const{
        return this->element_name;
    }

    // In this case, loads mesh from file path
    virtual void mesh(const std::vector<Force>& forces, 
                      const std::vector<Support>& supports,
                      std::vector<Spring>& springs) override;

    private:
    std::string element_name;
    std::string file_path;
};

}

#endif
