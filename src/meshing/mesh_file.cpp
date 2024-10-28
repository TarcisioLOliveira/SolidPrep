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

#include <fstream>
#include "meshing/mesh_file.hpp"

namespace meshing{

MeshFile::MeshFile(const std::vector<std::unique_ptr<Geometry>>& geometries,
         const MeshElementFactory* const elem_type,
         const ProjectData* const proj_data,
         double thickness, const std::string& file_path, 
         const std::string& elem_name, bool load):
    Meshing(geometries, elem_type, proj_data, thickness),
    element_name(elem_name), file_path(file_path)
{
    if(load){
        std::ifstream file(file_path, std::ios::in);
        std::getline(file, this->element_name);
        file.close();
    }
}

void MeshFile::save_mesh(const std::vector<std::unique_ptr<MeshNode>>& node_list,
                         const std::vector<MeshNode*>& boundary_node_list,
                         const std::vector<size_t>& geom_elem_mapping, 
                         const std::vector<size_t>& elem_node_tags, 
                         const std::vector<size_t>& bound_elem_node_tags,
                         std::unordered_map<size_t, MeshNode*>& id_map,
                         std::unordered_map<size_t, size_t>& duplicate_map) const{
    /**
     * Save order:
     *
     * - Element name
     * - Node positions (+ vector size)
     *   Boundary node IDs (+ vector size)
     * - Geometry element mapping (+ vector size)
     * - Element node tags (+ vector size)
     * - Boundary element node tags (+ vector size)
     * - ID map (old id/new id pairs)
     * - Duplicate map (id pairs)
     *
     * deduplicate is always true if number of geometries
     * is greater than 1
     * (remove variable, maybe?)
     * - from interface being rigid or not, this is then
     *   considered for deciding if deduplication will actually
     *   happen
     *
     * boundary_condition_inside is always false
     * (this is being phased out)
     */

    std::ofstream file(file_path, std::ios::out);
    file << this->element_name << std::endl;
    file << std::endl;

    file << "Nodes " << node_list.size() << std::endl;
    for(const auto& n:node_list){
        file << n->point.X() << " " << n->point.Y() << " " << n->point.Z() << std::endl;
    }
    file << std::endl;

    file << "Boundary nodes " << boundary_node_list.size() << std::endl;
    for(const auto& n:boundary_node_list){
        file << n->id << " ";
    }
    file << std::endl << std::endl;

    file << "Geom elem mapping " << geom_elem_mapping.size() << std::endl;
    for(const auto& n:geom_elem_mapping){
        file << n << " ";
    }
    file << std::endl << std::endl;

    file << "Elem node tags " << elem_node_tags.size() << std::endl;
    for(const auto& n:elem_node_tags){
        file << n << " ";
    }
    file << std::endl << std::endl;

    file << "Bound elem node tags " << bound_elem_node_tags.size() << std::endl;
    for(const auto& n:bound_elem_node_tags){
        file << n << " ";
    }
    file << std::endl << std::endl;

    file << "ID map " << id_map.size() << std::endl;
    for(const auto& n:id_map){
        file << n.first << " " << n.second->id << std::endl;
    }
    file << std::endl;

    file << "Duplicate map " << duplicate_map.size() << std::endl;
    for(const auto& n:duplicate_map){
        file << n.first << " " << n.second << std::endl;
    }
    file << std::endl;

    file.close();
}

void MeshFile::mesh(const std::vector<Force>& forces, 
                  const std::vector<Support>& supports,
                  std::vector<Spring>& springs){

    // TODO

}


}
