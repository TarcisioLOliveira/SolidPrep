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
#include <string>
#include "meshing/mesh_file.hpp"
#include "logger.hpp"
#include "project_data.hpp"

namespace meshing{

MeshFile::MeshFile(const std::vector<std::unique_ptr<Geometry>>& geometries,
         const MeshElementFactory* const elem_type,
         const ProjectData* const proj_data,
         double thickness, const std::string& file_path, 
         const std::string& elem_name):
    Meshing(geometries, elem_type, proj_data, thickness),
    element_name(elem_name), file_path(file_path)
{
}

void MeshFile::save_mesh(const std::vector<std::unique_ptr<MeshNode>>& node_list,
                         const std::vector<MeshNode*>& boundary_node_list,
                         const std::vector<size_t>& geom_elem_mapping, 
                         const std::vector<size_t>& elem_node_tags, 
                         const std::vector<size_t>& bound_elem_node_tags,
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
        file << n->id << " " << n->point.X() << " " << n->point.Y() << " " << n->point.Z() << std::endl;
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
    (void) forces;
    (void) supports;
    (void) springs;

    const size_t dof = this->proj_data->topopt_element->get_dof_per_node();

    std::vector<size_t> geom_elem_mapping, elem_node_tags, bound_elem_node_tags;
    std::unordered_map<size_t, size_t> duplicate_map;
    std::unordered_map<size_t, MeshNode*> id_map;

    std::ifstream file(file_path, std::ios::in);

    std::string line;
    std::string label;
    std::string num;
    std::string expected_label;

    std::getline(file, line);
    logger::log_assert(line == this->element_name, logger::ERROR, "element name was set incorrectly. File: {}, set: {}", line, this->element_name);

    std::getline(file, line);
    std::getline(file, line);
    expected_label = "Nodes";
    label = line.substr(0, expected_label.size());
    logger::log_assert(label == expected_label, logger::ERROR, "unexpected label: \"{}\", expected \"{}\"", label, expected_label);
    num = line.substr(expected_label.size()+1);
    this->node_list.resize(std::stoul(num));
    {
        size_t id;
        gp_Pnt p;
        for(size_t i = 0; i < this->node_list.size(); ++i){
            std::getline(file, num, ' ');
            id = std::stoul(num);
            std::getline(file, num, ' ');
            p.SetX(std::stod(num));
            std::getline(file, num, ' ');
            p.SetY(std::stod(num));
            std::getline(file, num);
            p.SetZ(std::stod(num));

            this->node_list[i] = std::make_unique<MeshNode>(p, id, dof);
            id_map.emplace(id, this->node_list[i].get());
        }
    }

    std::getline(file, line);
    std::getline(file, line);
    expected_label = "Boundary nodes";
    label = line.substr(0, expected_label.size());
    logger::log_assert(label == expected_label, logger::ERROR, "unexpected label: \"{}\", expected \"{}\"", label, expected_label);
    num = line.substr(expected_label.size()+1);
    this->boundary_node_list.resize(std::stoul(num));
    {
        size_t id;
        for(size_t i = 0; i < this->boundary_node_list.size(); ++i){
            std::getline(file, num, ' ');
            id = std::stoul(num);

            this->boundary_node_list[i] = id_map.at(id);
        }
    }
    
    std::getline(file, line);
    std::getline(file, line);
    std::getline(file, line);
    expected_label = "Geom elem mapping";
    label = line.substr(0, expected_label.size());
    logger::log_assert(label == expected_label, logger::ERROR, "unexpected label: \"{}\", expected \"{}\"", label, expected_label);
    num = line.substr(expected_label.size()+1);
    geom_elem_mapping.resize(std::stoul(num));
    {
        size_t id;
        for(size_t i = 0; i < geom_elem_mapping.size(); ++i){
            std::getline(file, num, ' ');
            id = std::stoul(num);

            geom_elem_mapping[i] = id;
        }
    }

    std::getline(file, line);
    std::getline(file, line);
    std::getline(file, line);
    expected_label = "Elem node tags";
    label = line.substr(0, expected_label.size());
    logger::log_assert(label == expected_label, logger::ERROR, "unexpected label: \"{}\", expected \"{}\"", label, expected_label);
    num = line.substr(expected_label.size()+1);
    elem_node_tags.resize(std::stoul(num));
    {
        size_t id;
        for(size_t i = 0; i < elem_node_tags.size(); ++i){
            std::getline(file, num, ' ');
            id = std::stoul(num);

            elem_node_tags[i] = id;
        }
    }

    std::getline(file, line);
    std::getline(file, line);
    std::getline(file, line);
    expected_label = "Bound elem node tags";
    label = line.substr(0, expected_label.size());
    logger::log_assert(label == expected_label, logger::ERROR, "unexpected label: \"{}\", expected \"{}\"", label, expected_label);
    num = line.substr(expected_label.size()+1);
    bound_elem_node_tags.resize(std::stoul(num));
    {
        size_t id;
        for(size_t i = 0; i < bound_elem_node_tags.size(); ++i){
            std::getline(file, num, ' ');
            id = std::stoul(num);

            bound_elem_node_tags[i] = id;
        }
    }

    std::getline(file, line);
    std::getline(file, line);
    std::getline(file, line);
    expected_label = "Duplicate map";
    label = line.substr(0, expected_label.size());
    logger::log_assert(label == expected_label, logger::ERROR, "unexpected label: \"{}\", expected \"{}\"", label, expected_label);
    num = line.substr(expected_label.size()+1);
    const size_t map_size = std::stoul(num); 
    duplicate_map.reserve(map_size);
    {
        size_t id1, id2;
        for(size_t i = 0; i < map_size; ++i){
            std::getline(file, num, ' ');
            id1 = std::stoul(num);
            std::getline(file, num);
            id2 = std::stoul(num);

            duplicate_map[id1] = id2;
        }
    }

    file.close();

    bool deduplicate = this->geometries.size() > 1;

    this->generate_elements(geom_elem_mapping, 
                            elem_node_tags, 
                            bound_elem_node_tags,
                            id_map,
                            duplicate_map,
                            deduplicate,
                            false);
}


}
