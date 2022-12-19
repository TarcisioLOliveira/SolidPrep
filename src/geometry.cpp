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

#include "geometry.hpp"

Geometry::Geometry(const std::string& path, double scale, utils::ProblemType type,
        MeshElementFactory* elem_type, bool do_topopt, bool with_void, Material* material, std::vector<Material*> alt_materials ):
    shape(utils::load_shape(path, scale)),
    material(material), alternative_materials(alt_materials.begin(), alt_materials.end()),
    element_type(elem_type), do_topopt(do_topopt), with_void(with_void),
    constitutive_matrices(this->init_constitutive_matrices(type)), mesh(), 
    type(type){}

Geometry::Geometry(TopoDS_Shape shape, utils::ProblemType type,
        MeshElementFactory* elem_type, bool do_topopt, bool with_void, Material* material, std::vector<Material*> alt_materials ):
    shape(std::move(shape)),
    material(material), alternative_materials(alt_materials.begin(), alt_materials.end()),
    element_type(elem_type), do_topopt(do_topopt), with_void(with_void),
    constitutive_matrices(this->init_constitutive_matrices(type)), mesh(), 
    type(type){}


