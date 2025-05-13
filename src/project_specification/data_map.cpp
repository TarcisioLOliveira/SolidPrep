/*
 *   Copyright (C) 2025 Tarcísio Ladeia de Oliveira.
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

#include "project_specification/data_map.hpp"

namespace projspec{

DataMap* DataArray::get_object(size_t key, DataMap* none) const{
    static DataMap empty_map(nullptr);

    if(none == nullptr){
        none = &empty_map;
    }
    auto found = this->object_map.find(key);
    if(found != this->object_map.end()){
        return found->second.get();
    } else {
        return none;
    }
}

}
