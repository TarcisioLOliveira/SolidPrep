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

#include "logger.hpp"
#include "project_data.hpp"
#include "utils.hpp"
#include "field/trabecular_cutoff.hpp"
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>

namespace field{

TrabecularCutoff::TrabecularCutoff(const projspec::DataMap& data):
    show(data.get_bool("display", true)),
    cutoff(data.get_double("cutoff")),
    bone_thickness(data.get_double("minimum_cortical_thickness"))
{
    if(!data.get_bool("STUB")){
        const auto geom_ids = data.get_array("geometries")->get_int_array();
        this->geoms.reserve(geom_ids.size());
        for(auto i:geom_ids){
            if(i >= 0){
                this->geoms.push_back(data.proj->geometries[i].get());
            }
        }
        Field* f = data.proj->fields[data.get_int("density_field")].get();
        logger::log_assert(f->get_type() == Field::Type::SCALAR, logger::ERROR, "trabecular_cutoff material type requires a scalar field");
        logger::log_assert(f->get_sub_type() == Field::SubType::DOMAIN || f->get_sub_type() == Field::SubType::PROJECTION, logger::ERROR, "field subtype for trabecular cutoff field material must be DOMAIN or PROJECTION");

        this->density = static_cast<ScalarField*>(f);
        this->shell = utils::load_shape(data.get_string("shell"), 1);
    }
}

void TrabecularCutoff::generate() {
    size_t elem_num = 0;
    for(auto& g:this->geoms){
        elem_num += g->mesh.size();
    }
    std::vector<double> map_tmp(elem_num, 0);
    size_t elem_offset = 0;
    for(auto& g:this->geoms){
        #pragma omp parallel for
        for(size_t i = 0; i < g->mesh.size(); ++i){
            const size_t map_pos = i + elem_offset;
            const auto& e = g->mesh[i];
            const double rho = this->density->get(e.get(), e->get_centroid());
            if(rho <= this->cutoff){
                const auto c = e->get_centroid();
                TopoDS_Vertex v = BRepBuilderAPI_MakeVertex(c);
                BRepExtrema_DistShapeShape dss(v, this->shell, Extrema_ExtFlag_MIN);
                double dist = dss.Value();
                if(dist > this->bone_thickness){
                    map_tmp[map_pos] = 1;
                }
            }
        }
        elem_offset += g->mesh.size();
    }
    size_t ei = 0;
    for(auto& g:this->geoms){
        for(auto& e:g->mesh){
            this->trabecular_map[e->id] = map_tmp[ei];
            ++ei;
        }
    }
}

void TrabecularCutoff::initialize_views(Visualization* viz) {
    this->distr = viz->add_view("Cortical/Trabecular", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::DENSITY);
}
void TrabecularCutoff::display_views() const {
    if(this->show){
        std::vector<size_t> gi;
        gi.reserve(geoms.size());
        std::vector<double> vals;
        size_t vals_size = 0;
        for(auto& g:geoms){
            gi.push_back(g->id);
            vals_size += g->mesh.size();
        }
        vals.reserve(vals_size);
        for(auto& g:geoms){
            for(auto& e:g->mesh){
                vals.push_back(this->trabecular_map.at(e->id));
            }
        }
        this->distr->update_view(vals, gi);
    }
}

using namespace projspec;
const bool TrabecularCutoff::reg = Factory<Field>::add(
    [](const DataMap& data){
        return std::make_unique<TrabecularCutoff>(data);
    },
    ObjectRequirements{
        "trabecular_cutoff",
        {
            DataEntry{.name = "display", .type = TYPE_BOOL, .required = false},
            DataEntry{.name = "cutoff", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "density_field", .type = TYPE_INT, .required = true},
            DataEntry{.name = "shell", .type = TYPE_RELATIVE_PATH, .required = true},
            DataEntry{.name = "minimum_cortical_thickness", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "geometries", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 0,
                                 .type = TYPE_INT
                             }
                         ),
                     },
        }
    }
);

}
