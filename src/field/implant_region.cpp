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

#include "project_data.hpp"
#include "project_specification/registry.hpp"
#include "field/implant_region.hpp"

namespace field{

ImplantRegion::ImplantRegion(const projspec::DataMap& data):
    elem_info(data.proj->topopt_element),
    r1(data.get_double("r1")),
    r2(data.get_double("r2")),
    decay_distance(data.get_double("decay_distance")),
    a(data.get_array("coefficients")->get_double_array()),
    min_str(a[0]),
    a_orig(a),
    a_len(a.size()),
    show(data.get_bool("display", true))
{
    
    if(!data.get_bool("STUB")){
        this->center_1 = gp_Pnt(
            data.get_array("center1")->get_double(0),
            data.get_array("center1")->get_double(1),
            data.get_array("center1")->get_double(2)
        );
        this->center_2 = gp_Pnt(
            data.get_array("center2")->get_double(0),
            data.get_array("center2")->get_double(1),
            data.get_array("center2")->get_double(2)
        );
        this->normal = gp_Vec(center_1, center_2);
        this->max_l = center_1.Distance(center_2);

        const auto geom_ids = data.get_array("geometries")->get_int_array();
        this->geoms.reserve(geom_ids.size());
        for(auto i:geom_ids){
            if(i >= 0){
                this->geoms.push_back(data.proj->geometries[i].get());
            }
        }

        Field* f = data.proj->fields[data.get_int("density_field")].get();
        logger::log_assert(f->get_type() == Field::Type::SCALAR, logger::ERROR, "bone material type requires a density field");
        logger::log_assert(f->get_sub_type() == Field::SubType::DOMAIN, logger::ERROR, "field subtype for density field material must be DOMAIN");

        this->density_field = static_cast<ScalarField*>(f);
    }
}

double ImplantRegion::get_implant_multiplier(const gp_Pnt& p) const{
    const gp_Vec v1(this->center_1, p);

    const double l = v1.Dot(this->normal);
    const gp_Vec v2 = v1 - l*this->normal;
    const double dist = v2.Magnitude();
    if(l >= 0 && l <= max_l){
        const double r = (r2 - r1)*l/max_l + r1;
        if(dist < r){
            return min_str;
        } else if(dist > r + decay_distance){
            return 1;
        } else {
            const double x = dist - r;
            return this->f(x);
        }
    } else if(l < 0){
        if(dist <= r1){
            if(-l > decay_distance){
                return 1;
            } else {
                return this->f(-l);
            }
        } else {
            const double dr = dist - r1;
            const double dc = std::sqrt(l*l + dr*dr);
            if(dc > decay_distance){
                return 1;
            } else {
                return this->f(dc);
            }
        }
    } else if(l > max_l){
        const double dl = l - max_l;
        if(dist <= r2){
            if(dl > decay_distance){
                return 1;
            } else {
                return this->f(dl);
            }
        } else {
            const double dr = dist - r2;
            const double dc = std::sqrt(dl*dl + dr*dr);
            if(dc > decay_distance){
                return 1;
            } else {
                return this->f(dc);
            }
        }
    }
    return 1;
}

void ImplantRegion::generate(){}
void ImplantRegion::initialize_views(Visualization* viz){
    if(this->show){
        this->density = viz->add_view("Apparent Density (With Implant)", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::DENSITY);
    }
}
void ImplantRegion::display_views() const{
    if(this->show){
        std::vector<size_t> ids(this->geoms.size(), 0);
        size_t elem_num = 0;
        for(size_t i = 0; i < ids.size(); ++i){
            ids[i] = this->geoms[i]->id;
            elem_num += this->geoms[i]->mesh.size();
        }
        std::vector<double> elem_rho;
        elem_rho.reserve(elem_num);
        for(const auto& g:geoms){
            for(const auto& e:g->mesh){
                const auto c = e->get_centroid();
                const auto m = this->get_implant_multiplier(c);
                elem_rho.push_back(m*this->density_field->get(e.get(), c));
            }
        }
        this->density->update_view(elem_rho, ids);
    }
}
double ImplantRegion::get(const MeshElement* e, const gp_Pnt& p) const{
    if(!this->frozen){
        const auto m = this->get_implant_multiplier(p);
        const auto rho = this->density_field->get(e, p);

        return m*rho;
    } else {
        return this->frozen_values[this->elem_id_to_value.at(e->id)];
    }
}


void ImplantRegion::freeze(std::vector<double>& current_values, std::vector<double>& maximum_values){
    this->frozen = true;

    size_t elem_num = 0;
    for(const auto& g:this->geoms){
        elem_num += g->mesh.size();
    }
    current_values.resize(elem_num);
    maximum_values.resize(elem_num);
    this->frozen_values.resize(elem_num);
    size_t id = 0;
    for(const auto& g:this->geoms){
        for(const auto& e:g->mesh){
            const auto c = e->get_centroid();
            maximum_values[id] = this->density_field->get(e.get(), c);
            const auto m = this->get_implant_multiplier(c);
            current_values[id] = m*maximum_values[id];
            this->frozen_values[id] = current_values[id];
            this->elem_id_to_value[e->id] = id;

            ++id;
        }
    }
}

using namespace projspec;
const bool ImplantRegion::reg = Factory<Field>::add(
    [](const DataMap& data){
        return std::make_unique<ImplantRegion>(data);
    },
    ObjectRequirements{
        "implant_region",
        {
            DataEntry{.name = "density_field", .type = TYPE_INT, .required = true},
            DataEntry{.name = "display", .type = TYPE_BOOL, .required = false},
            DataEntry{.name = "geometries", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 0,
                                 .type = TYPE_INT
                             }
                         ),
                     },
            DataEntry{.name = "center1", .type = TYPE_ARRAY, .required = true,
                 .array_data = std::shared_ptr<ArrayRequirements>(
                     new ArrayRequirements{
                         .size = 3,
                         .type = TYPE_DOUBLE
                     }
                 ),
            },
            DataEntry{.name = "center2", .type = TYPE_ARRAY, .required = true,
                 .array_data = std::shared_ptr<ArrayRequirements>(
                     new ArrayRequirements{
                         .size = 3,
                         .type = TYPE_DOUBLE
                     }
                 ),
            },
            DataEntry{.name = "r1", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "r2", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "decay_distance", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "coefficients", .type = TYPE_ARRAY, .required = true,
                 .array_data = std::shared_ptr<ArrayRequirements>(
                     new ArrayRequirements{
                         .size = 0,
                         .type = TYPE_DOUBLE
                     }
                 ),
            },
        }
    }
);
}
