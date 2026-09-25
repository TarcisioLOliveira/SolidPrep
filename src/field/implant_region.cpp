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
#include "math/matrix.hpp"
#include "project_data.hpp"
#include "project_specification/registry.hpp"
#include "field/implant_region.hpp"

namespace field{

ImplantRegion::ImplantRegion(const projspec::DataMap& data):
    elem_info(data.proj->topopt_element),
    r1(data.get_double("r1")),
    r2(data.get_double("r2")),
    r3(data.get_double("r3")),
    decay_distance(data.get_double("decay_distance")),
    a(data.get_array("coefficients")->get_double_array()),
    min_str(a[0]),
    a_orig(a),
    a_len(a.size()),
    NODES_PER_ELEM(elem_info->get_nodes_per_element()),
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
        this->center_3 = gp_Pnt(
            data.get_array("center3")->get_double(0),
            data.get_array("center3")->get_double(1),
            data.get_array("center3")->get_double(2)
        );
        this->normal = gp_Vec(center_1, center_3);
        this->max_l = center_1.Distance(center_3);

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

        this->is_shape_opt = data.proj->do_shape_opt;

        // Create radius interpolation function
        const std::array<double, 3> l{
            0,
            gp_Vec(this->center_1, this->center_2).Dot(this->normal),
            gp_Vec(this->center_1, this->center_3).Dot(this->normal)
        };
        math::Matrix A(3, 3);
        for(size_t i = 0; i < 3; ++i){
            for(size_t j = 0; j < 3; ++j){
                A(i,j) = std::pow(l[i], j);
            }
        }
        math::Vector r({this->r1, this->r2, this->r3});
        math::LU lu(A);
        lu.solve(r);
        for(size_t i = 0; i < 3; ++i){
            this->p[i] = r[i];
        }
    }
}

double ImplantRegion::get_implant_multiplier(const gp_Pnt& p) const{
    const gp_Vec v1(this->center_1, p);

    const double l = v1.Dot(this->normal);
    const gp_Vec v2 = v1 - l*this->normal;
    const double dist = v2.Magnitude();
    if(l >= 0 && l <= max_l){
        const double r = this->r(l);
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
        if(dist <= r3){
            if(dl > decay_distance){
                return 1;
            } else {
                return this->f(dl);
            }
        } else {
            const double dr = dist - r3;
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

void ImplantRegion::generate(){
    if(this->is_shape_opt){
        std::vector<double> tmp;
        this->freeze(tmp);
    }
}
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
        const auto Ni = e->get_Ni_1dof(p);
        double v = 0;
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            const auto n = this->id_pos_map.at(e->nodes[i]->id);
            v += Ni[i]*this->frozen_values[n];
        }
        return v;
    }
}
double ImplantRegion::get(const Node* n) const{
    if(!this->frozen){
        const auto m = this->get_implant_multiplier(n->point);
        const auto rho = this->density_field->get(n);

        return m*rho;
    } else {
        const auto ni = this->id_pos_map.at(n->id);
        return this->frozen_values[ni];
    }
}
double ImplantRegion::get_dsh(const MeshElement* e, const gp_Pnt& p, const size_t n, const size_t dof) const{
    if(this->frozen){
        const auto dNi = e->get_dN_1dof_sh(p, n, dof);
        double v = 0;
        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
            const auto n = this->id_pos_map.at(e->nodes[i]->id);
            v += dNi[i]*this->frozen_values[n];
        }
        return v;
    }
    return 0;
}
std::vector<double> ImplantRegion::get_dsh_full(const MeshElement* e, const gp_Pnt& p, const size_t n, const size_t dof) const{
    std::vector<double> dri(n*dof, 0);
    if(this->frozen){
        std::vector<size_t> npos(n, 0);
        for(size_t i = 0; i < n; ++i){
            npos[i] = this->id_pos_map.at(e->nodes[i]->id);
        }
        for(size_t i = 0; i < n; ++i){
            for(size_t j = 0; j < dof; ++j){
                const auto dNi = e->get_dN_1dof_sh(p, i, j);
                for(size_t k = 0; k < n; ++k){
                    dri[i*dof + j] += dNi[k]*this->frozen_values[npos[k]];
                }
            }
        }
    }
    return dri;
}


void ImplantRegion::freeze(std::vector<double>& current_values){
    for(const auto& g:this->geoms){
        for(const auto& n:g->node_list){
            this->id_pos_map[n->id] = 0;
        }
    }

    size_t new_id = 1;
    for(const auto& g:this->geoms){
        for(const auto& n:g->node_list){
            auto& idi = this->id_pos_map[n->id];
            if(idi == 0){
                idi = new_id;
                ++new_id;
            }
        }
    }
    for(auto& i:id_pos_map){
        --i.second;
    }
    --new_id;
    this->frozen_values.resize(new_id, 0);


    //std::map<const Node*, const MeshElement*> node_to_elem;
    //for(auto& g:this->geoms){
    //    for(auto& e:g->mesh){
    //        for(size_t i = 0; i < NODES_PER_ELEM; ++i){
    //            node_to_elem[e->nodes[i]] = e.get();
    //        }
    //    }
    //}

    current_values.resize(new_id, 0);
    this->frozen_values.resize(new_id, 0);
    //for(const auto& p:node_to_elem){
    //    const auto id = this->id_pos_map[p.first->id];
    //    const auto val = this->get(p.second, p.first->point);
    //    this->frozen_values[id] = val;
    //    current_values[id] = val;
    //}
    size_t id = 0;
    for(auto& g:this->geoms){
        for(auto& n:g->node_list){
            const auto val = this->get(n);
            this->frozen_values[id] = val;
            current_values[id] = val;

            ++id;
        }
    }

    this->frozen = true;
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
            DataEntry{.name = "center3", .type = TYPE_ARRAY, .required = true,
                 .array_data = std::shared_ptr<ArrayRequirements>(
                     new ArrayRequirements{
                         .size = 3,
                         .type = TYPE_DOUBLE
                     }
                 ),
            },
            DataEntry{.name = "r1", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "r2", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "r3", .type = TYPE_DOUBLE, .required = true},
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
