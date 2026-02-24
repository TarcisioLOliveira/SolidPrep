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
#include <BRepBuilderAPI_MakeEdge.hxx>

namespace field{

TrabecularCutoff::BoundaryApproximation::BoundaryApproximation(const projspec::DataMap* const data, const TopoDS_Shape& shell):
    defined(true),
    center1(data->get_array("center1")->get_double_array()),    
    center2(data->get_array("center2")->get_double_array()),
    long_dir(data->get_array("longitudinal_direction")->get_double_array()),
    v_dir(data->get_array("vertical_direction")->get_double_array()),
    h_dir(data->get_array("horizontal_direction")->get_double_array()),
    th_top(data->get_double("th_top")),
    th_bot(data->get_double("th_bot")),
    th_left(data->get_double("th_left")),
    th_right(data->get_double("th_right")),
    pts_long(data->get_int("points_long")),
    pts_v(data->get_int("points_vert"))
{
    this->long_dir.normalize();
    this->v_dir.normalize();
    this->h_dir.normalize();

    // Reposition center2 as the projection of center1 on the other face
    this->center2 = (center2 - center1).dot(long_dir)*long_dir + center1;

    this->bounds.resize(pts_long);
    this->base_points.resize(pts_long);

    const double long_l = (center2 - center1).dot(long_dir);
    this->spacing_l = long_l/(pts_long - 1);

    for(size_t i = 0; i < pts_long; ++i){
        // Generate vertical lines
        const auto cur_pnt_vec = center1 + (i*this->spacing_l)*long_dir;
        const auto cur_pnt = to_occt_point(cur_pnt_vec);

        bounds[i].resize(pts_v);

        gp_Lin lin(cur_pnt, to_occt_dir(v_dir));
        auto edge = BRepBuilderAPI_MakeEdge(lin);
        BRepExtrema_DistShapeShape dss;
        dss.SetMultiThread(true);
        dss.LoadS1(edge);
        dss.LoadS2(shell);
        dss.Perform();
        double dist = dss.Value();
        logger::log_assert(dss.NbSolution() == 2, logger::ERROR,
                            "each line must have exactly 2 intersections, found {}", dss.NbSolution());
        logger::log_assert(dist < Precision::Confusion(), logger::ERROR,
                            "especified vertical direction does not always intersect shell");

        const auto inter1 = to_vector(dss.PointOnShape1(1));
        const auto inter2 = to_vector(dss.PointOnShape1(2));

        const double min_pos = std::min((inter1 - cur_pnt_vec).dot(v_dir), (inter2 - cur_pnt_vec).dot(v_dir));

        this->base_points[i].p1 = (min_pos + th_bot)*v_dir + cur_pnt_vec;
        this->base_points[i].l = std::abs((inter1 - inter2).dot(v_dir)) - (th_bot + th_top);
        this->base_points[i].spacing_v = this->base_points[i].l/(pts_v - 1);

        // Generate horizontal lines
        const auto& cur_base = this->base_points[i];
        for(size_t j = 0; j < pts_v; ++j){
            const auto cur_pnt_vec_h = cur_base.p1 + (j*cur_base.spacing_v)*v_dir;
            const auto cur_pnt_h = to_occt_point(cur_pnt_vec_h);

            gp_Lin lin_h(cur_pnt_h, to_occt_dir(h_dir));
            auto edge_h = BRepBuilderAPI_MakeEdge(lin_h);
            BRepExtrema_DistShapeShape dss_h;
            dss_h.SetMultiThread(true);
            dss_h.LoadS1(edge_h);
            dss_h.LoadS2(shell);
            dss_h.Perform();
            double dist_h = dss_h.Value();
            logger::log_assert(dss_h.NbSolution() == 2, logger::ERROR,
                                "each horizontal line must have exactly 2 intersections, found {}", dss_h.NbSolution());
            logger::log_assert(dist_h < Precision::Confusion(), logger::ERROR,
                                "especified horizontal direction does not always intersect shell");

            const auto inter1_h = to_vector(dss_h.PointOnShape1(1));
            const auto inter2_h = to_vector(dss_h.PointOnShape1(2));

            const double d1 = (inter1_h - cur_pnt_vec_h).dot(h_dir);
            const double d2 = (inter2_h - cur_pnt_vec_h).dot(h_dir);
            bounds[i][j].d1 = std::max(d1, d2) - this->th_right;
            bounds[i][j].d2 = std::min(d1, d2) + this->th_left;
        }
    }
}

bool TrabecularCutoff::BoundaryApproximation::is_inside(const gp_Pnt& p) const{
    if(this->defined){
        const auto pv = to_vector(p);
        const auto relp = pv - center1;
        const auto relpos = relp.dot(long_dir)/spacing_l;
        const size_t id1 = std::floor(relpos);
        const size_t id2 = id1 + 1;
        const double w_l = relpos - id1;

        const auto relp2_1 = pv - this->base_points[id1].p1;
        const auto relp2_2 = pv - this->base_points[id2].p1;
        const auto relpos2 = ((1.0 - w_l)*relp2_1 + w_l*relp2_2).dot(v_dir);
        const double weigthed_height = (1.0 - w_l)*this->base_points[id1].l
                                        + w_l*this->base_points[id2].l;
        if(relpos2 < 0 || relpos2 > weigthed_height){
            return false;
        }

        const auto rp1_v = relp2_1.dot(v_dir)/this->base_points[id1].spacing_v;
        const auto rp2_v = relp2_2.dot(v_dir)/this->base_points[id2].spacing_v;
        const size_t id1_v = std::floor(rp1_v);
        const size_t id2_v = std::floor(rp2_v);
        const double w1_v = rp1_v - id1_v;
        const double w2_v = rp2_v - id2_v;

        const double w_d1 = (1.0 - w_l)*((1.0 - w1_v)*bounds[id1][id1_v].d1 + w1_v*bounds[id1][id1_v + 1].d1) +
                                    w_l*((1.0 - w2_v)*bounds[id2][id2_v].d1 + w2_v*bounds[id2][id2_v + 1].d1);
        const double w_d2 = (1.0 - w_l)*((1.0 - w1_v)*bounds[id1][id1_v].d2 + w1_v*bounds[id1][id1_v + 1].d2) +
                                    w_l*((1.0 - w2_v)*bounds[id2][id2_v].d2 + w2_v*bounds[id2][id2_v + 1].d2);

        const auto ref_v = (1.0 - w_l)*(relp2_1.dot(v_dir)*v_dir + this->base_points[id1].p1) +
                                   w_l*(relp2_2.dot(v_dir)*v_dir + this->base_points[id2].p1);

        const auto relpos_v = (pv - ref_v).dot(h_dir);

        if(relpos_v < w_d2 || relpos_v > w_d1){
            return false;
        }
    }
    return true;
}


TrabecularCutoff::TrabecularCutoff(const projspec::DataMap& data):
    show(data.get_bool("display", true)),
    cutoff(data.get_double("cutoff"))
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
        if(data.exists_object("inner_bounds")){
            this->inner_bounds = BoundaryApproximation(data.get_object("inner_bounds"), this->shell);
        }
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
        for(size_t i = 0; i < g->mesh.size(); ++i){
            const size_t map_pos = i + elem_offset;
            const auto& e = g->mesh[i];
            const auto p = e->get_centroid();
            const double rho = this->density->get(e.get(), p);
            if(rho <= this->cutoff && this->inner_bounds.is_inside(p)){
                map_tmp[map_pos] = 1;
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
    if(this->show){
        this->distr = viz->add_view("Cortical/Trabecular", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::DENSITY);
    }
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
            DataEntry{.name = "geometries", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 0,
                                 .type = TYPE_INT
                             }
                         ),
                     },
            DataEntry{.name = "inner_bounds", .type = TYPE_OBJECT, .required = false,
                .object_data = {
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
                    DataEntry{.name = "longitudinal_direction", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 3,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                    },
                    DataEntry{.name = "vertical_direction", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 3,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                    },
                    DataEntry{.name = "horizontal_direction", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 3,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                    },
                    DataEntry{.name = "th_top", .type = TYPE_DOUBLE, .required = true},
                    DataEntry{.name = "th_bot", .type = TYPE_DOUBLE, .required = true},
                    DataEntry{.name = "th_left", .type = TYPE_DOUBLE, .required = true},
                    DataEntry{.name = "th_right", .type = TYPE_DOUBLE, .required = true},
                    DataEntry{.name = "points_long", .type = TYPE_INT, .required = true},
                    DataEntry{.name = "points_vert", .type = TYPE_INT, .required = true},
                }
            }
        }
    }
);

}
