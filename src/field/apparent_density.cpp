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

#include "field/apparent_density.hpp"
#include "logger.hpp"
#include "math/matrix.hpp"
#include "project_data.hpp"
#include <algorithm>
#include <defs.hpp>

namespace field{

NRRDReader::NRRDReader():
    file(nrrdNew(), [](Nrrd* n){nrrdNuke(n);}){

}

void NRRDReader::load(const std::string& file_path){

    if(this->file != nullptr && this->file->data != NULL){
        this->file.reset(nrrdNew());
    }
    {
        auto err = nrrdLoad(this->file.get(), file_path.c_str(), NULL);
        logger::log_assert(err == 0, 
                logger::ERROR, "Error loading NRRD file (return {}): {}", std::string(biffGetDone(NRRD)), err);
    }

    logger::log_assert(file->dim == file->spaceDim, logger::ERROR, "NRRD dim != spaceDim");
    this->DIM = file->spaceDim;

    this->origin = math::Vector(DIM);
    this->R = math::Matrix(DIM, DIM);
    this->len.resize(DIM);

    {
        std::vector<double> origin_tmp(NRRD_SPACE_DIM_MAX);
        nrrdSpaceOriginGet(this->file.get(), origin_tmp.data());
        for(size_t i = 0; i < DIM; ++i){
            this->origin[i] = origin_tmp[i];
        }
    }
    this->lup = nrrdDLookup[this->file->type];

    for(size_t i = 0; i < file->dim; ++i){
        this->len[i] = file->axis[i].size;
        for(size_t j = 0; j < file->spaceDim; ++j){
            this->R(j, i) = file->axis[i].spaceDirection[j];
        }
    }
    this->Rinv = this->R.get_inverted_LU();

}

double NRRDReader::get_data(const std::vector<size_t>& p) const{
    size_t pos = 0;
    for(size_t i = 0; i < this->DIM; ++i){
        const size_t j = (this->DIM - 1) - i;
        pos *= this->len[j];
        pos += p[j];
    }

    return this->lup(this->file->data, pos);
}

double NRRDReader::get_averaged(const math::Vector& point) const{
    
    const auto dist = [](const math::Vector& v1, const math::Vector& v2)->double{
        double d = 0;
        for(size_t i = 0; i < std::min(v1.get_N(), v2.get_N()); ++i){
            const double dvi = v1[i] - v2[i];
            d += dvi*dvi;
        }

        return std::sqrt(d);
    };

    math::Vector ijk = this->Rinv*(point - this->origin);
    bool within = true;
    for(size_t i = 0; i < this->DIM; ++i){
        if(ijk[i] < 0 || ijk[i] > this->len[i] - 1){
            within = false;
            break;
        }
    }
    if(!within){
        logger::quick_log("NRRD Origin:");
        logger::quick_log(this->origin);
        logger::quick_log("NRRD R:");
        logger::quick_log(this->R);
        logger::quick_log("NRRD sizes:");
        logger::quick_log(this->len);
        logger::quick_log("NRRD ijk:");
        logger::quick_log(ijk);
        logger::log_assert(within, logger::ERROR, "position out of bounds in NRRD: {} {} {}", point[0], point[1], point[2]);
    }
    std::vector<size_t> ijk_uint(this->DIM);

    for(size_t i = 0; i < this->DIM; ++i){
        ijk_uint[i] = ijk[i];
    }
    std::vector<double> w(std::pow(2, this->DIM));
    double w_sum = 0;
    math::Vector pi(DIM);
    std::vector<size_t> pui(DIM);
    double avg = 0;
    size_t it_max = std::pow(2, DIM);
    for(size_t it = 0; it < it_max; ++it){
        size_t i = it % 2;
        size_t j = (it / 2) % 2;
        size_t k = (it / 4) % 2;
        pui[0] = ijk_uint[0] + i;
        pui[1] = ijk_uint[1] + j;
        pui[2] = ijk_uint[2] + k;
        pi[0] = pui[0];
        pi[1] = pui[1];
        pi[2] = pui[2];
        const double wi = dist(pi, ijk);
        w_sum += wi;
        avg += wi*this->get_data(pui); 
    }

    return avg/w_sum;
}

ApparentDensity::ApparentDensity(const projspec::DataMap& data):
    elem_info(data.proj->topopt_element),
    thickness(data.proj->thickness),
    R(data.get_matrix("R")),
    t(data.get_array("t")->get_double_array()),
    radius(data.get_double("smoothing_radius")),
    rho_func(data.get_array("HU_to_density")->get_double_array_fixed<2>()),
    show(data.get_bool("display", true)),
    DIM((elem_info->get_problem_type() == utils::PROBLEM_TYPE_2D) ? 2 : 3),
    NODES_PER_ELEM(elem_info->get_nodes_per_element())
{
    if(!data.get_bool("STUB")){
        const auto geom_ids = data.get_array("geometries")->get_int_array();
        this->geoms.reserve(geom_ids.size());
        for(auto i:geom_ids){
            if(i >= 0){
                this->geoms.push_back(data.proj->geometries[i].get());
            }
        }
        this->nrrd.load(data.get_string("nrrd_path"));
    }
}

void ApparentDensity::generate(){
    std::set<size_t> node_ids;
    for(auto& g:this->geoms){
        for(auto& n:g->node_list){
            node_ids.insert(n->id);
        }
    }
    size_t new_id = 0;
    for(auto& n:node_ids){
        id_pos_map[n] = new_id;
        ++new_id;
    }
    this->nodal_densities.resize(new_id, 0);
    // // Test run
    // for(auto& g:this->geoms){
    //     for(auto& n:g->node_list){
    //         const size_t ni = this->id_pos_map[n->id];
    //         const auto pi = n->point;
    //         const math::Vector pv({pi.X(), pi.Y(), pi.Z()});
    //         const math::Vector Rp(R.T()*(pv - t));
    //         nodal_densities[ni] = this->rho_func[0]*this->nrrd.get_averaged(Rp) + this->rho_func[1];
    //     }
    // }
    
    // Generate problem
    const math::Matrix I({1, 0, 0,
                          0, 1, 0,
                          0, 0, 1}, 3, 3);
    general_solver::MUMPSGeneral solver;
    solver.initialize_matrix(true, new_id);
    std::vector<long> pos(NODES_PER_ELEM);
    double rho_o_min = 1e100, rho_o_max = -1e100;
    for(auto& g:this->geoms){
        for(auto& e:g->mesh){
            const auto pi = e->get_centroid();
            const math::Vector pv({pi.X(), pi.Y(), pi.Z()});
            const math::Vector Rp(R.T()*(pv - t));
            const double rho_i = (this->rho_func[0]*this->nrrd.get_averaged(Rp) + this->rho_func[1]);
            rho_o_min = std::min(rho_o_min, rho_i);
            rho_o_max = std::max(rho_o_max, rho_i);
            const math::Vector Ni = rho_i*e->source_1dof(this->thickness);

            // Source vector
            for(size_t i = 0; i < NODES_PER_ELEM; ++i){
                const size_t ni = this->id_pos_map[e->nodes[i]->id];
                nodal_densities[ni] += Ni[i];
                pos[i] = ni;
            }

            const auto M = radius*radius*e->diffusion_1dof(thickness, I) +
                           e->absorption_1dof(thickness);
            solver.add_element(M, pos);
        }
    }
    solver.compute();
    solver.solve(this->nodal_densities);

    math::Matrix A(5, 5, 0);
    math::Vector b(5, 0);

    const double rho_n_min = *std::min_element(this->nodal_densities.cbegin(), this->nodal_densities.cend());
    const double rho_n_max = *std::max_element(this->nodal_densities.cbegin(), this->nodal_densities.cend());

    b[0] = rho_o_min;
    b[1] = rho_o_max;
    for(size_t j = 0; j < 5; ++j){
        A(0, j) = std::pow(rho_n_min, j);
        A(1, j) = std::pow(rho_n_max, j);
    }
    for(size_t j = 1; j < 5; ++j){
        A(2, j) = j*std::pow(rho_n_min, j-1);
        A(3, j) = j*std::pow(rho_n_max, j-1);
    }
    math::Vector rho_iv(NODES_PER_ELEM);
    for(auto& g:this->geoms){
        for(auto& e:g->mesh){
            const auto pi = e->get_centroid();
            for(size_t i = 0; i < NODES_PER_ELEM; ++i){
                rho_iv[i] = this->nodal_densities[id_pos_map[e->nodes[i]->id]];
            }
            const double vi = e->get_volume(this->thickness);
            const math::Vector Nc = e->get_Ni_1dof(pi);
            const double rho_i = Nc.T()*rho_iv;
            const math::Vector Ni = e->source_1dof(this->thickness);

            b[4] += Ni.T()*rho_iv;
            for(size_t j = 0; j < 5; ++j){
                A(4, j) += vi*std::pow(rho_i, j);
            }
        }
    }

    math::LU Mlu(A);
    Mlu.solve(b);

    this->scaler = std::move(b);
}
void ApparentDensity::initialize_views(Visualization* viz){
    this->density = viz->add_view("Apparent Density", spview::defs::ViewType::NODAL, spview::defs::DataType::DENSITY);
}
void ApparentDensity::display_views() const{
    if(this->show){
        std::vector<size_t> gi;
        gi.reserve(geoms.size());
        for(auto& g:geoms){
            gi.push_back(g->id);
        }
        this->density->update_view(this->nodal_densities, gi);
    }
}
double ApparentDensity::get(const MeshElement* e, const gp_Pnt& p) const{
    const auto Ni = e->get_Ni_1dof(p);
    double v = 0;
    for(size_t i = 0; i < NODES_PER_ELEM; ++i){
        const auto n = this->id_pos_map.at(e->nodes[i]->id);
        v += Ni[i]*this->nodal_densities[n];
    }
    double vr = 0;
    for(size_t j = 0; j < 5; ++j){
        vr += scaler[j]*std::pow(v, j);
    }

    return vr;
}

using namespace projspec;
const bool ApparentDensity::reg = Factory<Field>::add(
    [](const DataMap& data){
        return std::make_unique<ApparentDensity>(data);
    },
    ObjectRequirements{
        "apparent_density",
        {
            DataEntry{.name = "nrrd_path", .type = TYPE_RELATIVE_PATH, .required = true},
            DataEntry{.name = "display", .type = TYPE_BOOL, .required = false},
            DataEntry{.name = "smoothing_radius", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "geometries", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 0,
                                 .type = TYPE_INT
                             }
                         ),
                     },
            DataEntry{.name = "R", .type = TYPE_MATRIX, .required = true,
                         .matrix_data = std::shared_ptr<MatrixRequirements>(
                             new MatrixRequirements{
                                 .W = 3,
                                 .H = 3
                             }
                         ),
                     },
            DataEntry{.name = "t", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 3,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
            DataEntry{.name = "HU_to_density", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = 2,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
        }
    }
);

}
