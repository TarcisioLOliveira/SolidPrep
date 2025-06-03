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

#include "simulation/marginal_bone_loss.hpp"
#include "logger.hpp"
#include "material/mandible.hpp"
#include "math/matrix.hpp"
#include <algorithm>
#include "project_data.hpp"

namespace simulation{

MarginalBoneLoss::MarginalBoneLoss(const projspec::DataMap& data):
    mesh(data.proj->topopt_mesher.get()),
    fem(data.proj->topopt_fea.get()),
    pc(data.get_double("pc")),
    t(data.get_array("traction")->get_double_array_fixed<RANGE_NUM>()),
    c(data.get_array("compression")->get_double_array_fixed<RANGE_NUM>()),
    s(data.get_array("shear")->get_double_array_fixed<RANGE_NUM>()),
    K_e1(this->get_K_e1(t, c)),
    K_g(this->get_K_g(s)),
    K_e2(this->get_K_e2(t, c)),
    problem_type(data.proj->type),
    //a0(data.get_array("a0")->get_double_array()),
    //a1(data.get_array("a1")->get_double_array()),
    //a2(data.get_array("a2")->get_double_array()),
    //a3(data.get_array("a3")->get_double_array()),
    eps_to_x(this->make_eps_to_x()),
    geom_id(data.get_int("mandible_geometry")),
    mandible(data.proj->geometries[geom_id].get()),
    time_step(data.get_double("time_step")),
    time_limit(data.get_double("time_limit")),
    dv1_orig(data.get_double("volume_variation_rate")),
    dv2_orig(data.get_double("overload_volume_variation_rate"))
{
    logger::log_assert(mandible->materials.get_materials()[0]->get_type() == Material::MANDIBLE,
            logger::ERROR,
            "geometry selected for marginal bone loss simulation must have a material of type \"mandible\"");

    this->dv0 = dv1_orig*(this->c[0] - 0.0);
    this->dv1 = 0;
    this->dv2 = dv1_orig*(this->c[2] - this->c[1]);
    this->dv3 = -dv2_orig*(this->c[3] - this->c[2]);

    StrainVector3D eps(6);

    eps[0] = -this->c[1];
    this->lhs2_offset = this->LHS_3D(2, eps);
    eps[0] = -this->c[2];
    this->lhs3_offset = this->LHS_3D(3, eps);

    this->max_dv2 = this->dv2*(1 - this->lhs2_offset);
}


MarginalBoneLoss::Range MarginalBoneLoss::get_K_e1(const Range& t, const Range& c) const{
    Range K;
    for(size_t i = 0; i < this->RANGE_NUM; ++i){
        const double Ki_tmp = (t[i]+c[i])/(2*t[i]*c[i]);
        K[i] = 0.5*Ki_tmp*Ki_tmp;
    }
    return K;
}
MarginalBoneLoss::Range MarginalBoneLoss::get_K_e2(const Range& t, const Range& c) const{
    Range K;
    for(size_t i = 0; i < this->RANGE_NUM; ++i){
        K[i] = -(t[i]-c[i])/(2*t[i]*c[i]);
    }
    return K;
}
MarginalBoneLoss::Range MarginalBoneLoss::get_K_g (const Range& s) const{
    Range K;
    for(size_t i = 0; i < this->RANGE_NUM; ++i){
        K[i] = 1.0/(2*s[i]*s[i]);
    }
    return K;
}

void MarginalBoneLoss::initialize_views(Visualization* viz){
    this->stress_view = viz->add_view("Von Mises Stress", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
    this->density_view = viz->add_view("Elemental Density", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::OTHER);
    this->growth_view = viz->add_view("Growth (Element Density)", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::OTHER);
}
void MarginalBoneLoss::initialize(){

}
void MarginalBoneLoss::run(){
    size_t elem_num = 0;
    size_t mand_geom_offset = 0;
    size_t density_num = 0;
    std::vector<double> v;
    {
        size_t id = 0;
        for(const auto& g:mesh->geometries){
            if(id == this->geom_id){
                mand_geom_offset = elem_num;
                density_num = g->mesh.size();
                v.resize(density_num);
                for(size_t i = 0; i < density_num; ++i){
                    v[i] = g->mesh[i]->get_volume(1);
                }
            }
            elem_num += g->mesh.size();
            ++id;
        }
    }
    std::vector<double> density(density_num, 1);
    std::vector<double> density_vector_view(elem_num, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> stress(elem_num, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> growth(elem_num, std::numeric_limits<double>::quiet_NaN());

    std::vector<double> u(this->mesh->max_dofs);

    //material::Mandible* material = static_cast<material::Mandible*>(this->mandible->materials.get_materials()[0].get());

    // TODO: make this affect CT scan based material
    double time_count = 0;
    while(time_count < this->time_limit){

        logger::quick_log("");
        logger::quick_log("Time:", time_count);

        //material->set_maturation_alpha(alpha_count);
        this->fem->generate_matrix(this->mesh, density, this->pc);
        this->fem->calculate_displacements_global(this->mesh, this->mesh->load_vector, u);

        std::copy(density.begin(), density.end(), density_vector_view.begin() + mand_geom_offset);
        this->density_view->update_view(density_vector_view);

        {
            auto s_it = stress.begin();
            size_t D_offset = 0;
            for(const auto& g:this->mesh->geometries){
                g->get_stresses(u, false, D_offset, this->fem->D_vec, s_it);
            }
            for(size_t i = 0; i < density_num; ++i){
                stress[i + mand_geom_offset] *= std::pow(density[i], this->pc);
            }
            this->stress_view->update_view(stress);
        }

        for(size_t i = 0; i < density_num; ++i){
            const auto& e = this->mandible->mesh[i];
            const auto strain = 1e6*e->get_strain_vector(e->get_centroid(), u);
            const double gi = this->get_density_variation(strain, v[i])*this->time_step;
            growth[mand_geom_offset + i] = gi;
            density[i] += gi;
            if(density[i] < 0){
                density[i] = 0;
            } else if(density[i] > 1){
                density[i] = 1;
            }
        }
        this->growth_view->update_view(growth);

        time_count += this->time_step;
    }
}

math::Vector MarginalBoneLoss::make_eps_to_x() const{
    math::Matrix M(RANGE_NUM, RANGE_NUM);
    math::Vector b(RANGE_NUM);

    StrainVector3D eps(6);

    for(size_t i = 0; i < RANGE_NUM; ++i){
        b[i] = i + 1;
        eps[0] = -this->c[i];
        for(size_t j = 0; j < RANGE_NUM; ++j){
            const double LHS = this->LHS_3D(j, eps);
            M(i,j) = LHS;
        } 
    }

    math::LU LU(M);
    LU.solve(b);

    return b;
}

double MarginalBoneLoss::get_density_variation(const StrainVector3D& eps, const double v) const{
    StrainVector3D lhs(RANGE_NUM);
    for(size_t j = 0; j < RANGE_NUM; ++j){
        lhs[j] = this->LHS_3D(j, eps);
    }

    if(lhs[0] < 1){
        return dv0*(1 - lhs[0])/v;
    } else if(lhs[1] < 1){
        return 0;
    } else if(lhs[2] < 1){
        return dv2*(lhs[2] - this->lhs2_offset)/v;
    } else {
        return max_dv2 + dv3*(lhs[3] - this->lhs3_offset)/v;
    }
}

using namespace projspec;
const bool MarginalBoneLoss::reg = Factory<Simulation>::add(
    [](const DataMap& data){
        return std::make_unique<MarginalBoneLoss>(data);
    },
    ObjectRequirements{
        "marginal_bone_loss",
        {
            DataEntry{.name = "time_step", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "volume_variation_rate", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "overload_volume_variation_rate", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "time_limit", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "pc", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "mandible_geometry", .type = TYPE_INT, .required = true},
            DataEntry{.name = "traction", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = MarginalBoneLoss::RANGE_NUM,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
            DataEntry{.name = "compression", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = MarginalBoneLoss::RANGE_NUM,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
            DataEntry{.name = "shear", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = MarginalBoneLoss::RANGE_NUM,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
        }
    }
);

}
