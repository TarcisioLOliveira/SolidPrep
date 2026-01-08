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
#include "field/implant_region.hpp"
#include "field/principal_stress.hpp"
#include "logger.hpp"
#include "math/matrix.hpp"
#include <algorithm>
#include "project_data.hpp"

namespace simulation{

MarginalBoneLoss::MarginalBoneLoss(const projspec::DataMap& data):
    mesh(data.proj->topopt_mesher.get()),
    fem(data.proj->topopt_fea.get()),
    proj(data.proj),
    base_range(data.get_array("intervals")->get_double_array_fixed<RANGE_NUM>()),
    t(this->interval_mult(data.get_double("traction"), base_range)),
    c(this->interval_mult(data.get_double("compression"), base_range)),
    s(this->interval_mult(data.get_double("shear"), base_range)),
    K_e1(this->get_K_e1(t, c)),
    K_g(this->get_K_g(s)),
    K_e2(this->get_K_e2(t, c)),
    problem_type(data.proj->type),
    geom_id(data.get_int("mandible_geometry")),
    mandible(data.proj->geometries[geom_id].get()),
    time_step(data.get_double("time_step")),
    time_limit(data.get_double("time_limit"))
{
    logger::log_assert(
            mandible->materials.get_materials()[0]->get_type() == Material::POWER_LAW_ORTHOTROPIC ||
            mandible->materials.get_materials()[0]->get_type() == Material::POWER_LAW_ISOTROPIC ||
            mandible->materials.get_materials()[0]->get_type() == Material::BONE,
            logger::ERROR,
            "geometry selected for marginal bone loss simulation must have a material of type \"bone\"");

    //this->dv0 = dv1_orig*this->c[3];

    StrainVector3D eps(6);

    this->coeff = math::Vector(3);
    this->coeff[2] = rho_d;

    eps[0] = -this->c[0];
    this->x0 = this->LHS_3D(3, eps);
    eps[0] = -this->c[1];
    this->x1 = this->LHS_3D(3, eps);
    eps[0] = -this->c[2];
    this->x2 = this->LHS_3D(3, eps);

    const double xm = (x1 + x2)/2.0;
    math::Matrix A(
        {1, x1, x1*x1,
         1, x2, x2*x2,
         1, xm, xm*xm},
        3, 3);

    math::LU LU(A);
    LU.solve(this->coeff);

    this->dv0 = coeff[1] + 2*coeff[2]*x1;
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
    this->strain_view = viz->add_view("Microstrains", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::OTHER);
}
void MarginalBoneLoss::initialize(){

}
void MarginalBoneLoss::run(){
    if(this->elem_num == 0){
        size_t id = 0;
        for(const auto& g:mesh->geometries){
            if(id == this->geom_id){
                mand_geom_offset = elem_num;
                density_num = g->mesh.size();
            }
            elem_num += g->mesh.size();
            ++id;
        }
    }

    std::vector<double> density, max_density;
    std::vector<double> density_vector_view(this->elem_num, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> stress(this->elem_num, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> growth(this->elem_num, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> strains(this->elem_num, std::numeric_limits<double>::quiet_NaN());

    std::vector<double> u(this->mesh->max_dofs);

    field::ImplantRegion* field = nullptr;
    field::PrincipalStress* dir_field = nullptr;
    for(auto& f:this->proj->fields){
        if(f->get_type() == Field::Type::SCALAR){
            ScalarField* sf = static_cast<ScalarField*>(f.get());
            if(sf->get_class() == ScalarField::Class::IMPLANT_REGION){
                field = static_cast<field::ImplantRegion*>(sf);
            }
        } else if(f->get_type() == Field::Type::COORDINATE){
            CoordinateField* cf = static_cast<CoordinateField*>(f.get());
            if(cf->get_class() == CoordinateField::Class::PRINCIPAL_STRESS){
                dir_field = static_cast<field::PrincipalStress*>(cf);
            }
        }
    }
    logger::log_assert(field != nullptr, logger::ERROR, "marginal_bone_loss simulation requires an implant_region field");
    field->freeze(density, max_density);
    max_density.clear();
    std::copy(density.begin(), density.end(), density_vector_view.begin() + mand_geom_offset);
    std::fill(max_density.begin(), max_density.end(), MAX_RHO);
    this->density_view->update_view(density_vector_view);

    dir_field->set_max_it(1);

    const double MIN_DENSITY = 1e-1;

    double time_count = 0;
    while(time_count < this->time_limit){

        logger::quick_log("");
        logger::quick_log("Time:", time_count);

        this->fem->update_materials();
        this->mesh->apply_boundary_conditions(this->proj->forces, this->proj->supports, this->proj->springs, this->proj->internal_loads, this->proj->sub_problems);
        dir_field->generate();
        this->fem->generate_matrix(this->mesh);
        this->fem->calculate_displacements_global(this->mesh, this->mesh->load_vector, u);

        {
            auto s_it = stress.begin();
            size_t D_offset = 0;
            for(const auto& g:this->mesh->geometries){
                g->get_stresses(u, false, D_offset, this->fem->D_vec, s_it);
            }
            this->stress_view->update_view(stress);
        }

        for(size_t i = 0; i < density_num; ++i){
            const auto& e = this->mandible->mesh[i];
            const auto strain = 1e6*e->get_strain_vector(e->get_centroid(), u);
            const double gi = this->get_density_variation(strain)*this->time_step;
            const double old_d = density[i];
            density[i] += gi;
            if(density[i] < MIN_DENSITY){
                density[i] = MIN_DENSITY;
            } else if(density[i] > this->MAX_RHO){
                density[i] = max_density[i];
            }
            growth[this->mand_geom_offset + i] = density[i] - old_d;
            density_vector_view[i + this->mand_geom_offset] = density[i];
            strains[i + this->mand_geom_offset] = this->c[3]*this->LHS_3D(3, strain);
        }
        field->set_values(density);
        this->density_view->update_view(density_vector_view);

        this->growth_view->update_view(growth);
        this->strain_view->update_view(strains);

        time_count += this->time_step;
    }
}

double MarginalBoneLoss::get_density_variation(const StrainVector3D& eps) const{
    const double lhs = this->LHS_3D(3, eps);

    if(lhs < x0){
        return dv0*(lhs - x0);
    } else if(lhs < x1){
        return 0;
    } else {
        return (this->coeff[0] + this->coeff[1]*lhs + this->coeff[2]*lhs*lhs);
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
            DataEntry{.name = "time_limit", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "mandible_geometry", .type = TYPE_INT, .required = true},
            DataEntry{.name = "intervals", .type = TYPE_ARRAY, .required = true,
                         .array_data = std::shared_ptr<ArrayRequirements>(
                             new ArrayRequirements{
                                 .size = MarginalBoneLoss::RANGE_NUM,
                                 .type = TYPE_DOUBLE
                             }
                         ),
                     },
            DataEntry{.name = "compression", .type = TYPE_DOUBLE, .required = true,},
            DataEntry{.name = "traction", .type = TYPE_DOUBLE, .required = true,},
            DataEntry{.name = "shear", .type = TYPE_DOUBLE, .required = true,},
        }
    }
);

}
