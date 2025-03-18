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

namespace simulation{

MarginalBoneLoss::MarginalBoneLoss(Meshing* mesh, SolverManager* fem, Geometry* mandible, size_t geom_id, double time_step, double maximum_volume_variation, double time_limit, double maturation_rate, double pc, Range traction, Range compression, Range shear, std::vector<double> a0, std::vector<double> a1, std::vector<double> a2, std::vector<double> a3, utils::ProblemType type)
    :
    mesh(mesh), fem(fem), pc(pc),
    t(traction), c(compression), s(shear),
    K_e1(this->get_K_e1(t, c)),
    K_g(this->get_K_g(s)),
    K_e2(this->get_K_e2(t, c)),
    problem_type(type),
    a0(a0), a1(a1), a2(a2), a3(a3),
    eps_to_x(this->make_eps_to_x()),
    mandible(mandible), geom_id(geom_id),
    time_step(time_step), maximum_volume_variation(maximum_volume_variation), time_limit(time_limit),
    maturation_rate(maturation_rate)
{
    logger::log_assert(mandible->materials.get_materials()[0]->get_type() == Material::MANDIBLE,
            logger::ERROR,
            "geometry selected for marginal bone loss simulation must have a material of type \"mandible\"");
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

    material::Mandible* material = static_cast<material::Mandible*>(this->mandible->materials.get_materials()[0]);

    double time_count = 0;
    double alpha_count = 0;
    while(time_count < this->time_limit){

        logger::quick_log("");
        logger::quick_log("Time:", time_count);

        material->set_maturation_alpha(alpha_count);
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
            const double gi = this->get_density_variation(strain)*this->time_step*this->maturation_rate/v[i];
            growth[mand_geom_offset + i] = gi;
            density[i] += gi;
            if(density[i] < 0){
                density[i] = 0;
            } else if(density[i] > 1){
                density[i] = 1;
            }
        }
        this->growth_view->update_view(growth);

        alpha_count += this->maturation_rate;
        time_count += this->time_step;
    }
}

math::Vector MarginalBoneLoss::make_eps_to_x() const{
    math::Matrix M(RANGE_NUM, RANGE_NUM);
    math::Vector b(RANGE_NUM);

    StrainVector3D eps(6);

    for(size_t i = 0; i < RANGE_NUM; ++i){
        b[i] = i + 1;
        eps[0] = 1.0/(3*K_e2[i]);
        eps[1] = 1.0/(3*K_e2[i]);
        eps[2] = 1.0/(3*K_e2[i]);
        for(size_t j = 0; j < RANGE_NUM; ++j){

            const double LHS = this->LHS_3D(j, eps);
            M(i,j) = LHS;
        } 
    }

    math::LU LU(M);
    LU.solve(b);

    return b;
}

double MarginalBoneLoss::get_density_variation(const StrainVector3D& eps) const{
    StrainVector3D lhs(RANGE_NUM);
    for(size_t j = 0; j < RANGE_NUM; ++j){
        lhs[j] = this->LHS_3D(j, eps);
    }

    const auto poly_calc = [](const std::vector<double>& a, const double x)->double{
        double y = 0;
        for(size_t i = 0; i < a.size(); ++i){
            y += a[i]*std::pow(x, i);
        }
        return y;
    };

    const double x = eps_to_x.T()*lhs;
    if(x < 0){
        return poly_calc(this->a0, 0);
    } else if(x < 1){
        return poly_calc(this->a0, x);
    } else if(x < 2){
        return poly_calc(this->a1, x);
    } else if(x < 3){
        return poly_calc(this->a2, x);
    } else if(x < 4){
        return poly_calc(this->a3, x);
    } else {
        return poly_calc(this->a3, 4);
    }
}

}
