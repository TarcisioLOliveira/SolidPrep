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
#include "field/trabecular_cutoff.hpp"
#include "field/apparent_density.hpp"
#include "logger.hpp"
#include "math/matrix.hpp"
#include <algorithm>
#include <numeric>
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
    MAX_RHO(data.get_double("rho_max_cort")),
    MAX_RHO_TRAB(data.get_double("rho_max_trab")),
    rho_d(data.get_double("growth_cort")),
    rho_d_t(data.get_double("growth_trab")),
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
    this->coeff[2] = this->rho_d;
    this->coeff2 = math::Vector(3);
    this->coeff2[2] = this->rho_d_t;

    //eps[0] = -this->c[0];
    //this->x0 = this->LHS_3D(3, eps);
    //eps[0] = -this->c[1];
    //this->x1 = this->LHS_3D(3, eps);
    //eps[0] = -this->c[2];
    //this->x2 = this->LHS_3D(3, eps);
    this->x0 = this->c[0]/this->c[3];
    this->x1 = this->c[1]/this->c[3];
    this->x2 = this->c[2]/this->c[3];

    const double xm = (x1 + x2)/2.0;
    math::Matrix A(
        {1, x1, x1*x1,
         1, x2, x2*x2,
         1, xm, xm*xm},
        3, 3);

    math::LU LU(A);
    LU.solve(this->coeff);
    LU.solve(this->coeff2);

    this->dv0 = coeff[1] + 2*coeff[2]*x1;
    this->dv0t = coeff2[1] + 2*coeff2[2]*x1;
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
    this->density_view = viz->add_view("Density", spview::defs::ViewType::NODAL, spview::defs::DataType::OTHER);
    this->growth_view = viz->add_view("Growth (Density)", spview::defs::ViewType::NODAL, spview::defs::DataType::OTHER);
    this->strain_view = viz->add_view("Microstrains", spview::defs::ViewType::NODAL, spview::defs::DataType::OTHER);
    this->reactions_view = viz->add_view("Reactions", spview::defs::ViewType::VECTOR, spview::defs::DataType::OTHER);
}
void MarginalBoneLoss::initialize(){

}
void MarginalBoneLoss::run(){
    //size_t total_node_num = this->mesh->node_list.size();
    if(this->elem_num == 0){
        size_t id = 0;
        for(const auto& g:mesh->geometries){
            if(id == this->geom_id){
                mand_geom_offset = elem_num;
                density_num = g->node_list.size();
                //this->v.resize(density_num);
                //for(size_t i = 0; i < density_num; ++i){
                //    this->v[i] = g->mesh[i]->get_volume(this->mesh->thickness);
                //}
            }
            elem_num += g->mesh.size();
            ++id;
        }
        //this->v_ref = std::accumulate(this->v.begin(), this->v.end(), 0.0)/density_num;
    }

    const auto NODES_PER_ELEM = this->proj->topopt_element->get_nodes_per_element();
    const auto elem_shape = this->proj->topopt_element->get_shape_type();
    const auto integ_order = this->proj->topopt_element->get_element_order();
    const auto S_SIZE = this->proj->topopt_element->get_D_dimension();

    field::ImplantRegion* field = nullptr;
    field::PrincipalStress* dir_field = nullptr;
    const field::ApparentDensity* ref_field = nullptr;
    const field::TrabecularCutoff* trab_field = nullptr;
    for(auto& f:this->proj->fields){
        if(f->get_type() == Field::Type::SCALAR){
            ScalarField* sf = static_cast<ScalarField*>(f.get());
            if(sf->get_class() == ScalarField::Class::IMPLANT_REGION){
                field = static_cast<field::ImplantRegion*>(sf);
            } else if(sf->get_class() == ScalarField::Class::APPARENT_DENSITY){
                ref_field = static_cast<field::ApparentDensity*>(sf);
            } else if(sf->get_class() == ScalarField::Class::TRABECULAR_CUTOFF){
                trab_field = static_cast<field::TrabecularCutoff*>(sf);
            }
        } else if(f->get_type() == Field::Type::COORDINATE){
            CoordinateField* cf = static_cast<CoordinateField*>(f.get());
            if(cf->get_class() == CoordinateField::Class::PRINCIPAL_STRESS){
                dir_field = static_cast<field::PrincipalStress*>(cf);
            }
        }
    }

    std::vector<double> density;
    std::vector<double> stress(this->elem_num, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> growth (density_num, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> strains(density_num, std::numeric_limits<double>::quiet_NaN());

    std::vector<double> u(this->mesh->max_dofs);

    logger::log_assert(field != nullptr, logger::ERROR, "marginal_bone_loss simulation requires an implant_region field");
    field->freeze(density);
    this->density_view->update_view(density, {this->mandible->id});


    std::map<const Node*, const MeshElement*> node_to_elem;
    std::map<size_t, size_t> id_pos_map;
    std::vector<math::Vector> nodal_strains;
    std::vector<double> nodal_divisors;
    std::vector<bool> nodal_trab;
    {
        size_t new_id = 0;
        for(auto& n:this->mandible->node_list){
            id_pos_map[n->id] = new_id;
            ++new_id;
        }
    }
    if(!(elem_shape == Element::Shape::TRI && integ_order == 1)){
        for(auto& e:this->mandible->mesh){
            for(size_t i = 0; i < NODES_PER_ELEM; ++i){
                node_to_elem[e->nodes[i]] = e.get();
            }
        }
    } else {
        nodal_strains.resize(density_num, math::Vector(S_SIZE));
        nodal_divisors.resize(density_num);
        nodal_trab.resize(density_num);
        for(size_t i = 0; i < density_num; ++i){
            const auto n = this->mandible->node_list[i];
            node_to_elem[n] = nullptr;
            const double ref_rho = ref_field->get(n);
            nodal_trab[i] = (ref_rho <= MAX_RHO_TRAB);
        }
        for(auto& e:this->mandible->mesh){
            const auto c = e->get_centroid();
            const bool center_trab = (trab_field->get(e.get(), c) > 0.5);
            for(size_t i = 0; i < NODES_PER_ELEM; ++i){
                const auto& n = e->nodes[i];
                auto& ei = node_to_elem[n];
                const auto id = id_pos_map[n->id];
                nodal_divisors[id] += 1;

                if(ei == nullptr && center_trab == nodal_trab[id]){
                    ei = e.get();
                }
            }
        }
    }

    dir_field->set_max_it(1);

    const double MIN_DENSITY = 0.05;

    double time_count = 0;
    while(time_count < this->time_limit){

        logger::quick_log("");
        logger::quick_log("Time:", time_count);

        //this->mesh->apply_boundary_conditions(this->proj->forces, this->proj->supports, this->proj->springs, this->proj->internal_loads, this->proj->contact_springs, this->proj->internal_loads2, this->proj->internal_loads3, this->proj->sub_problems);
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

        if(elem_shape == Element::Shape::TRI && integ_order == 1){
            #pragma omp parallel for
            for(auto& v:nodal_strains){
                v.fill(0);
            }
            for(auto& e:this->mandible->mesh){
                const auto s = e->get_strain_vector(e->get_centroid(), u);
                for(size_t i = 0; i < NODES_PER_ELEM; ++i){
                    const auto& n = e->nodes[i];
                    const auto id = id_pos_map[n->id];
                    nodal_strains[id] += 1e6*s;
                }
            }
            #pragma omp parallel for
            for(size_t i = 0; i < density_num; ++i){
                //const auto& n = this->mandible->node_list[i];
                //const auto& e = node_to_elem[n];
                nodal_strains[i] /= nodal_divisors[i];

                const auto& strain = nodal_strains[i];
                const double old_d = density[i];

                //const auto& D = this->mandible->materials.get_D(e, n->point);
                //const auto& D0 = this->mandible->materials.get_D_orig(e, n->point);

                const bool is_trab = nodal_trab[i];

                //const double gi = this->get_density_variation(strain, D0, D, is_trab)*this->time_step;
                const double gi = this->get_density_variation(strain, is_trab)*this->time_step;

                density[i] += gi;
                if(density[i] < MIN_DENSITY){
                    density[i] = MIN_DENSITY;
                } else {
                    if(is_trab){
                        density[i] = std::min(density[i], this->MAX_RHO_TRAB);
                    } else {
                        density[i] = std::min(density[i], this->MAX_RHO);
                    }
                }
                growth[i] = density[i] - old_d;
                //strains[i] = this->c[3]*this->Ud(3, strain, D0, D);
                strains[i] = this->c[3]*this->LHS_3D(3, strain);
            }
        } else {
            #pragma omp parallel for
            for(size_t i = 0; i < density_num; ++i){
                const auto& n = this->mandible->node_list[i];
                const auto& e = node_to_elem[n];
                const auto strain = 1e6*e->get_strain_vector(n->point, u);
                const double old_d = density[i];

                //const auto& D = this->mandible->materials.get_D(e, n->point);
                //const auto& D0 = this->mandible->materials.get_D_orig(e, n->point);

                const auto trab_num = trab_field->get(e, n->point);
                const bool is_trab = (trab_num > 0.5);

                //const double gi = this->get_density_variation(strain, D0, D, is_trab)*this->time_step;
                const double gi = this->get_density_variation(strain, is_trab)*this->time_step;

                density[i] += gi;
                if(density[i] < MIN_DENSITY){
                    density[i] = MIN_DENSITY;
                } else {
                    if(is_trab){
                        density[i] = std::min(density[i], this->MAX_RHO_TRAB);
                    } else {
                        density[i] = std::min(density[i], this->MAX_RHO);
                    }
                }
                growth[n->id] = density[i] - old_d;
                //strains[n->id] = this->c[3]*this->Ud(3, strain, D0, D);
                strains[i] = this->c[3]*this->LHS_3D(3, strain);
            }
        }
        field->set_values(density);
        this->density_view->update_view(density, {this->mandible->id});

        this->growth_view->update_view(growth, {this->mandible->id});
        this->strain_view->update_view(strains, {this->mandible->id});

        time_count += this->time_step;
        
        sleep(5);
    }
    
    auto reactions = proj->topopt_fea->calculate_reactions(proj->topopt_mesher.get(), u);
    this->reactions_view->update_view(reactions);
}

double MarginalBoneLoss::get_density_variation(const StrainVector3D& eps, bool trabecular) const{
    const double lhs = this->LHS_3D(3, eps);

    if(!trabecular){
        if(lhs < x0){
            return dv0*(lhs - x0);
        } else if(lhs < x1){
            return 0;
        } else {
            return (this->coeff[0] + this->coeff[1]*lhs + this->coeff[2]*lhs*lhs);
        }
    } else {
        if(lhs < x0){
            return dv0t*(lhs - x0);
        } else if(lhs < x1){
            return 0;
        } else {
            return (this->coeff2[0] + this->coeff2[1]*lhs + this->coeff2[2]*lhs*lhs);
        }
    }
}
//double MarginalBoneLoss::get_density_variation(const StrainVector3D& eps, const math::Matrix& D0, const math::Matrix& D, bool trabecular) const{
//    const double U = this->Ud(3, eps, D0, D);
//
//    if(!trabecular){
//        if(U < x0){
//            return dv0*(U - x0);
//        } else if(U < x1){
//            return 0;
//        } else {
//            return (this->coeff[0] + this->coeff[1]*U + this->coeff[2]*U*U);
//        }
//    } else {
//        if(U < x0){
//            return dv0t*(U - x0);
//        } else if(U < x1){
//            return 0;
//        } else {
//            return (this->coeff2[0] + this->coeff2[1]*U + this->coeff2[2]*U*U);
//        }
//    }
//}

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
            DataEntry{.name = "rho_max_cort", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "rho_max_trab", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "growth_cort", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "growth_trab", .type = TYPE_DOUBLE, .required = true},
        }
    }
);

}
