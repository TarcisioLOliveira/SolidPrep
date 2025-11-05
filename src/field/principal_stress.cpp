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

#include "field/principal_stress.hpp"
#include "logger.hpp"
#include "math/matrix.hpp"
#include "project_data.hpp"
#include "utils/delayed_pointer.hpp"
#include <functional>

namespace field{

PrincipalStress::PrincipalStress(const projspec::DataMap& data):
    proj_data(data.proj),
    initial_material(data.proj->get_material(data.get_string("initial_material"))),
    elem_info(data.proj->topopt_element),
    geoms(),
    thickness(data.proj->thickness),
    show(data.get_bool("display", true)),
    DIM((elem_info->get_problem_type() == utils::PROBLEM_TYPE_2D) ? 2 : 3),
    NODES_PER_ELEM(elem_info->get_nodes_per_element()),
    max_it(data.get_int("max_it")),
    use_strain(data.get_bool("use_strain"))
    {
   
    if(!data.get_bool("STUB")){ 
        logger::log_assert(this->initial_material.get() != nullptr &&
                initial_material->get_class() == Material::Class::ISOTROPIC,
                logger::ERROR,
                "initial material chosen for principal stress field is not isotropic");
    }
}

void PrincipalStress::initialize_views(Visualization* viz){
    if(this->show && this->geoms.size() > 0){
        this->R_view = viz->add_view("Radial Vector", spview::defs::ViewType::VECTOR, spview::defs::DataType::OTHER);
        this->T_view = viz->add_view("Tangential Vector", spview::defs::ViewType::VECTOR, spview::defs::DataType::OTHER);
        this->L_view = viz->add_view("Longitudinal Vector", spview::defs::ViewType::VECTOR, spview::defs::DataType::OTHER);
    }
}

void PrincipalStress::generate(){
    const size_t num_nodes = this->elem_info->get_nodes_per_element();

    this->geoms.clear();

    auto& proj_geoms = this->proj_data->geometries;
    std::set<Geometry*> picked_geoms;
    for(auto& g:proj_geoms){
        auto mats = g->materials.get_materials();
        for(auto& m:mats){
            if(m->get_field() == this){
                picked_geoms.insert(g.get());
                break;
            }
        }
    }
    this->geoms.insert(this->geoms.begin(), picked_geoms.begin(), picked_geoms.end());

    if(this->geoms.size() == 0){
        return;
    }

    if(this->show){
        if(!(this->elem_info->get_element_order() == 1 && 
             this->elem_info->get_shape_type() == Element::Shape::TRI)){

            std::set<UniqueNodeForView> uniques;
            for(const auto& g:this->geoms){
                for(const auto& e:g->mesh){
                    for(size_t i = 0; i < num_nodes; ++i){
                        uniques.insert(UniqueNodeForView{e->nodes[i], e.get()});
                    }
                }
            }
            this->node_view_list.insert(this->node_view_list.begin(), uniques.begin(), uniques.end());
        }

        std::set<size_t> unique_ids;

        for(const auto& g:this->geoms){
            for(const auto& e:g->mesh){
                for(size_t i = 0; i < num_nodes; ++i){
                    unique_ids.insert(e->nodes[i]->id);
                }
            }
        }

        size_t new_id = 0;
        for(const auto& id:unique_ids){
            this->id_pos_map[id] = new_id;
            ++new_id;
        }
        unique_ids.clear();
    }

    Meshing* mesh = this->proj_data->topopt_mesher.get();
    SolverManager* fem = this->proj_data->topopt_fea.get();

    std::vector<std::vector<utils::DelayedPointerView<Material>>> materials_backup(this->geoms.size());

    std::vector<utils::DelayedPointerView<Material>> iso_material{this->initial_material};
    std::vector<std::vector<utils::DelayedPointerView<Material>>> ortho_material(this->geoms.size());
    for(size_t i = 0; i < this->geoms.size(); ++i){
        auto& g = this->geoms[i];
        auto mats = g->materials.get_materials();
        materials_backup[i] = mats;
        for(auto& m:mats){
            if(m->get_field() == this){
                // Assumes there is a single material that uses this field
                // for each material list because it's much easier this
                // way
                ortho_material[i] = std::vector<utils::DelayedPointerView<Material>>{m};
                break;
            }
        }
        g->set_materials(iso_material);
    }
    // Check if it's necessary to update BCs at each iteration
    bool update_boundary_conditions = false;
    if(this->proj_data->internal_loads.size() > 0){
        const size_t int_size = this->proj_data->internal_loads.size();
        for(size_t i = 0; i < int_size; ++i){
            auto& inti = this->proj_data->internal_loads[i];
            if(inti.mat->get_field() == this){
                update_boundary_conditions = true;
                break;
            }
        }
    }
    if(this->proj_data->springs.size() > 0 && !update_boundary_conditions){
        const size_t int_size = this->proj_data->springs.size();
        for(size_t i = 0; i < int_size; ++i){
            auto& inti = this->proj_data->springs[i];
            if(inti.mat->get_field() == this){
                update_boundary_conditions = true;
                break;
            }
        }
    }


    this->u.resize(mesh->max_dofs, 0);
    auto u_old = this->u;
    double c_old = 0;

    double du = 1.0;
    size_t it = 0;
    logger::quick_log("Generating principal stress field....");
    logger::quick_log("");
    while(it < max_it){
        fem->update_materials();
        fem->generate_matrix(mesh);

        fem->calculate_displacements_global(mesh, mesh->load_vector, this->u);

        du = 0;
        for(size_t i = 0; i < mesh->max_dofs; ++i){
            du = std::max(du, std::abs(u[i] - u_old[i]));
            u_old[i] = u[i];
        }

        double c = 0;
        for(size_t i = 0; i < this->u.size(); ++i){
            c += this->u[i]*mesh->global_load_vector[i];
        }

        if(it == 0){
            this->first_run = false;
        }
        logger::quick_log("Iteration:", it, "du:", du, "c:", c, "dc:", std::abs(c - c_old), "dc/c:", std::abs(c - c_old)/std::abs(c));
        c_old = c;

        if(update_boundary_conditions){
            mesh->apply_boundary_conditions(this->proj_data->forces, this->proj_data->supports, this->proj_data->springs, this->proj_data->internal_loads, this->proj_data->sub_problems);
        }

        ++it;
    }

    for(size_t i = 0; i < this->geoms.size(); ++i){
        auto& g = this->geoms[i];
        g->set_materials(materials_backup[i]);
    }
    logger::quick_log("");
    logger::quick_log("Finished.");
}

void PrincipalStress::display_views() const{
    if(this->geoms.size() == 0){
        return;
    }
    if(this->show){
        const size_t nodes_per_elem = elem_info->get_nodes_per_element();
        const size_t node_num = this->id_pos_map.size();

        std::vector<size_t> ids(this->geoms.size(), 0);
        for(size_t i = 0; i < ids.size(); ++i){
            ids[i] = this->geoms[i]->id;
        }
        std::vector<double> L_vecs(node_num*3, 0);
        std::vector<double> R_vecs(node_num*3, 0);
        std::vector<double> T_vecs(node_num*3, 0);
        if(this->node_view_list.size() > 0){
            for(auto& n:this->node_view_list){
                const auto dirs = this->get_array(n.element, n.node->point);
                const size_t id = this->id_pos_map.at(n.node->id);
                for(size_t i = 0; i < 3; ++i){
                    R_vecs[id*3 + i] = dirs[0].Coord(1+i);
                    T_vecs[id*3 + i] = dirs[1].Coord(1+i);
                    L_vecs[id*3 + i] = dirs[2].Coord(1+i);
                }
            }
        } else {
            for(const auto& g:geoms){
                for(const auto& e:g->mesh){
                    const auto dirs = this->get_array(e.get(), e->get_centroid());
                    for(size_t i = 0; i < nodes_per_elem; ++i){
                        const auto n = e->nodes[i];
                        const size_t id = this->id_pos_map.at(n->id);
                        for(size_t j = 0; j < 3; ++j){
                            R_vecs[id*3 + j] = dirs[0].Coord(1+j);
                            T_vecs[id*3 + j] = dirs[1].Coord(1+j);
                            L_vecs[id*3 + j] = dirs[2].Coord(1+j);
                        }
                    }
                }
            }
            for(size_t i = 0; i < node_num; ++i){
                double L = 0, R = 0, T = 0;
                for(size_t j = 0; j < 3; ++j){
                    L += L_vecs[i*3 + j]*L_vecs[i*3 + j];
                    R += R_vecs[i*3 + j]*R_vecs[i*3 + j];
                    T += T_vecs[i*3 + j]*T_vecs[i*3 + j];
                }
                L = std::sqrt(L);
                R = std::sqrt(R);
                T = std::sqrt(T);
                for(size_t j = 0; j < 3; ++j){
                    L_vecs[i*3 + j] /= L;
                    R_vecs[i*3 + j] /= R;
                    T_vecs[i*3 + j] /= T;
                }
            }
        }
        this->L_view->update_view(L_vecs, ids);
        this->R_view->update_view(R_vecs, ids);
        this->T_view->update_view(T_vecs, ids);
    }
}

std::array<gp_Dir, 3> PrincipalStress::get_array(const MeshElement* e, const gp_Pnt& p) const{
    if(this->first_run){
        return {
            gp_Dir(1, 0, 0),
            gp_Dir(0, 1, 0),
            gp_Dir(0, 0, 1)
        };
    } else {
        auto M = this->get_matrix(e, p);
        return {
            gp_Dir(M(0,0), M(1,0), M(2,0)),
            gp_Dir(M(0,1), M(1,1), M(2,1)),
            gp_Dir(M(0,2), M(1,2), M(2,2)),
        };
    }
}
math::Matrix PrincipalStress::get_matrix(const MeshElement* e, const gp_Pnt& p) const{
    if(this->first_run){
        return 
            math::Matrix(
                {
                    1, 0, 0,
                    0, 1, 0,
                    0, 0, 1
                }, 3, 3);
    } else {
        math::Matrix tensor;
        if(this->use_strain){
            tensor = e->get_strain_tensor(p, this->u);
        } else {
            if(this->DIM == 3){
                const auto D = this->initial_material->stiffness_3D(e, p);
                tensor = e->get_stress_tensor(D, p, this->u);
            } else {
                const auto D = this->initial_material->stiffness_2D(e, p);
                tensor = e->get_stress_tensor(D, p, this->u);
            }
        }
        math::Eigen eigen(tensor);
        logger::log_assert(std::abs(eigen.eigenvectors.determinant()) > 1e-7,
                logger::ERROR,
                "eigen vector singular");
        logger::log_assert(
                (eigen.eigenvectors*eigen.eigenvectors.T()).is_equal(math::Matrix::identity(DIM)),
                logger::ERROR,
                "eigen vector not orthonormal");
        std::vector<double> abs_sort(DIM);
        for(size_t i = 0; i < DIM; ++i){
            abs_sort[i] = std::abs(eigen.eigenvalues[i]);
        }
        std::sort(abs_sort.begin(), abs_sort.end(), std::less<double>());
        std::vector<size_t> pos(DIM);
        for(size_t i = 0; i < DIM; ++i){
            for(size_t j = 0; j < DIM; ++j){
                if(abs_sort[i] == std::abs(eigen.eigenvalues[j])){
                    pos[i] = j;
                    break;
                }
            }
        }
        math::Matrix dirs(DIM, DIM);
        for(size_t i = 0; i < DIM; ++i){
            for(size_t j = 0; j < DIM; ++j){
                dirs(i, j) = eigen.eigenvectors(i, pos[j]);
            }
        }
        return math::Matrix::identity(DIM);
        return dirs;
    }
}

using namespace projspec;
const bool PrincipalStress::reg = Factory<Field>::add(
    [](const DataMap& data){
        return std::make_unique<PrincipalStress>(data);
    },
    ObjectRequirements{
        "principal_stress",
        {
            DataEntry{.name = "initial_material", .type = TYPE_STRING, .required = true},
            DataEntry{.name = "display", .type = TYPE_BOOL, .required = false},
            DataEntry{.name = "use_strain", .type = TYPE_BOOL, .required = false},
            DataEntry{.name = "max_it", .type = TYPE_INT, .required = true},
        }
    }
);
}
