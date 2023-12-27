/*
 *   Copyright (C) 2023 Tarcísio Ladeia de Oliveira.
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

#include <set>
#include "field/orthotropic_flow.hpp"
#include "general_solver/mumps_general.hpp"
#include "utils/sparse_matrix.hpp"
#include "logger.hpp"
#include "utils.hpp"

namespace field{

OrthotropicFlow::OrthotropicFlow(const MeshElementFactory* elem_info, std::vector<Geometry*> geoms, std::vector<CrossSection> entries, std::vector<double> coeffs, double thickness, bool show):
    elem_info(elem_info), geoms(std::move(geoms)), 
    entries(std::move(entries)), coeffs(std::move(coeffs)),
    thickness(thickness), show(show),
    DIM((elem_info->get_problem_type() == utils::PROBLEM_TYPE_2D) ? 2 : 3),
    NODES_PER_ELEM(elem_info->get_nodes_per_element()){

}

void OrthotropicFlow::initialize_views(Visualization* viz){
    if(this->show){
        this->L_view = viz->add_view("Longitudinal Vector", spview::defs::ViewType::VECTOR, spview::defs::DataType::OTHER);
        this->R_view = viz->add_view("Radial Vector", spview::defs::ViewType::VECTOR, spview::defs::DataType::OTHER);
        this->T_view = viz->add_view("Tangential Vector", spview::defs::ViewType::VECTOR, spview::defs::DataType::OTHER);
    }
}

void OrthotropicFlow::display_views() const{
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
                    L_vecs[id*3 + i] = dirs[0].Coord(1+i);
                    R_vecs[id*3 + i] = dirs[1].Coord(1+i);
                    T_vecs[id*3 + i] = dirs[2].Coord(1+i);
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
                            L_vecs[id*3 + j] = dirs[0].Coord(1+j);
                            R_vecs[id*3 + j] = dirs[1].Coord(1+j);
                            T_vecs[id*3 + j] = dirs[2].Coord(1+j);
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

void OrthotropicFlow::generate(){
    logger::quick_log("Generating orthotropic flow...");
    const size_t num_nodes = this->elem_info->get_nodes_per_element();
    const size_t bound_num = elem_info->get_boundary_nodes_per_element();

    size_t elem_num = 0;
    for(const auto& g:this->geoms){
        elem_num += g->boundary_mesh.size();
    }
    this->elem_mult_long.resize(elem_num, 0);
    this->elem_mult_rad.resize(elem_num, -1);

    size_t mult_offset = 0;
    for(const auto& g:this->geoms){
        #pragma omp parallel for
        for(size_t i = 0; i < g->boundary_mesh.size(); ++i){
            const auto c = g->boundary_mesh[i]->get_centroid(bound_num);
            for(size_t j = 0; j < this->entries.size(); ++j){
                if(this->entries[j].is_inside(c)){
                    this->elem_mult_long[mult_offset + i] = this->coeffs[j];
                    elem_mult_rad[mult_offset+i] = 0;
                    break;
                }
            }
        }

        mult_offset += g->boundary_mesh.size();
    }

    if(this->show){
        std::set<UniqueNodeForView> uniques;
        if(!(this->elem_info->get_element_order() == 1 && 
             this->elem_info->get_shape_type() == Element::Shape::TRI)){
            for(const auto& g:this->geoms){
                for(const auto& e:g->mesh){
                    for(size_t i = 0; i < bound_num; ++i){
                        uniques.insert(UniqueNodeForView{e->nodes[i], e.get()});
                    }
                }
            }
        }
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

    const size_t phi_size = new_id;
    this->longitudinal.resize(phi_size);
    this->radial.resize(phi_size);

    general_solver::MUMPSGeneral solver;
    solver.initialize_matrix(true, phi_size);

    std::vector<double> b(phi_size, 0);;

    std::vector<double> A{1, 0, 0,
                          0, 1, 0,
                          0, 0, 1};

    for(const auto& g:this->geoms){
        for(const auto& e:g->mesh){
            const Eigen::MatrixXd M_e = e->diffusion_1dof(thickness, A) + 1e-4*e->absorption_1dof(thickness);
            std::vector<double> M_ev(num_nodes*num_nodes, 0);
            std::copy(M_e.data(), M_e.data()+M_ev.size(), M_ev.begin());
            std::vector<long> pos(num_nodes);

            for(size_t i = 0; i < num_nodes; ++i){
                const long id1 = id_pos_map.at(e->nodes[i]->id);
                pos[i] = id1;
            }

            solver.add_element(M_ev, pos);
        }
    }

    solver.compute();

    size_t it = 0;
    for(const auto& g:this->geoms){
        for(const auto& e:g->boundary_mesh){
            const Eigen::VectorXd N = this->elem_mult_long[it]*e->parent->flow_1dof(thickness, e->nodes);

            for(size_t i = 0; i < num_nodes; ++i){
                const long id1 = id_pos_map.at(e->parent->nodes[i]->id);
                b[id1] += N[i];
            }

            ++it;
        }
    }

    solver.solve(b);
    std::copy(b.begin(), b.end(), this->longitudinal.begin());

    std::fill(b.begin(), b.end(), 0);
    it = 0;
    for(const auto& g:this->geoms){
        for(const auto& e:g->boundary_mesh){
            const Eigen::VectorXd N = this->elem_mult_rad[it]*e->parent->flow_1dof(thickness, e->nodes);

            for(size_t i = 0; i < num_nodes; ++i){
                const long id1 = id_pos_map.at(e->parent->nodes[i]->id);
                b[id1] += N[i];
            }

            ++it;
        }
    }

    solver.solve(b);
    std::copy(b.begin(), b.end(), this->radial.begin());

    this->elem_mult_long.clear();
    this->elem_mult_rad.clear();

    logger::quick_log("Done.");
}

std::array<gp_Dir, 3> OrthotropicFlow::get_array(const MeshElement* e, const gp_Pnt& p) const{
    const auto grad = e->get_nodal_density_gradient(p);
    std::vector<double> L(3, 0);
    std::vector<double> R(3, 0);
    for(size_t i = 0; i < DIM; ++i){
        for(size_t j = 0; j < NODES_PER_ELEM; ++j){
            const size_t id = this->id_pos_map.at(e->nodes[j]->id);
            L[i] += grad[i*NODES_PER_ELEM + j]*this->longitudinal[id];
            R[i] += grad[i*NODES_PER_ELEM + j]*this->radial[id];
        }
    }

    gp_Dir Ld(L[0], L[1], L[2]);
    gp_Dir Rd(R[0], R[1], R[2]);

    return {Ld, Rd, Ld.Crossed(Rd)};
}
Eigen::Matrix<double, 3, 3> OrthotropicFlow::get_matrix(const MeshElement* e, const gp_Pnt& p) const{
    const auto d = this->get_array(e, p);

    return Eigen::Matrix<double, 3, 3>
            {
                {d[0].X(), d[1].X(), d[2].X()},
                {d[0].Y(), d[1].Y(), d[2].Y()},
                {d[0].Z(), d[1].Z(), d[2].Z()}
            };
}

}
