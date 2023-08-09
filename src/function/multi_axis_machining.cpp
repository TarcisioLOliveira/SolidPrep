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

#include "function/multi_axis_machining.hpp"
#include "logger.hpp"
#include "projection/threshold.hpp"
#include <Eigen/src/Core/util/Constants.h>
#include <numeric>
#include <set>

namespace function{

MultiAxisMachining::MultiAxisMachining(const Meshing* const mesh, const DensityFilter* const filter, gp_Pnt center, gp_Dir axis, double v_norm, double L, double beta)
    : mesh(mesh), filter(filter), center(center), axis(axis), v_norm(v_norm), L(L), beta(beta), proj(projection::Threshold::Parameter{50, 50, 0, 0}, 0.65){}

void MultiAxisMachining::initialize_views(Visualization* viz){
    this->shadow_view = viz->add_view("Shadows", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::DENSITY);
}

void MultiAxisMachining::initialize(const Optimizer* const op){
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    const size_t num_nodes_bound = mesh->elem_info->get_boundary_nodes_per_element();
    size_t N = 0;
    if(this->mesh->elem_info->get_problem_type() == utils::PROBLEM_TYPE_2D){
        N = 2;
    } else if(this->mesh->elem_info->get_problem_type() == utils::PROBLEM_TYPE_3D){
        N = 3;
    }
    size_t elem_num = 0;
    size_t phi_size = 0;
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n = e->nodes[i];
                    this->id_mapping[n->id] = 0;
                }
            }
            elem_num += g->mesh.size();
        }
    }
    for(auto& e:mesh->boundary_elements){
        const gp_Pnt p = e.get_centroid(num_nodes_bound);
        const gp_Vec pc(this->center, p);
        const gp_Vec vv = pc - pc.Dot(this->axis)*this->axis;
        const gp_Dir vn(vv);
        if(vn.Dot(e.normal) > 0){
            for(size_t i = 0; i < num_nodes_bound; ++i){
                this->id_mapping[e.nodes[i]->id] = -1;
            }
        }
    }
    size_t id = 0;
    for(auto& n:this->mesh->node_list){
        if(this->id_mapping.count(n->id) && this->id_mapping.at(n->id) > -1){
            this->id_mapping[n->id] = id;
            ++id;
        }
    }
    phi_size = id;
    this->b_grad.resize(phi_size);
    this->b.resize(phi_size);
    this->diff.resize(elem_num,0);
    this->Hgrad.resize(elem_num,0);
    this->shadows.resize(elem_num, 0);
    this->shadows_grad.resize(elem_num, 0);
    this->Phi = Eigen::SparseMatrix<double>(phi_size, phi_size);
}

double MultiAxisMachining::calculate(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x){
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    if(!this->first_time){
        std::fill(this->Phi.valuePtr(), this->Phi.valuePtr() + this->Phi.nonZeros(), 0);
    }
    std::fill(this->b.begin(), this->b.end(), 0);
    auto x_it = x.cbegin();
    std::vector<double> v(3,0);
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            auto N = g->mesh.front()->helmholtz_vector(this->mesh->thickness);
            for(const auto& e:g->mesh){
                const gp_Pnt p = e->get_centroid();
                const gp_Vec pc(this->center, p);
                const gp_Vec vv = pc - pc.Dot(this->axis)*this->axis;
                const double vmag = vv.Magnitude();
                const gp_Dir vn(vv);
                v[0] = this->v_norm*vn.X();
                v[1] = this->v_norm*vn.Y();
                v[2] = this->v_norm*vn.Z();

                double dv = 0;
                for(size_t i = 0; i < 3; ++i){
                    const double ai = this->axis.Coord(1+i);
                    dv += (1.0/vmag - v[i]*v[i]/(vmag*vmag*vmag))*(1 - ai*ai);
                }
                dv *= this->v_norm;
                auto psi_e = e->get_phi_radial(this->mesh->thickness, this->beta, this->L, v, dv, *x_it);
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& ni = e->nodes[i];
                    const long id1 = this->id_mapping[ni->id];
                    if(id1 < 0){
                        continue;
                    }
                    for(size_t j = 0; j < num_nodes; ++j){
                        const auto& nj = e->nodes[j];
                        const long id2 = this->id_mapping[nj->id];
                        if(id2 < 0){
                            continue;
                        }
                        // Is the ordering correct?
                        this->Phi.coeffRef(id1, id2) += psi_e[i*num_nodes + j];
                    }
                }

                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& ni = e->nodes[i];
                    const long id1 = this->id_mapping[ni->id];
                    if(id1 < 0){
                        continue;
                    }
                    this->b[id1] -= this->beta*(*x_it)*N[i];
                }

                x_it += num_den;
            }
        }
    }

    if(this->first_time){
        this->Phi.makeCompressed();
        this->first_time = false;
        this->solver.analyzePattern(Phi);
    }

    this->solver.factorize(this->Phi);
    Eigen::VectorXd psi = this->solver.solve(this->b);
    x_it = x.cbegin();
    auto d_it = this->diff.begin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& ni = e->nodes[i];
                    const long id1 = this->id_mapping[ni->id];
                    if(id1 < 0){
                        continue;
                    }
                    *d_it += psi[id1]/num_nodes;
                }
                *d_it -= *x_it;
                ++d_it;
                x_it += num_den;
            }
        }
    }

    this->proj.project_densities(this->diff);

    return std::accumulate(this->diff.begin(), this->diff.end(), 0);
}

double MultiAxisMachining::calculate_with_gradient(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    std::fill(this->b.begin(), this->b.end(), 0);
    std::fill(this->b_grad.begin(), this->b_grad.end(), 0);
    auto x_it = x.cbegin();
    std::vector<double> v(3,0);
    // b vector, which is the same for every axis
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            auto N = g->mesh.front()->helmholtz_vector(this->mesh->thickness);
            for(const auto& e:g->mesh){
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& ni = e->nodes[i];
                    const long id1 = this->id_mapping[ni->id];
                    if(id1 < 0){
                        continue;
                    }
                    this->b[id1] -= this->beta*(*x_it)*N[i];
                }

                x_it += num_den;
            }
        }
    }
    // Make global matrix for each axis then solve it
    for(size_t i = 0; i < this->axes.size(); ++i){
        if(!this->first_time){
            std::fill(this->Phi.valuePtr(), this->Phi.valuePtr() + this->Phi.nonZeros(), 0);
        }
        x_it = x.cbegin();
        v[0] = this->axes[i].X();
        v[1] = this->axes[i].Y();
        v[2] = this->axes[i].Z();
        for(const auto& g:mesh->geometries){
            if(g->do_topopt){
                const size_t num_den = g->number_of_densities_needed();
                auto N = g->mesh.front()->helmholtz_vector(this->mesh->thickness);
                for(const auto& e:g->mesh){
                    auto psi_e = e->get_phi_unidirectional(this->mesh->thickness, this->beta*(*x_it), this->L, v, this->v_norm);
                    for(size_t i = 0; i < num_nodes; ++i){
                        const auto& ni = e->nodes[i];
                        const long id1 = this->id_mapping[ni->id];
                        if(id1 < 0){
                            continue;
                        }
                        for(size_t j = 0; j < num_nodes; ++j){
                            const auto& nj = e->nodes[j];
                            const long id2 = this->id_mapping[nj->id];
                            if(id2 < 0){
                                continue;
                            }
                            this->Phi.coeffRef(id1, id2) += psi_e[i*num_nodes + j];
                        }
                    }

                    x_it += num_den;
                }
            }
        }

        if(this->first_time){
            this->Phi.makeCompressed();
            this->first_time = false;
            this->solver.analyzePattern(Phi);
        }

        this->solver.factorize(this->Phi);
        logger::log_assert(this->solver.info() == Eigen::Success, logger::ERROR, "matrix decomposition failed");
        Eigen::VectorXd psi = this->solver.solve(this->b);

        x_it = x.cbegin();
        auto d_it = this->diff[i].begin();
        for(const auto& g:mesh->geometries){
            if(g->do_topopt){
                const size_t num_den = g->number_of_densities_needed();
                for(const auto& e:g->mesh){
                    for(size_t i = 0; i < num_nodes; ++i){
                        const auto& ni = e->nodes[i];
                        const long id1 = this->id_mapping[ni->id];
                        if(id1 < 0){
                            continue;
                        }
                        *d_it += psi[id1]/num_nodes;
                    }
                    *d_it -= *x_it;
                    ++d_it;
                    x_it += num_den;
                }
            }
        }

        std::fill(this->Hgrad[i].begin(), this->Hgrad[i].end(), 1);
        this->proj.project_gradient(this->Hgrad[i], this->diff[i]);
        this->proj.project_densities(this->diff[i]);
    }

    std::copy(this->diff[0].begin(), this->diff[0].end(), this->shadows.begin());
    std::copy(this->Hgrad[0].begin(), this->Hgrad[0].end(), this->shadows_grad.begin());
    if(this->axes.size() > 1){
        for(size_t i = 1; i < this->axes.size(); ++i){
            for(size_t j = 0; j < this->shadows.size(); ++j){
                this->shadows[j] *= this->diff[i][j];
                this->shadows_grad[j] *= this->diff[i][j];
            }
        }
        for(size_t j = 0; j < this->shadows.size(); ++j){
            for(size_t i = 1; i < this->axes.size(); ++i){
                double tmp = this->Hgrad[i][j];
                for(size_t k = 0; k < this->axes.size(); ++k){
                    if(k != i){
                        tmp *= this->diff[k][j];
                    }
                }
                this->shadows_grad[j] += tmp;
            }
        }
    }


    this->shadow_view->update_view(this->shadows);

    const double result = std::accumulate(this->shadows.begin(), this->shadows.end(), 0);

    auto dg_it = this->shadows_grad.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& ni = e->nodes[i];
                    const long id1 = this->id_mapping[ni->id];
                    if(id1 < 0){
                        continue;
                    }
                    b_grad[id1] += (*dg_it)/num_nodes;
                }
                ++dg_it;
            }
        }
    }

    Eigen::VectorXd psi_tilde = this->solver.solve(this->b_grad);

    auto g_it = grad.begin();
    auto hg_it = this->Hgrad.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            auto phi_e = g->mesh.front()->get_phi_grad(this->mesh->thickness, this->beta);
            auto b_e = g->mesh.front()->helmholtz_vector(this->mesh->thickness);
            for(auto& i:b_e){
                i *= this->beta;
            }
            for(const auto& e:g->mesh){
                *g_it = 0;
                double psi_tmp = 0;
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& ni = e->nodes[i];
                    const long id1 = this->id_mapping[ni->id];
                    if(id1 < 0){
                        continue;
                    }
                    psi_tmp = 0;
                    for(size_t j = 0; j < num_nodes; ++j){
                        const auto& nj = e->nodes[j];
                        const long id2 = this->id_mapping[nj->id];
                        if(id2 < 0){
                            continue;
                        }
                        // Is the ordering correct?
                        psi_tmp += phi_e[i*num_nodes + j]*psi[id2];
                    }
                    *g_it += psi_tilde[id1]*psi_tmp;
                    *g_it += psi_tilde[id1]*b_e[i];
                    psi_tmp = 0;
                }
                *g_it -= *hg_it;
                ++hg_it;
                g_it += num_den;
            }
        }
    }

    return result;
}

}
