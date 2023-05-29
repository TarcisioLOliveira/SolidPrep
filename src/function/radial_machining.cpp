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

#include "function/radial_machining.hpp"
#include "logger.hpp"
#include "projection/heaviside.hpp"
#include "projection/threshold.hpp"
#include <Eigen/src/Core/util/Constants.h>
#include <limits>
#include <mpi.h>
#include <numeric>
#include <set>

namespace function{

RadialMachining::RadialMachining(const Meshing* const mesh, const DensityFilter* const filter, gp_Pnt center, gp_Dir axis, double v_norm, double beta)
        : mesh(mesh), filter(filter), center(center), axis(axis), v_norm(v_norm), beta(beta), proj(Projection::Parameter{25, 25, 0, 0}, 0.5){}

void RadialMachining::initialize_views(Visualization* viz){
    this->shadow_view = viz->add_view("Shadows", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::DENSITY);
    this->shadow_view_continuous = viz->add_view("Shadows Continuous", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::DENSITY);
    this->grad_view = viz->add_view("Shadows Gradient", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::DENSITY);
}

void RadialMachining::initialize(const Optimizer* const op){
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
        if(std::abs(vn.Dot(e.normal)) < 1.0-Precision::Confusion()){
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
    this->Hgrad.resize(elem_num,1);
    this->Phi = Eigen::SparseMatrix<double>(phi_size, phi_size);
}

double RadialMachining::calculate(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return 0;
    }

    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    if(!this->first_time){
        std::fill(this->Phi.valuePtr(), this->Phi.valuePtr() + this->Phi.nonZeros(), 0);
    }
    std::fill(this->b.begin(), this->b.end(), 0);
    std::fill(this->diff.begin(), this->diff.end(), 0);

    std::vector<double> v(3,0);
    double dv = -3;
    for(size_t i = 0; i < 3; ++i){
        const double ai = this->axis.Coord(1+i);
        dv += ai*ai;
    }

    const double BETA_RHO = 100;

    auto x_it = x.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            auto N = g->mesh.front()->helmholtz_vector(this->mesh->thickness);
            for(const auto& e:g->mesh){
                const gp_Pnt p = e->get_centroid();
                const gp_Vec pc(p, this->center);
                const gp_Vec vv = pc - pc.Dot(this->axis)*this->axis;
                v[0] = vv.X();
                v[1] = vv.Y();
                v[2] = vv.Z();

                const double Hx = this->heaviside(*x_it, BETA_RHO, 0.9);
                const auto phi_e = e->get_phi_radial(this->mesh->thickness, this->beta, this->v_norm, v, dv, Hx);
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
                        this->Phi.coeffRef(id1, id2) += phi_e[i*num_nodes + j];
                    }
                }

                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& ni = e->nodes[i];
                    const long id1 = this->id_mapping[ni->id];
                    if(id1 < 0){
                        continue;
                    }
                    this->b[id1] += this->beta*Hx*N[i];
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
                *d_it *= 1.0 - *x_it;
                ++d_it;
                x_it += num_den;
            }
        }
    }

    this->proj.project_densities(this->diff);

    return std::accumulate(this->diff.begin(), this->diff.end(), 0.0);// /this->diff.size();
}

double RadialMachining::calculate_with_gradient(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return 0;
    }

    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    if(!this->first_time){
        std::fill(this->Phi.valuePtr(), this->Phi.valuePtr() + this->Phi.nonZeros(), 0);
    }
    std::fill(this->b.begin(), this->b.end(), 0);
    std::fill(this->b_grad.begin(), this->b_grad.end(), 0);
    std::fill(this->diff.begin(), this->diff.end(), 0);

    std::vector<double> v(3,0);
    double dv = -3;
    for(size_t i = 0; i < 3; ++i){
        const double ai = this->axis.Coord(1+i);
        dv += ai*ai;
    }

    const double BETA_RHO = 100;

    auto x_it = x.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            auto N = g->mesh.front()->helmholtz_vector(this->mesh->thickness);
            for(const auto& e:g->mesh){
                const gp_Pnt p = e->get_centroid();
                const gp_Vec pc(p, this->center);
                const gp_Vec vv = pc - pc.Dot(this->axis)*this->axis;
                v[0] = vv.X();
                v[1] = vv.Y();
                v[2] = vv.Z();

                const double Hx = this->heaviside(*x_it, BETA_RHO, 0.9);
                const auto phi_e = e->get_phi_radial(this->mesh->thickness, this->beta, this->v_norm, v, dv, Hx);
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
                        this->Phi.coeffRef(id1, id2) += phi_e[i*num_nodes + j];
                    }
                }

                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& ni = e->nodes[i];
                    const long id1 = this->id_mapping[ni->id];
                    if(id1 < 0){
                        continue;
                    }
                    this->b[id1] += this->beta*Hx*N[i];
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
                *d_it *= 1.0 - *x_it;
                ++d_it;
                x_it += num_den;
            }
        }
    }

    this->shadow_view_continuous->update_view(this->diff);
    std::fill(this->Hgrad.begin(), this->Hgrad.end(), 1);
    this->proj.project_gradient(this->Hgrad, this->diff);
    this->proj.project_densities(this->diff);

    this->shadow_view->update_view(this->diff);

    const double result = std::accumulate(this->diff.begin(), this->diff.end(), 0.0);// /this->diff.size();

    auto dg_it = this->Hgrad.cbegin();
    x_it = x.cbegin();
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
                    b_grad[id1] += (*dg_it)*(1.0 - *x_it)/num_nodes;
                }
                ++dg_it;
                x_it += num_den;
            }
        }
    }

    Eigen::VectorXd psi_tilde = this->solver.solve(this->b_grad);

    x_it = x.cbegin();
    auto g_it = grad.begin();
    auto hg_it = this->Hgrad.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            auto phi_e = g->mesh.front()->get_phi_grad(this->mesh->thickness, this->beta);
            auto b_e = g->mesh.front()->helmholtz_vector(this->mesh->thickness);
            const double dHx = this->heaviside_grad(*x_it, BETA_RHO, 0.9);
            for(auto& i:b_e){
                i *= dHx*this->beta;
            }
            for(auto& p:phi_e){
                p *= dHx;
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
                    *g_it += psi[id1]/num_nodes;
                }
                *g_it *= -(*hg_it);

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
                        psi_tmp += phi_e[i*num_nodes + j]*psi[id2];
                    }
                    *g_it -= psi_tilde[id1]*psi_tmp;
                    *g_it += psi_tilde[id1]*b_e[i];
                }
                //*g_it /= this->diff.size();
                ++hg_it;
                g_it += num_den;
                x_it += num_den;
            }
        }
    }

    this->grad_view->update_view(grad);

    return result;
}

}
