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

#include "function/am_support.hpp"
#include "logger.hpp"
#include "projection/threshold.hpp"
#include <Eigen/src/Core/util/Constants.h>
#include <numeric>
#include <set>
#include <mpich-x86_64/mpi.h>

namespace function{

AMSupport::AMSupport(const Meshing* const mesh, const DensityFilter* const filter, const Projection* const global_proj, gp_Dir axis, double v_norm, double L, double beta, double angle)
    : mesh(mesh), filter(filter), global_proj(global_proj), axis(axis), v_norm(v_norm), L(L), beta(beta), support_angle(angle), proj(projection::Threshold::Parameter{100, 50, 0, 0}, 0.75){}

void AMSupport::initialize_views(Visualization* viz){
    this->shadow_view = viz->add_view("AM Supports", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::DENSITY);
}

void AMSupport::initialize(const Optimizer* const op){
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
        if(this->axis.Dot(e.normal) < -1e-7){
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
    this->gradx.resize(elem_num*N,0);
    this->proj_grad.resize(elem_num,1);
    this->Phi = Eigen::SparseMatrix<double>(phi_size, phi_size);
}

double AMSupport::calculate(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return 0;
    }

    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    size_t N = 0;
    if(this->mesh->elem_info->get_problem_type() == utils::PROBLEM_TYPE_2D){
        N = 2;
    } else if(this->mesh->elem_info->get_problem_type() == utils::PROBLEM_TYPE_3D){
        N = 3;
    }
    if(!this->first_time){
        std::fill(this->Phi.valuePtr(), this->Phi.valuePtr() + this->Phi.nonZeros(), 0);
    }
    std::fill(this->b.begin(), this->b.end(), 0);
    std::fill(this->b_grad.begin(), this->b_grad.end(), 0);
    this->filter->get_gradient(this->gradx);
    std::fill(this->proj_grad.begin(), this->proj_grad.end(), 1);
    this->global_proj->project_gradient(this->proj_grad, op->get_filtered_densities());
    for(size_t i = 0; i < this->proj_grad.size(); ++i){
        for(size_t j = 0; j < N; ++j){
            this->gradx[i*N + j] *= this->proj_grad[i];
        }
    }
    std::vector<double> v(3,0);
    v[0] = this->axis.X();
    v[1] = this->axis.Y();
    v[2] = this->axis.Z();
    const std::vector<double> A{std::sqrt(v[0]*v[0] + 1e-14), 0.0, 0.0,
                                0.0, std::sqrt(v[1]*v[1] + 1e-14), 0.0,
                                0.0, 0.0, std::sqrt(v[2]*v[2] + 1e-14)};
    auto x_it = x.cbegin();
    auto gx_it = this->gradx.cbegin();
    std::vector<double> p(num_nodes,0);
    const double cos_a = std::cos(this->support_angle);
    auto filter_mapping = this->filter->get_id_mapping();
    auto xn = this->filter->get_nodal_densities();
    size_t geom_id = 0;
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                double beta_switch = 0;
                double norm = 0;

                for(size_t i = 0; i < N; ++i){
                    beta_switch += v[i]*(*gx_it);
                    norm += (*gx_it)*(*gx_it);
                    ++gx_it;
                }
                norm = std::sqrt(norm + 1e-14);
                beta_switch = this->H3(beta_switch/norm - cos_a, 1000, 0.0);

                const double Hx = this->H3(*x_it, 1000, 0.9);
                const auto Ne = (1.0 - Hx)*this->beta*beta_switch*e->source_1dof(this->mesh->thickness);
                const auto psi_e = (1.0 - Hx)*(this->beta*beta_switch*e->absorption_1dof(this->mesh->thickness) +
                                   this->L*this->L*e->diffusion_1dof(this->mesh->thickness, A) +
                                   this->v_norm*this->L*e->advection_1dof(this->mesh->thickness, v).transpose());

                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& ni = e->nodes[i];
                    const long id1 = this->id_mapping[ni->id];
                    if(id1 < 0){
                        continue;
                    }
                    this->b[id1] += Ne[i];
                }

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
                        this->Phi.coeffRef(id1, id2) += psi_e(i, j);
                    }
                }

                x_it += num_den;
            }
        }
        ++geom_id;
    }

    if(this->first_time){
        this->Phi.makeCompressed();
        this->first_time = false;
        this->solver.analyzePattern(Phi);
    }

    this->solver.factorize(this->Phi);
    logger::log_assert(this->solver.info() == Eigen::Success, logger::ERROR, "matrix decomposition failed");
    Eigen::VectorXd psi = this->solver.solve(this->b);

    auto d_it = this->diff.begin();
    geom_id = 0;
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
                    *d_it += psi[id1];
                    *d_it -= xn[filter_mapping[std::make_pair(geom_id, ni->id)]];
                }
                *d_it /= num_nodes;
                ++d_it;
            }
        }
        ++geom_id;
    }

    std::fill(this->Hgrad.begin(), this->Hgrad.end(), 1);
    this->proj.project_gradient(this->Hgrad, this->diff);
    this->proj.project_densities(this->diff);

    this->shadow_view->update_view(this->diff);

    //const double result = std::accumulate(this->diff.begin(), this->diff.end(), 0);
    double result = 0;
    for(auto& d:this->diff){
        result += d;
    }
    return result;
}

double AMSupport::calculate_with_gradient_nodal(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id != 0){
        return 0;
    }

    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    size_t N = 0;
    if(this->mesh->elem_info->get_problem_type() == utils::PROBLEM_TYPE_2D){
        N = 2;
    } else if(this->mesh->elem_info->get_problem_type() == utils::PROBLEM_TYPE_3D){
        N = 3;
    }
    if(!this->first_time){
        std::fill(this->Phi.valuePtr(), this->Phi.valuePtr() + this->Phi.nonZeros(), 0);
    }
    std::fill(this->b.begin(), this->b.end(), 0);
    std::fill(this->b_grad.begin(), this->b_grad.end(), 0);
    this->filter->get_gradient(this->gradx);
    std::fill(this->proj_grad.begin(), this->proj_grad.end(), 1);
    this->global_proj->project_gradient(this->proj_grad, op->get_filtered_densities());
    for(size_t i = 0; i < this->proj_grad.size(); ++i){
        for(size_t j = 0; j < N; ++j){
            this->gradx[i*N + j] *= this->proj_grad[i];
        }
    }
    std::vector<double> v(3,0);
    v[0] = this->axis.X();
    v[1] = this->axis.Y();
    v[2] = this->axis.Z();
    const std::vector<double> A{std::sqrt(v[0]*v[0] + 1e-14), 0.0, 0.0,
                                0.0, std::sqrt(v[1]*v[1] + 1e-14), 0.0,
                                0.0, 0.0, std::sqrt(v[2]*v[2] + 1e-14)};
    auto x_it = x.cbegin();
    auto gx_it = this->gradx.cbegin();
    std::vector<double> p(num_nodes,0);
    const double cos_a = std::cos(this->support_angle);
    auto filter_mapping = this->filter->get_id_mapping();
    auto xn = this->filter->get_nodal_densities();
    size_t geom_id = 0;
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                double beta_switch = 0;
                double norm = 0;

                for(size_t i = 0; i < N; ++i){
                    // Opposite of gradient is the boundary normal
                    beta_switch -= v[i]*(*gx_it);
                    norm += (*gx_it)*(*gx_it);
                    ++gx_it;
                }
                // TODO: fix this
                // It is not generating for the belowmost part of the geometry,
                // for some reason.
                // After fixing it, update gradient calculation and remove
                // std::fill at the end of the function.
                norm = std::sqrt(norm + 1e-14);
                beta_switch = this->H3(beta_switch/norm - cos_a, 1000, 0.0);

                const double Hx = this->H3(*x_it, 1000, 0.9);
                const auto Ne = (1.0 - Hx)*this->beta*beta_switch*e->source_1dof(this->mesh->thickness);
                const auto psi_e = (1.0 - Hx)*(this->beta*beta_switch*e->absorption_1dof(this->mesh->thickness) +
                                   this->L*this->L*e->diffusion_1dof(this->mesh->thickness, A) +
                                   this->v_norm*this->L*e->advection_1dof(this->mesh->thickness, v).transpose());
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& ni = e->nodes[i];
                    const long id1 = this->id_mapping[ni->id];
                    if(id1 < 0){
                        continue;
                    }
                    this->b[id1] += Ne[i];
                }

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
                        this->Phi.coeffRef(id1, id2) += psi_e(i, j);
                    }
                }

                x_it += num_den;
            }
        }
        ++geom_id;
    }
    for(long i = 0; i < this->b.size(); ++i){
        this->Phi.coeffRef(i, i) += 1e-6;
    }

    if(this->first_time){
        this->Phi.makeCompressed();
        this->first_time = false;
        this->solver.analyzePattern(Phi);
    }

    this->solver.factorize(this->Phi);
    logger::log_assert(this->solver.info() == Eigen::Success, logger::ERROR, "matrix decomposition failed");
    Eigen::VectorXd psi = this->solver.solve(this->b);

    auto d_it = this->diff.begin();
    geom_id = 0;
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
                    *d_it += psi[id1];
                    *d_it -= xn[filter_mapping[std::make_pair(geom_id, ni->id)]];
                }
                *d_it /= num_nodes;
                ++d_it;
            }
        }
        ++geom_id;
    }

    std::fill(this->Hgrad.begin(), this->Hgrad.end(), 1);
    this->proj.project_gradient(this->Hgrad, this->diff);
    this->proj.project_densities(this->diff);

    this->shadow_view->update_view(this->diff);

    //const double result = std::accumulate(this->diff.begin(), this->diff.end(), 0);
    double result = 0;
    for(auto& d:this->diff){
        result += d;
    }

    auto dg_it = this->Hgrad.cbegin();
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

    std::fill(grad.begin(), grad.end(), 0);

    auto g_it = grad.begin();
    auto hg_it = this->Hgrad.cbegin();
    x_it = x.cbegin();
    gx_it = this->gradx.cbegin();
    geom_id = 0;
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                double psi_tilde_sum = 0;
                auto Ne = e->get_nodal_density_gradient(e->get_centroid());
                double Ne_sum = e->get_volume(this->mesh->thickness);
                //auto Me = e->get_phi_grad(this->mesh->thickness, 1.0);
                const auto Me = e->absorption_1dof(this->mesh->thickness);
                auto Me_dot = 0;
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& ni = e->nodes[i];
                    const long id1 = this->id_mapping[ni->id];
                    if(id1 < 0){
                        continue;
                    }
                    double Me_psi = 0;
                    for(size_t j = 0; j < num_nodes; ++j){
                        const auto& nj = e->nodes[j];
                        const long id2 = this->id_mapping[nj->id];
                        if(id2 < 0){
                            continue;
                        }
                        Me_psi += Me(i, j)*psi[id2];
                    }
                    Me_dot += psi_tilde[id1]*Me_psi;
                }
                const auto gradN = e->get_nodal_density_gradient(e->get_centroid());
                double h1 = 0, h2 = 0, h3 = 0, norm = 0, dh1 = 0, dh2 = 0, dh3 = 0;
                for(size_t i = 0; i < num_nodes; ++i){
                    grad[filter_mapping[std::make_pair(geom_id, e->nodes[i]->id)]] -= *hg_it/num_nodes;
                    p[i] = xn[filter_mapping[std::make_pair(geom_id, e->nodes[i]->id)]];
                    psi_tilde_sum += psi_tilde[this->id_mapping[e->nodes[i]->id]];
                }
                for(size_t i = 0; i < N; ++i){
                    h2 += v[i]*(*(gx_it+i));
                    norm += (*(gx_it+i))*(*(gx_it+i));
                }
                norm = std::sqrt(norm + 1e-14);
                h1 = this->H3(norm, 1000, 0.1);
                dh1 = this->dH3(norm, 1000, 0.1);
                h3 = this->H3(h2, 1000, 0.02);
                dh3 = this->dH3(h2, 1000, 0.02);
                //h2 = this->H(cos_a - h2/norm, 50, 0, 1.0, 0.0);
                //dh2 = this->dH(cos_a - h2/norm, 50, 0, 1.0);
                h2 = this->H3(cos_a - h2/norm, 1000, 0);
                dh2 = this->dH3(cos_a - h2/norm, 1000, 0);
                double h2l = this->H2(h2/norm - cos_a, 1000, 0.0);
                double dh2l = this->dH2(h2/norm - cos_a, 1000, 0.0);
                for(size_t i = 0; i < num_nodes; ++i){
                    double switch_deriv = 0;
                    double norm_deriv = 0;
                    double vDp_deriv = 0;
                    for(size_t j = 0; j < N; ++j){
                        norm_deriv += Ne[j*num_nodes + i]*Ne[j*num_nodes + i];
                        switch_deriv += v[j]*(Ne[j*num_nodes + i]/norm - norm_deriv*(*(gx_it + j))/(norm*norm));
                        vDp_deriv += v[j]*Ne[j*num_nodes + i];
                    }
                    norm_deriv *= p[i]/norm;
                    //const double mult = 1.0/num_nodes - (dh1*norm_deriv*h2l + h1*dh2l*(-switch_deriv));
                    //const double mult = this->beta*(h1*h2*h3/num_nodes + (*x_it)*(dh1*norm_deriv*h2*h3 + h1*dh2*switch_deriv*h3 + h1*h2*dh3*vDp_deriv));
                    //const double mult = this->beta*((dh1*norm_deriv*(h2 -h3) + h1*(dh2*switch_deriv*h3 + h2*dh3*vDp_deriv)));
                    const double mult = this->beta*((dh1*norm_deriv*h2*h3 + h1*dh2*(-switch_deriv)*h3 + h1*h2*dh3*vDp_deriv));
                    for(size_t j = 0; j < num_nodes; ++j){
                        grad[filter_mapping[std::make_pair(geom_id, e->nodes[j]->id)]] -= mult*Ne_sum*psi_tilde_sum;
                        grad[filter_mapping[std::make_pair(geom_id, e->nodes[j]->id)]] += mult*Me_dot;
                    }
                }
                gx_it += N;
                ++hg_it;
                x_it += num_den;
                g_it += num_den;
            }
        }
        ++geom_id;
    }

    std::fill(grad.begin(), grad.end(), 0);

    return result;
}

}
