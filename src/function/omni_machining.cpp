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

#include "function/omni_machining.hpp"
#include "logger.hpp"
#include "projection/heaviside.hpp"
#include "projection/threshold.hpp"
#include <Eigen/src/Core/util/Constants.h>
#include <limits>
#include <mpi.h>
#include <numeric>
#include <set>

namespace function{

OmniMachining::OmniMachining(const Meshing* const mesh, const DensityFilter* const filter, gp_Pnt center, gp_Dir axis, double v_norm, double beta1, double beta2, double L)
        : mesh(mesh), filter(filter), center(center), axis(axis), v_norm(v_norm), beta1(beta1), beta2(beta2), L(L), proj(Projection::Parameter{20, 20, 0, 0}, 0.5){}

void OmniMachining::initialize_views(Visualization* viz){
    this->shadow_view = viz->add_view("Shadows", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::DENSITY);
    this->shadow_view_continuous = viz->add_view("Shadows Continuous", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::DENSITY);
    this->grad_view = viz->add_view("Shadows Gradient", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::DENSITY);
}

void OmniMachining::initialize(const Optimizer* const op){
    (void)op;
    const size_t num_nodes = mesh->elem_info->get_nodes_per_element();
    const size_t num_nodes_bound = mesh->elem_info->get_boundary_nodes_per_element();
    std::map<size_t, long> id_mapping;
    size_t elem_num = 0;
    size_t phi_size = 0;
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& n = e->nodes[i];
                    id_mapping[n->id] = 0;
                }
            }
            elem_num += g->mesh.size();
        }
    }
    for(auto& e:mesh->boundary_elements){
        if(std::abs(this->axis.Dot(e.normal)) < 1.0-Precision::Confusion()){
            for(size_t i = 0; i < num_nodes_bound; ++i){
                if(id_mapping.count(e.nodes[i]->id)){
                    id_mapping[e.nodes[i]->id] = -1;
                }
            }
        }
    }
    size_t id = 0;
    for(auto& n:this->mesh->node_list){
        if(id_mapping.count(n->id) && id_mapping.at(n->id) > -1){
            id_mapping[n->id] = id;
            ++id;
        }
    }
    this->id_mapping_linear.resize(elem_num*num_nodes,0);
    auto id_it = this->id_mapping_linear.begin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            for(const auto& e:g->mesh){
                for(size_t i = 0; i < num_nodes; ++i){
                    const auto& ni = e->nodes[i];
                    const long id1 = id_mapping[ni->id];
                    *id_it = id1;
                    ++id_it;
                }
            }
        }
    }
    phi_size = id;
    this->b_grad.resize(phi_size);
    this->b.resize(phi_size);
    this->diff.resize(elem_num,0);
    this->Hgrad.resize(elem_num,1);
    this->Phi = Eigen::SparseMatrix<double>(phi_size, phi_size);
}

double OmniMachining::calculate(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x){
    (void)op;
    (void)u;
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

    std::vector<double> A{axis.X(), axis.Y(), axis.Z()};
    std::vector<double> C{center.X(), center.Y(), center.Z()};

    const double BETA_RHO = 75;
    const double HX_ETA = 0.9;

    auto x_it = x.cbegin();
    auto id_it = this->id_mapping_linear.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                const double Hx = this->heaviside(*x_it, BETA_RHO, HX_ETA);
                const auto p = e->get_centroid();
                std::vector<double> v{C[0] - p.X(), C[1] - p.Y(), C[2] - p.Z()};
                const double vnorm = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + 1e-14);
                v[0] = v[0]/vnorm;
                v[1] = v[1]/vnorm;
                v[2] = v[2]/vnorm;
                const std::vector<double> a{std::sqrt(v[0]*v[0] + 1e-14), 0.0, 0.0,
                                            0.0, std::sqrt(v[1]*v[1] + 1e-14), 0.0,
                                            0.0, 0.0, std::sqrt(v[2]*v[2] + 1e-14)};

                const auto N = Hx*beta2*e->source_1dof(this->mesh->thickness);
                const auto phi_e = (beta1*(1.0-Hx) + beta2*Hx)*e->absorption_1dof(this->mesh->thickness) +
                                   L*L*e->diffusion_1dof(this->mesh->thickness, a) +
                                   v_norm*L*e->advection_1dof(this->mesh->thickness, v).transpose();
                for(size_t i = 0; i < num_nodes; ++i){
                    const long id1 = *(id_it + i);
                    if(id1 < 0){
                        continue;
                    }
                    for(size_t j = 0; j < num_nodes; ++j){
                        const long id2 = *(id_it + j);
                        if(id2 < 0){
                            continue;
                        }
                        this->Phi.coeffRef(id1, id2) += phi_e(i, j);
                    }
                }

                for(size_t i = 0; i < num_nodes; ++i){
                    const long id1 = *(id_it + i);
                    if(id1 < 0){
                        continue;
                    }
                    this->b[id1] += N[i];
                }

                x_it += num_den;
                id_it += num_nodes;
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
    id_it = this->id_mapping_linear.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(size_t it = 0; it < g->mesh.size(); ++it){
                for(size_t i = 0; i < num_nodes; ++i){
                    const long id1 = *(id_it + i);
                    if(id1 < 0){
                        continue;
                    }
                    *d_it += psi[id1]/num_nodes;
                }
                *d_it *= 1.0 - *x_it;
                ++d_it;
                x_it += num_den;
                id_it += num_nodes;
            }
        }
    }

    this->proj.project_densities(this->diff);

    return std::accumulate(this->diff.begin(), this->diff.end(), 0.0);
}

double OmniMachining::calculate_with_gradient(const Optimizer* const op, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& grad){
    (void)op;
    (void)u;
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

    const double BETA_RHO = 50;
    const double HX_ETA = 0.9;

    std::vector<double> A{axis.X(), axis.Y(), axis.Z()};
    std::vector<double> C{center.X(), center.Y(), center.Z()};

    auto x_it = x.cbegin();
    auto id_it = this->id_mapping_linear.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                const double Hx = this->heaviside(*x_it, BETA_RHO, HX_ETA);
                const auto p = e->get_centroid();
                std::vector<double> v{C[0] - p.X(), C[1] - p.Y(), C[2] - p.Z()};
                const double vnorm = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + 1e-14);
                v[0] = v[0]/vnorm;
                v[1] = v[1]/vnorm;
                v[2] = v[2]/vnorm;
                const std::vector<double> a{std::sqrt(v[0]*v[0] + 1e-14), 0.0, 0.0,
                                            0.0, std::sqrt(v[1]*v[1] + 1e-14), 0.0,
                                            0.0, 0.0, std::sqrt(v[2]*v[2] + 1e-14)};

                const auto N = Hx*beta2*e->source_1dof(this->mesh->thickness);
                const auto phi_e = (beta1*(1.0-Hx) + beta2*Hx)*e->absorption_1dof(this->mesh->thickness) +
                                   L*L*e->diffusion_1dof(this->mesh->thickness, a) +
                                   v_norm*L*e->advection_1dof(this->mesh->thickness, v).transpose();
                for(size_t i = 0; i < num_nodes; ++i){
                    const long id1 = *(id_it + i);
                    if(id1 < 0){
                        continue;
                    }
                    for(size_t j = 0; j < num_nodes; ++j){
                        const long id2 = *(id_it + j);
                        if(id2 < 0){
                            continue;
                        }
                        this->Phi.coeffRef(id1, id2) += phi_e(i, j);
                    }
                }

                for(size_t i = 0; i < num_nodes; ++i){
                    const long id1 = *(id_it + i);
                    if(id1 < 0){
                        continue;
                    }
                    this->b[id1] += N[i];
                }

                x_it += num_den;
                id_it += num_nodes;
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
    id_it = this->id_mapping_linear.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(size_t it = 0; it < g->mesh.size(); ++it){
                for(size_t i = 0; i < num_nodes; ++i){
                    const long id1 = *(id_it + i);
                    if(id1 < 0){
                        continue;
                    }
                    *d_it += psi[id1]/num_nodes;
                }
                *d_it *= 1.0 - *x_it;
                ++d_it;
                x_it += num_den;
                id_it += num_nodes;
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
    id_it = this->id_mapping_linear.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(size_t it = 0; it < g->mesh.size(); ++it){
                for(size_t i = 0; i < num_nodes; ++i){
                    const long id1 = *(id_it + i);
                    if(id1 < 0){
                        continue;
                    }
                    b_grad[id1] += (*dg_it)*(1.0 - *x_it)/num_nodes;
                }
                ++dg_it;
                x_it += num_den;
                id_it += num_nodes;
            }
        }
    }

    Eigen::VectorXd psi_tilde = this->solver.adjoint().solve(this->b_grad);

    x_it = x.cbegin();
    auto g_it = grad.begin();
    auto hg_it = this->Hgrad.cbegin();
    id_it = this->id_mapping_linear.cbegin();
    for(const auto& g:mesh->geometries){
        if(g->do_topopt){
            const size_t num_den = g->number_of_densities_needed();
            for(const auto& e:g->mesh){
                const double dHx = this->heaviside_grad(*x_it, BETA_RHO, HX_ETA);
                const auto b_e = dHx*beta2*e->source_1dof(this->mesh->thickness);
                const auto phi_e = ((beta1*(1.0-dHx)/10 + beta2*dHx)*e->absorption_1dof(this->mesh->thickness));
                *g_it = 0;
                double psi_tmp = 0;
                for(size_t i = 0; i < num_nodes; ++i){
                    const long id1 = *(id_it + i);
                    if(id1 < 0){
                        continue;
                    }
                    *g_it += psi[id1]/num_nodes;
                }
                *g_it *= -(*hg_it);

                for(size_t i = 0; i < num_nodes; ++i){
                    const long id1 = *(id_it + i);
                    if(id1 < 0){
                        continue;
                    }
                    psi_tmp = 0;
                    for(size_t j = 0; j < num_nodes; ++j){
                        const long id2 = *(id_it + j);
                        if(id2 < 0){
                            continue;
                        }
                        psi_tmp += phi_e(i,j)*psi[id2];
                    }
                    *g_it -= psi_tilde[id1]*psi_tmp;
                    *g_it += psi_tilde[id1]*b_e[i];
                }
                //*g_it /= this->diff.size();
                ++hg_it;
                g_it += num_den;
                x_it += num_den;
                id_it += num_nodes;
            }
        }
    }

    this->grad_view->update_view(grad);

    return result;
}

}
