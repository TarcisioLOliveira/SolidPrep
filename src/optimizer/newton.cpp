/*
 *   Copyright (C) 2024 Tarcísio Ladeia de Oliveira.
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

#include "optimizer/newton.hpp"
#include "logger.hpp"
#include "project_data.hpp"
#include "optimization/MMASolver.hpp"
#include <cblas.h>
#include <chrono>
#include <mpich-x86_64/mpi.h>

namespace optimizer{

Newton::Newton(DensityFilter* filter, Projection* projection, ProjectData* data, std::vector<Constraint> functions, double pc, double psi, double rho_init, double xtol_abs, double ftol_rel, double result_threshold, bool save):
    data(data), rho_init(rho_init), xtol_abs(xtol_abs), ftol_rel(ftol_rel), pc(pc), psi(psi), result_threshold(result_threshold), save_result(save), functions(std::move(functions)), filter(filter), projection(projection), viz(nullptr)
    {}

void Newton::initialize_views(Visualization* viz){
    this->viz = viz;

    this->stress_view = viz->add_view("Von Mises Stress", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
    this->density_view = viz->add_view("Elemental Density", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::DENSITY);

    for(auto& f:this->functions){
        f.fun->initialize_views(viz);
    }
}

TopoDS_Shape Newton::optimize(SolverManager* fem, Meshing* mesh){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id == 0){
        logger::quick_log("Preparing for optimization...");
    }

    if(mpi_id == 0){
        this->initialize_optimizer(mesh);
    }
    auto stress_render = this->stresses;

    size_t x_size = 0;
    size_t x_size_view = 0;
    if(mpi_id == 0){
        for(auto& f:this->functions){
            f.fun->initialize(this);
        }

        for(const auto& g:mesh->geometries){
            if(g->do_topopt){
                const size_t num_den = g->number_of_densities_needed();
                x_size += g->mesh.size()*num_den;
                x_size_view += g->mesh.size();
            }
        }
    }

    auto start_to = std::chrono::high_resolution_clock::now();

    std::vector<double> x(x_size, this->rho_init);
    this->filtered_densities.resize(x_size);
    std::fill(this->filtered_densities.begin(), this->filtered_densities.end(), this->rho_init);
    std::vector<double>& x_fil = this->filtered_densities;
    std::vector<double> new_x(x_size, this->rho_init);
    std::vector<double> x_view(x_size_view, this->rho_init);

    if(mpi_id == 0){
        this->filter->initialize(mesh, x_size);
    }

    size_t M = 2*this->functions.size();

    std::vector<double> dgtmp1(x.size());
    std::vector<double> dgtmp2(x.size());
    std::vector<double> dgtmp1_nodal(this->filter->get_nodal_density_size());
    std::vector<double> dg(x.size()*M);
    std::vector<double> df(x.size(),0);
    std::vector<double> g(M);

    // std::vector<double> xnew(x);
    std::vector<double> xold(x);

    std::vector<double> dx(x.size());
    std::vector<double> dxold(x.size());
    //std::vector<double> Jacobi(x.size());
    //std::vector<double> r0(x.size());
    //std::vector<double> r1(x.size());
    //std::vector<double> p(x.size());
    //std::vector<double> JTg(x.size());
    
    std::vector<double> xmin = std::vector<double>(x.size(), -0.1);
    std::vector<double> xmax = std::vector<double>(x.size(), 0.1);

    if(mpi_id == 0){
        logger::quick_log("Done.");
    }

    if(mpi_id == 0){
        this->projection->project_densities(new_x);
    }
    auto loads = mesh->load_vector;
    std::vector<double> u(mesh->max_dofs, 0);

    // int max_innerit = 30;
    double ch = 1.0;
    int iter;
	for (iter = 0; (ch > this->xtol_abs); ++iter){

        fem->generate_matrix(mesh, new_x, pc, this->psi);
        for(size_t i = 0; i < mesh->node_positions.size(); ++i){
            auto& lv = mesh->load_vector[i];
            auto& l = loads[i];
            std::copy(lv.begin(), lv.end(), l.begin());
        }
        fem->calculate_displacements_global(mesh, loads, u);
        if(mpi_id == 0){
            this->get_stresses(mesh->geometries, u, fem->D_vec, this->stresses);
            std::copy(stresses.begin(), stresses.end(), stress_render.begin());
            auto xit = new_x.cbegin();
            auto xvit = x_view.begin();
            for(const auto& g:mesh->geometries){
                if(g->do_topopt){
                    auto xvit2 = xvit;
                    xvit += g->mesh.size();
                    const size_t num_den = g->number_of_densities_needed();
                    while(xvit2 < xvit){
                        *xvit2 = *xit;
                        xit += num_den;
                        ++xvit2;
                    }
                }
            }

            this->density_view->update_view(x_view);
            this->stress_view->update_view(stress_render);
        }
        size_t g_id = 0;
        for(size_t i = 0; i < this->functions.size(); ++i){
            auto& c = this->functions[i];
            if(c.fun->filter_gradient_type() == DensityFilter::FilterGradient::ELEMENTAL){
                double val = c.fun->calculate_with_gradient(this, u, new_x, dgtmp1);
                for(size_t k = 0; k < c.types.size(); ++k){
                    g[2*g_id + 0] =   val - c.bounds[k];
                    g[2*g_id + 1] = -(val - c.bounds[k]);
                    if(mpi_id == 0){
                        this->projection->project_gradient(dgtmp1, x_fil);
                        this->filter->filter_gradient(dgtmp1, dgtmp2);
                        for(size_t j = 0; j < x_size; ++j){
                            dg[M*j + 2*g_id+0] =  dgtmp2[j];
                            dg[M*j + 2*g_id+1] = -dgtmp2[j];
                        }
                    }
                    ++g_id;
                }
            } else {
                double val = c.fun->calculate_with_gradient_nodal(this, u, new_x, dgtmp1_nodal);
                for(size_t k = 0; k < c.types.size(); ++k){
                    g[2*g_id + 0] =   val - c.bounds[k];
                    g[2*g_id + 1] = -(val - c.bounds[k]);
                    if(mpi_id == 0){
                        //this->projection->project_gradient(dgtmp1, x_fil);
                        this->filter->filter_gradient_nodal(dgtmp1_nodal, dgtmp2);
                        for(size_t j = 0; j < x_size; ++j){
                            dg[M*j + 2*g_id+0] =  dgtmp2[j];
                            dg[M*j + 2*g_id+1] = -dgtmp2[j];
                        }
                    }
                    ++g_id;
                }
            }
        }

        if(mpi_id == 0){
            /* Crude attempt using PCG
            const double EPS = 1e-7;
            double max_diff = 1.0;
            std::fill(dx.begin(), dx.end(), 0);
            std::fill(r1.begin(), r1.end(), 0);
            double alpha = 0;
            double beta = 0;
            double beta_tmp = 0;
            double alpha_div = 0;
            // Matrix-free Newton method, because the generalized inverse
            // can get quite large
            logger::quick_log("Updating density...");
            #pragma omp parallel
            {
                #pragma omp for
                for(size_t i = 0; i < x.size(); ++i){
                    Jacobi[i] = 0;
                    JTg[i] = 0;
                    for(size_t k = 0; k < M; ++k){
                        Jacobi[i] += dg[M*i + k]*dg[M*i + k];
                    }
                    Jacobi[i] *= Jacobi[i];
                    for(size_t k = 0; k < M; ++k){
                        JTg[i] -= dg[M*i + k]*g[k];
                    }
                }
                #pragma omp for
                for(size_t i = 0; i < x.size(); ++i){
                    r0[i] = JTg[i];
                    for(size_t j = 0; j < x.size(); ++j){
                        double tmp = 0;
                        for(size_t k = 0; k < M; ++k){
                            tmp += dg[M*i + k]*dg[M*j + k];
                        }
                        r0[i] -= tmp*dx[j];
                    }
                    p[i] = r0[i]/Jacobi[i];
                }
            }
            logger::quick_log("Starting optimization");
            while(max_diff > EPS){
                alpha = 0;
                beta = 0;
                beta_tmp = 0;
                alpha_div = 0;
                max_diff = 0;
                #pragma omp parallel
                {
                    #pragma omp for reduction(+:alpha,alpha_div)
                    for(size_t i = 0; i < x.size(); ++i){
                        alpha += r0[i]*r0[i]/Jacobi[i];
                        double alpha_tmp = 0;
                        for(size_t j = 0; j < x.size(); ++j){
                            double tmp = 0;
                            for(size_t k = 0; k < M; ++k){
                                tmp += dg[M*i + k]*dg[M*j + k];
                            }
                            alpha_tmp += tmp*p[j];
                        }
                        alpha_div += alpha_tmp*p[i];
                    }
                    #pragma omp single
                    {
                        beta_tmp = alpha;
                        alpha /= alpha_div;
                    }
                    #pragma omp for reduction(+:beta)
                    for(size_t i = 0; i < x.size(); ++i){
                        dx[i] += alpha*p[i];
                        r1[i] = JTg[i];
                        for(size_t j = 0; j < x.size(); ++j){
                            double tmp = 0;
                            for(size_t k = 0; k < M; ++k){
                                tmp += dg[M*i + k]*dg[M*j + k];
                            }
                            r1[i] -= tmp*dx[j];
                        }
                        beta += r1[i]*r1[i]/Jacobi[i];
                    }
                    #pragma omp single
                    {
                        beta /= beta_tmp;
                    }
                    #pragma omp for reduction(max:max_diff)
                    for(size_t i = 0; i < x.size(); ++i){
                        p[i] = r1[i]/Jacobi[i] + beta*p[i];
                        const double diff = std::abs(r1[i]);
                        if(diff > max_diff){
                            max_diff = diff;
                        }
                        r0[i] = r1[i];
                    }
                }
                logger::quick_log(max_diff);
            }
            */

            // Update x
            logger::quick_log("Updating density...");
            optimization::MMASolver mma(x.size(), M, 0, 1e-15, 1);
            mma.SetAsymptotes(0.03, 0.3, 1.5);
            mma.SetBoundFactors(1e-5, 0.1);

            std::fill(dx.begin(), dx.end(), 0);
            std::fill(dxold.begin(), dxold.end(), 0);

            const double EPS = 1e-3;
            double max_diff = 1.0;
            while(max_diff > EPS){
                max_diff = 0;
                mma.Update(dx.data(), df.data(), g.data(), dg.data(), xmin.data(), xmax.data());
                for(size_t i = 0; i < dx.size(); ++i){
                    const double diff = std::abs(dx[i] - dxold[i]);
                    if(diff > max_diff){
                        max_diff = diff;
                    }
                    dxold[i] = dx[i];
                }
                logger::quick_log(max_diff);
            }
            #pragma omp parallel for
            for(size_t i = 0; i < x.size(); ++i){
                x[i] = std::max(std::min(x[i] + dx[i], 1.0), 0.0);
            }
            logger::quick_log("Done");

            // Update functions
            for(auto& f:this->functions){
                f.fun->update();
            }
            this->projection->update(iter);

            ch = 0.0;
            for (size_t i = 0; i < x.size(); ++i) {
                ch = std::max(ch, std::abs(xold[i] - x[i]));
                xold[i] = x[i];
            }

            // Filter densities for new iteration
            this->filter->filter_densities(x, x_fil);
            std::copy(x_fil.begin(), x_fil.end(), new_x.begin());
            this->projection->project_densities(new_x);
        }
        MPI_Bcast(&ch, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(mpi_id == 0){
            logger::quick_log("");
            logger::quick_log("");
            logger::quick_log("Iteration: ", iter);
            logger::quick_log("Results: ");
            logger::quick_log(g);
            logger::quick_log("");
            logger::quick_log("Design var change: ", ch);
            logger::quick_log("");
        }
	}

    if(mpi_id == 0){
        logger::quick_log("");
        auto stop_to = std::chrono::high_resolution_clock::now();
        logger::quick_log("Final result: ");
        logger::quick_log(g);
        auto to_duration = std::chrono::duration_cast<std::chrono::seconds>(stop_to-start_to);
        double to_time = to_duration.count();
        double it_time = to_time/iter;
        logger::quick_log("Time (topology optimization): ", to_time, " seconds");
        logger::quick_log("Time per iteration (topology optimization): ", it_time, " seconds");
        logger::quick_log("Number of iterations (topology optimization): ", iter);

        if(this->save_result){
            return this->make_shape(new_x, mesh->geometries, this->result_threshold, this->data->type);
        }
    }

    return TopoDS_Shape();
}

}
