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

#include "optimizer/mma.hpp"
#include "logger.hpp"
#include "optimization/MMASolver.hpp"
#include "project_data.hpp"
#include <chrono>
#include <mpich-x86_64/mpi.h>

namespace optimizer{

MMA::MMA(DensityFilter* filter, Projection* projection, ProjectData* data, std::vector<std::unique_ptr<DensityBasedFunction>> objective, std::vector<double> weights, std::vector<Constraint> constraints, double pc, double psi, double rho_init, double xtol_abs, double ftol_rel, double result_threshold, bool save):
    data(data), rho_init(rho_init), xtol_abs(xtol_abs), ftol_rel(ftol_rel), pc(pc), psi(psi), result_threshold(result_threshold), save_result(save), objective(std::move(objective)), objective_weights(std::move(weights)), constraints(std::move(constraints)), filter(filter), projection(projection), viz(nullptr)
    {}

void MMA::initialize_views(Visualization* viz){
    this->viz = viz;

    this->stress_view = viz->add_view("Von Mises Stress", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::STRESS);
    this->density_view = viz->add_view("Elemental Density", spview::defs::ViewType::ELEMENTAL, spview::defs::DataType::DENSITY);

    for(auto& f:this->objective){
        f->initialize_views(viz);
    }
    for(auto& f:this->constraints){
        f.fun->initialize_views(viz);
    }
}

TopoDS_Shape MMA::optimize(FiniteElement* fem, Meshing* mesh){
    int mpi_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    if(mpi_id == 0){
        logger::quick_log("Preparing for optimization...");
    }

    size_t fem_steps = 1;
    for(auto& f:this->objective){
        fem_steps += f->additional_steps();
    }
    for(auto& f:this->constraints){
        fem_steps += f.fun->additional_steps();
    }

    fem->set_steps(fem_steps);

    if(mpi_id == 0){
        this->initialize_optimizer(mesh);
    }
    auto stress_render = this->stresses;

    size_t x_size = 0;
    size_t x_size_view = 0;
    if(mpi_id == 0){
        for(auto& f:this->objective){
            f->initialize(this);
        }
        for(auto& f:this->constraints){
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

    size_t M = 0;
    for(auto& f:this->constraints){
        for(auto& t:f.types){
            if(t == Constraint::Type::GREATER_THAN || t == Constraint::Type::LESS_THAN){
                ++M;
            } else if(t == Constraint::Type::EQUAL){
                M += 2;
            }
        }
    }

    optimization::MMASolver mma(x_size, M, 0, 1e5, 1); //1e5
    mma.SetAsymptotes(0.005, 0.7, 1.2);

    double ff = 0;
    std::vector<double> dftmp(x.size());
    std::vector<double> dftmp_nodal(this->filter->get_nodal_density_size());
    std::vector<double> dgtmp1(x.size());
    std::vector<double> dgtmp2(x.size());
    std::vector<double> dgtmp1_nodal(this->filter->get_nodal_density_size());
    std::vector<double> df(x.size());
    std::vector<double> dg(x.size()*M);
    std::vector<double> g(M);

    // std::vector<double> xnew(x);
    std::vector<double> xold(x);

    std::vector<double> xmin;
    std::vector<double> xmax;

    xmin = std::vector<double>(x.size(), 0.0);
    xmax = std::vector<double>(x.size(), 1.0);

    if(mpi_id == 0){
        logger::quick_log("Done.");
    }

    double fnew = 0;

    if(mpi_id == 0){
        this->projection->project_densities(new_x);
    }
    auto u = fem->calculate_displacements(mesh, mesh->load_vector, new_x, this->pc);
    if(mpi_id == 0){
        this->get_stresses(mesh->geometries, u, new_x, this->stresses, this->pc, this->psi);
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
    if(this->objective.size() == 1){
        if(this->objective[0]->filter_gradient_type() == DensityFilter::FilterGradient::ELEMENTAL){
            ff = this->objective_weights.front()*this->objective.front()->calculate_with_gradient(this, u, new_x, dftmp);
            if(mpi_id == 0){
                this->projection->project_gradient(dftmp, new_x);
                this->filter->filter_gradient(dftmp, df);
                if(this->objective_weights.front() != 1.0){
                    for(auto& v:df){
                        v *= this->objective_weights.front();
                    }
                }
            }
        } else {
            ff = this->objective_weights.front()*this->objective.front()->calculate_with_gradient_nodal(this, u, new_x, dftmp_nodal);
            if(mpi_id == 0){
                this->filter->filter_gradient_nodal(dftmp_nodal, df);
                if(this->objective_weights.front() != 1.0){
                    for(auto& v:df){
                        v *= this->objective_weights.front();
                    }
                }
            }
        }
    } else {
        ff = 0;
        std::fill(df.begin(), df.end(), 0);
        for(size_t i = 0; i < this->objective.size(); ++i){
            if(this->objective[i]->filter_gradient_type() == DensityFilter::FilterGradient::ELEMENTAL){
                ff += this->objective_weights[i]*this->objective[i]->calculate_with_gradient(this, u, new_x, dftmp);
                if(mpi_id == 0){
                    this->projection->project_gradient(dftmp, new_x);
                    this->filter->filter_gradient(dftmp, df);
                    for(size_t j = 0; j < dftmp.size(); ++j){
                        df[j] += this->objective_weights[i]*dftmp[j];
                    }
                }
            } else {
                ff += this->objective_weights[i]*this->objective[i]->calculate_with_gradient_nodal(this, u, new_x, dftmp_nodal);
                if(mpi_id == 0){
                    this->filter->filter_gradient_nodal(dftmp_nodal, dftmp);
                    for(size_t j = 0; j < dftmp.size(); ++j){
                        df[j] += this->objective_weights[i]*dftmp[j];
                    }
                }
            }
        }
    }
    MPI_Bcast(&ff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(mpi_id == 0){
        this->projection->project_gradient(dftmp, x_fil);
        this->filter->filter_gradient(dftmp, df);
    }
    if(M == 1){
        auto& c = this->constraints[0];
        if(c.fun->filter_gradient_type() == DensityFilter::FilterGradient::ELEMENTAL){
            if(c.types[0] == Constraint::Type::LESS_THAN){
                g[0] = c.fun->calculate_with_gradient(this, u, new_x, dgtmp1)
                       - c.bounds[0];
            } else if(c.types[0] == Constraint::Type::GREATER_THAN){
                g[0] = - c.fun->calculate_with_gradient(this, u, new_x, dgtmp1)
                       + c.bounds[0];
            }
            if(mpi_id == 0){
                this->projection->project_gradient(dgtmp1, x_fil);
                this->filter->filter_gradient(dgtmp1, dg);
                if(c.types[0] == Constraint::Type::GREATER_THAN){
                    for(auto& d:dg){
                        d *= -1.0;
                    }
                }
            }
        } else {
            if(c.types[0] == Constraint::Type::LESS_THAN){
                g[0] = c.fun->calculate_with_gradient_nodal(this, u, new_x, dgtmp1_nodal)
                       - c.bounds[0];
            } else if(c.types[0] == Constraint::Type::GREATER_THAN){
                g[0] = - c.fun->calculate_with_gradient_nodal(this, u, new_x, dgtmp1_nodal)
                       + c.bounds[0];
            }
            if(mpi_id == 0){
                //this->projection->project_gradient(dgtmp1_nodal, new_x);
                this->filter->filter_gradient_nodal(dgtmp1_nodal, dg);
                if(c.types[0] == Constraint::Type::GREATER_THAN){
                    for(auto& d:dg){
                        d *= -1.0;
                    }
                }
            }
        }
    } else if(M > 1){
        size_t g_id = 0;
        for(size_t i = 0; i < this->constraints.size(); ++i){
            auto& c = this->constraints[i];
            if(c.fun->filter_gradient_type() == DensityFilter::FilterGradient::ELEMENTAL){
                double val = c.fun->calculate_with_gradient(this, u, new_x, dgtmp1);
                for(size_t k = 0; k < c.types.size(); ++k){
                    if(c.types[k] == Constraint::Type::LESS_THAN){
                        g[g_id] = val - c.bounds[k];
                        if(mpi_id == 0){
                            this->projection->project_gradient(dgtmp1, x_fil);
                            this->filter->filter_gradient(dgtmp1, dgtmp2);
                            for(size_t j = 0; j < x_size; ++j){
                                dg[M*j + g_id] = dgtmp2[j];
                            }
                        }
                        ++g_id;
                    } else if(c.types[k] == Constraint::Type::GREATER_THAN){
                        g[g_id] = - val + c.bounds[k];
                        if(mpi_id == 0){
                            this->projection->project_gradient(dgtmp1, x_fil);
                            this->filter->filter_gradient(dgtmp1, dgtmp2);
                            for(size_t j = 0; j < x_size; ++j){
                                dg[M*j + g_id] = -dgtmp2[j];
                            }
                        }
                        ++g_id;
                    } else if(c.types[k] == Constraint::Type::EQUAL){
                        g[g_id] = - val + c.bounds[k];
                        g[g_id + 1] = -g[g_id];
                        if(mpi_id == 0){
                            this->projection->project_gradient(dgtmp1, x_fil);
                            this->filter->filter_gradient(dgtmp1, dgtmp2);
                            for(size_t j = 0; j < x_size; ++j){
                                dg[M*j + g_id] = dgtmp2[j];
                                dg[M*j + g_id + 1] = -dgtmp2[j];
                            }
                        }
                        g_id += 2;
                    }
                }
            } else {
                double val = c.fun->calculate_with_gradient_nodal(this, u, new_x, dgtmp1_nodal);
                for(size_t k = 0; k < c.types.size(); ++k){
                    if(c.types[k] == Constraint::Type::LESS_THAN){
                        g[g_id] = val - c.bounds[k];
                        if(mpi_id == 0){
                            //this->projection->project_gradient(dgtmp1, x_fil);
                            this->filter->filter_gradient_nodal(dgtmp1_nodal, dgtmp2);
                            for(size_t j = 0; j < x_size; ++j){
                                dg[M*j + g_id] = dgtmp2[j];
                            }
                        }
                        ++g_id;
                    } else if(c.types[k] == Constraint::Type::GREATER_THAN){
                        g[g_id] = - val + c.bounds[k];
                        if(mpi_id == 0){
                            //this->projection->project_gradient(dgtmp1, x_fil);
                            this->filter->filter_gradient_nodal(dgtmp1_nodal, dgtmp2);
                            for(size_t j = 0; j < x_size; ++j){
                                dg[M*j + g_id] = -dgtmp2[j];
                            }
                        }
                        ++g_id;
                    } else if(c.types[k] == Constraint::Type::EQUAL){
                        g[g_id] = val - c.bounds[k];
                        g[g_id + 1] = -g[g_id];
                        if(mpi_id == 0){
                            //this->projection->project_gradient(dgtmp1, x_fil);
                            this->filter->filter_gradient_nodal(dgtmp1_nodal, dgtmp2);
                            for(size_t j = 0; j < x_size; ++j){
                                dg[M*j + g_id] = dgtmp2[j];
                                dg[M*j + g_id + 1] = -dgtmp2[j];
                            }
                        }
                        g_id += 2;
                    }
                }
            }
        }
    }

    // int max_innerit = 30;
    double ch = 1.0;
    int iter;
	for (iter = 0; (ch > this->xtol_abs && std::abs(ff-fnew)/ff > this->ftol_rel); ++iter){

        if(mpi_id == 0){
            for(auto& f:this->objective){
                f->update();
            }
            for(auto& f:this->constraints){
                f.fun->update();
            }
        }

        fnew = ff;

        if(mpi_id == 0){
            mma.Update(x.data(), df.data(), g.data(), dg.data(), xmin.data(), xmax.data());
            this->projection->update(iter);

            ch = 0.0;
            for (size_t i = 0; i < x.size(); ++i) {
                ch = std::max(ch, std::abs(xold[i] - x[i]));
                xold[i] = x[i];
            }

            this->filter->filter_densities(x, x_fil);
            std::copy(x_fil.begin(), x_fil.end(), new_x.begin());
            this->projection->project_densities(new_x);
        }
        MPI_Bcast(&ch, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        u = fem->calculate_displacements(mesh, mesh->load_vector, new_x, this->pc);
        if(mpi_id == 0){
            this->get_stresses(mesh->geometries, u, new_x, this->stresses, this->pc, this->psi);
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
        if(this->objective.size() == 1){
            if(this->objective[0]->filter_gradient_type() == DensityFilter::FilterGradient::ELEMENTAL){
                ff = this->objective_weights.front()*this->objective.front()->calculate_with_gradient(this, u, new_x, dftmp);
                if(mpi_id == 0){
                    this->projection->project_gradient(dftmp, new_x);
                    this->filter->filter_gradient(dftmp, df);
                    if(this->objective_weights.front() != 1.0){
                        for(auto& v:df){
                            v *= this->objective_weights.front();
                        }
                    }
                }
            } else {
                ff = this->objective_weights.front()*this->objective.front()->calculate_with_gradient_nodal(this, u, new_x, dftmp_nodal);
                if(mpi_id == 0){
                    this->filter->filter_gradient_nodal(dftmp_nodal, df);
                    if(this->objective_weights.front() != 1.0){
                        for(auto& v:df){
                            v *= this->objective_weights.front();
                        }
                    }
                }
            }
        } else {
            ff = 0;
            std::fill(df.begin(), df.end(), 0);
            for(size_t i = 0; i < this->objective.size(); ++i){
                if(this->objective[i]->filter_gradient_type() == DensityFilter::FilterGradient::ELEMENTAL){
                    ff += this->objective_weights[i]*this->objective[i]->calculate_with_gradient(this, u, new_x, dftmp);
                    if(mpi_id == 0){
                        this->projection->project_gradient(dftmp, new_x);
                        this->filter->filter_gradient(dftmp, df);
                        for(size_t j = 0; j < dftmp.size(); ++j){
                            df[j] += this->objective_weights[i]*dftmp[j];
                        }
                    }
                } else {
                    ff += this->objective_weights[i]*this->objective[i]->calculate_with_gradient_nodal(this, u, new_x, dftmp_nodal);
                    if(mpi_id == 0){
                        this->filter->filter_gradient_nodal(dftmp_nodal, dftmp);
                        for(size_t j = 0; j < dftmp.size(); ++j){
                            df[j] += this->objective_weights[i]*dftmp[j];
                        }
                    }
                }
            }
        }
        MPI_Bcast(&ff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if(M == 1){
            auto& c = this->constraints[0];
            if(c.fun->filter_gradient_type() == DensityFilter::FilterGradient::ELEMENTAL){
                if(c.types[0] == Constraint::Type::LESS_THAN){
                    g[0] = c.fun->calculate_with_gradient(this, u, new_x, dgtmp1)
                           - c.bounds[0];
                } else if(c.types[0] == Constraint::Type::GREATER_THAN){
                    g[0] = - c.fun->calculate_with_gradient(this, u, new_x, dgtmp1)
                           + c.bounds[0];
                }
                if(mpi_id == 0){
                    this->projection->project_gradient(dgtmp1, x_fil);
                    this->filter->filter_gradient(dgtmp1, dg);
                    if(c.types[0] == Constraint::Type::GREATER_THAN){
                        for(auto& d:dg){
                            d *= -1.0;
                        }
                    }
                }
            } else {
                if(c.types[0] == Constraint::Type::LESS_THAN){
                    g[0] = c.fun->calculate_with_gradient_nodal(this, u, new_x, dgtmp1_nodal)
                           - c.bounds[0];
                } else if(c.types[0] == Constraint::Type::GREATER_THAN){
                    g[0] = - c.fun->calculate_with_gradient_nodal(this, u, new_x, dgtmp1_nodal)
                           + c.bounds[0];
                }
                if(mpi_id == 0){
                    //this->projection->project_gradient(dgtmp1_nodal, new_x);
                    this->filter->filter_gradient_nodal(dgtmp1_nodal, dg);
                    if(c.types[0] == Constraint::Type::GREATER_THAN){
                        for(auto& d:dg){
                            d *= -1.0;
                        }
                    }
                }
            }
        } else if(M > 1){
            size_t g_id = 0;
            for(size_t i = 0; i < this->constraints.size(); ++i){
                auto& c = this->constraints[i];
                if(c.fun->filter_gradient_type() == DensityFilter::FilterGradient::ELEMENTAL){
                    double val = c.fun->calculate_with_gradient(this, u, new_x, dgtmp1);
                    for(size_t k = 0; k < c.types.size(); ++k){
                        if(c.types[k] == Constraint::Type::LESS_THAN){
                            g[g_id] = val - c.bounds[k];
                            if(mpi_id == 0){
                                this->projection->project_gradient(dgtmp1, x_fil);
                                this->filter->filter_gradient(dgtmp1, dgtmp2);
                                for(size_t j = 0; j < x_size; ++j){
                                    dg[M*j + g_id] = dgtmp2[j];
                                }
                            }
                            ++g_id;
                        } else if(c.types[k] == Constraint::Type::GREATER_THAN){
                            g[g_id] = - val + c.bounds[k];
                            if(mpi_id == 0){
                                this->projection->project_gradient(dgtmp1, x_fil);
                                this->filter->filter_gradient(dgtmp1, dgtmp2);
                                for(size_t j = 0; j < x_size; ++j){
                                    dg[M*j + g_id] = -dgtmp2[j];
                                }
                            }
                            ++g_id;
                        } else if(c.types[k] == Constraint::Type::EQUAL){
                            g[g_id] = val - c.bounds[k];
                            g[g_id + 1] = -g[g_id];
                            if(mpi_id == 0){
                                this->projection->project_gradient(dgtmp1, x_fil);
                                this->filter->filter_gradient(dgtmp1, dgtmp2);
                                for(size_t j = 0; j < x_size; ++j){
                                    dg[M*j + g_id] = dgtmp2[j];
                                    dg[M*j + g_id + 1] = -dgtmp2[j];
                                }
                            }
                            g_id += 2;
                        }
                    }
                } else {
                    double val = c.fun->calculate_with_gradient_nodal(this, u, new_x, dgtmp1_nodal);
                    for(size_t k = 0; k < c.types.size(); ++k){
                        if(c.types[k] == Constraint::Type::LESS_THAN){
                            g[g_id] = val - c.bounds[k];
                            if(mpi_id == 0){
                                //this->projection->project_gradient(dgtmp1, x_fil);
                                this->filter->filter_gradient_nodal(dgtmp1_nodal, dgtmp2);
                                for(size_t j = 0; j < x_size; ++j){
                                    dg[M*j + g_id] = dgtmp2[j];
                                }
                            }
                            ++g_id;
                        } else if(c.types[k] == Constraint::Type::GREATER_THAN){
                            g[g_id] = - val + c.bounds[k];
                            if(mpi_id == 0){
                                //this->projection->project_gradient(dgtmp1, x_fil);
                                this->filter->filter_gradient_nodal(dgtmp1_nodal, dgtmp2);
                                for(size_t j = 0; j < x_size; ++j){
                                    dg[M*j + g_id] = -dgtmp2[j];
                                }
                            }
                            ++g_id;
                        } else if(c.types[k] == Constraint::Type::EQUAL){
                            g[g_id] = val - c.bounds[k];
                            g[g_id + 1] = -g[g_id];
                            if(mpi_id == 0){
                                //this->projection->project_gradient(dgtmp1, x_fil);
                                this->filter->filter_gradient_nodal(dgtmp1_nodal, dgtmp2);
                                for(size_t j = 0; j < x_size; ++j){
                                    dg[M*j + g_id] = dgtmp2[j];
                                    dg[M*j + g_id + 1] = -dgtmp2[j];
                                }
                            }
                            g_id += 2;
                        }
                    }
                }
            }
        }

        if(mpi_id == 0){
            logger::quick_log("");
            logger::quick_log("");
            logger::quick_log("Iteration: ", iter);
            logger::quick_log("Results: ", ff);
            logger::quick_log(g);
            logger::quick_log("");
            logger::quick_log("Design var change: ", ch);
            logger::quick_log("fobj change: ", std::abs(ff-fnew)/ff);
            logger::quick_log("");
        }
	}

    if(mpi_id == 0){
        logger::quick_log("");
        auto stop_to = std::chrono::high_resolution_clock::now();
        logger::quick_log("Final result: ", ff);
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
