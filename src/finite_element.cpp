/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
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

#include <numeric>
#include <set>
#include <cblas.h>
#include "finite_element.hpp"
#include "logger.hpp"
#include "project_data.hpp"

FiniteElement::FiniteElement(NonlinearSolver* nl)
    :prob_type(nl->get_class()), nl_type(nl->get_type()), nl_solver(nl){

}

void FiniteElement::generate_matrix(const Meshing* const mesh, const size_t u_size, const size_t l_num, const std::vector<long>& node_positions, bool topopt, const std::vector<std::vector<double>>& D_cache){
    MatrixType mtype = MatrixType::RIGID;
    switch(this->prob_type){
        case NonlinearSolver::SolverClass::LINEAR:
            mtype = MatrixType::RIGID;
            break;
        case NonlinearSolver::SolverClass::GRADIENT_BASED:
            mtype = MatrixType::LAMBDA_SLIDING;
            break;
        case NonlinearSolver::SolverClass::HESSIAN_BASED:
            mtype = MatrixType::LAMBDA_HESSIAN;
            break;
    }
    this->generate_matrix_base(mesh, u_size, l_num, node_positions, topopt, D_cache, mtype);
}

void FiniteElement::calculate_displacements(const Meshing* const mesh, std::vector<double>& load, std::vector<double>& lambda, const bool topopt, const std::vector<std::vector<double>>& D_cache){
    switch(this->prob_type){
        case NonlinearSolver::SolverClass::LINEAR:
            this->solve_rigid(load, lambda);
            break;
        case NonlinearSolver::SolverClass::GRADIENT_BASED:
            this->solve_opt(mesh, load, lambda, topopt, D_cache);
            break;
        case NonlinearSolver::SolverClass::HESSIAN_BASED:
            this->solve_newton(mesh, load, lambda);
            break;
    }
}

void FiniteElement::solve_rigid(std::vector<double>& load, std::vector<double>& lambda){
    this->solve(load, lambda);
}
void FiniteElement::solve_opt(const Meshing* const mesh, std::vector<double>& load, std::vector<double>& lambda, const bool topopt, const std::vector<std::vector<double>>& D_cache){
    const size_t l_num = mesh->lambda_elements.size();

    this->nl_solver->setup(l_num);
    std::vector<double> grad(l_num, 0);
    double E = 0;
    size_t it = 0;

    std::vector<double> f(load.size() + 2*l_num);
    std::vector<double> u_ext(mesh->global_load_vector.size(), 0);
    std::vector<double> f_ext(mesh->global_load_vector.size(), 0);
    mesh->extend_vector(0, load, f_ext);

    do{
        // Solve linear subproblem
        std::copy(load.begin(), load.end(), f.begin());
        std::fill(f.begin() + load.size(), f.end(), 0);
        this->apply_lambda_force(mesh, f, lambda, topopt, D_cache);

        this->solve(f, lambda);
        mesh->extend_vector(0, f, u_ext);
        mesh->apply_lambda(lambda, u_ext);

        E = -0.5*cblas_ddot(u_ext.size(), f_ext.data(), 1, u_ext.data(), 1); 

        this->calculate_gradient(mesh, grad, u_ext, f_ext, lambda, topopt, D_cache);

        logger::quick_log("Iteration:", it);
        logger::quick_log("Potential energy:", E);
        logger::quick_log("");

        ++it;

        // Update
    } while(!this->nl_solver->update(lambda.data() + 2*l_num, E, grad.data()));

    std::copy(load.begin(), load.end(), f.begin());
    std::fill(f.begin() + load.size(), f.end(), 0);
    this->apply_lambda_force(mesh, f, lambda, topopt, D_cache);
    this->solve(f, lambda);
    std::copy(f.begin(), f.begin()+load.size(), load.begin());
}
void FiniteElement::solve_newton(const Meshing* const mesh, std::vector<double>& load, std::vector<double>& lambda){
    this->solve(load, lambda);
}

void FiniteElement::apply_lambda_force(const Meshing* const mesh, std::vector<double>& f, const std::vector<double>& lambda, bool topopt, const std::vector<std::vector<double>>& D_cache) const{
    const size_t l_num = mesh->lambda_elements.size();
    const size_t kw = mesh->elem_info->get_k_dimension();
    const size_t lw = mesh->elem_info->get_dof_per_node();
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t t = mesh->thickness;
    const size_t u_size = mesh->load_vector[0].size();

    std::vector<size_t> lv_it(node_num);
    std::vector<const std::vector<size_t>*> lv_vec(node_num);

    const std::vector<long>& npos = mesh->node_positions[0];
    const auto& lambda_list = mesh->lambda_elements;
    std::vector<double> k;
    std::vector<double> R(kw);
    std::vector<double> u_e(kw);
    size_t geom = 0;
    const Geometry* g = mesh->geometries[geom];
    bool use_D_cache = (topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous();
    size_t D_offset = 0;
    std::vector<double> D;
    if(!use_D_cache){
        D = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
    }
    for(const auto& l_e:mesh->lambda_affected_elements){
        if(l_e.geom_id != geom){
            if(!use_D_cache){
                D_offset += g->mesh.size();
            }
            geom = l_e.geom_id;
            const Geometry* g = mesh->geometries[geom];
            use_D_cache = (topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous();
            if(!use_D_cache){
                D = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
            }
        }
        if(use_D_cache){
            k = l_e.e->get_k(D_cache[l_e.e->id - D_offset], t);
        } else {
            k = l_e.e->get_k(D, t);
        }
        const size_t ln = l_e.lambdas.size();
        std::vector<double> R(kw, 0);
        std::vector<double> Rp(kw*ln*2, 0);
        std::vector<double> KR(kw, 0);
        std::vector<double> RpKR(ln*2, 0);

        std::fill(lv_it.begin(), lv_it.end(), 0);
        for(size_t i = 0; i < node_num; ++i){
            const auto& vi = mesh->node_lambda_map.find(l_e.e->nodes[i]->id);
            if(vi != mesh->node_lambda_map.end()){
                lv_vec[i] = &vi->second;
            } else {
                lv_vec[i] = nullptr;
            }
        }
        for(size_t l_it = 0; l_it < ln; ++l_it){
            const size_t l_i = l_e.lambdas[l_it];
            const auto& ll_e = lambda_list[l_i];
            const double ld = lambda[2*l_num + l_i];
            const std::vector<double> R0 =
            {
                ll_e.n.X()*(ld*ld),
                ll_e.n.Y()*(ld*ld),
                ll_e.n.Z()*(ld*ld)
            };
            const std::vector<double> Rp0 =
            {
                ll_e.p1.X(), ll_e.p2.X(),
                ll_e.p1.Y(), ll_e.p2.Y(),
                ll_e.p1.Z(), ll_e.p2.Z()
            };

            for(size_t i = 0; i < node_num; ++i){
                if(lv_vec[i] != nullptr && lv_it[i] < lv_vec[i]->size() && lv_vec[i]->at(lv_it[i]) == l_i){
                    for(size_t ii = 0; ii < lw; ++ii){
                        R[i*lw + ii] += R0[ii];
                        for(size_t jj = 0; jj < 2; ++jj){
                            Rp[(i*lw + ii)*2*ln + (l_it*2 + jj)] = -Rp0[ii*2 + jj];
                        }
                    }
                    ++lv_it[i];
                }
            }
        }

        for(size_t i = 0; i < kw; ++i){
            for(size_t j = 0; j < kw; ++j){
                KR[i] += k[i*kw + j]*R[j];
            }
        }
        for(size_t i = 0; i < 2*ln; ++i){
            for(size_t j = 0; j < kw; ++j){
                RpKR[i] += Rp[j*2*ln + i]*KR[j];
            }
        }

        for(size_t i = 0; i < node_num; ++i){
            const auto& n = l_e.e->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                const size_t p = n->id*dof + j;
                const long fpos = npos[p];
                if(fpos > -1){
                    f[fpos] += KR[(i*dof + j)];
                }
            }
        }
        for(size_t l_it2 = 0; l_it2 < ln; ++l_it2){
            const size_t l_i = l_e.lambdas[l_it2];
            f[u_size + l_i] += RpKR[2*l_it2];
            f[u_size + l_i + l_num] += RpKR[(2*l_it2 + 1)];
        }
    }
}

void FiniteElement::calculate_gradient(const Meshing* const mesh, std::vector<double>& grad, const std::vector<double>& u_ext, const std::vector<double>& f_ext, const std::vector<double>& lambda, bool topopt, const std::vector<std::vector<double>>& D_cache) const{
    const size_t l_num = mesh->lambda_elements.size();
    const size_t kw = mesh->elem_info->get_k_dimension();
    const size_t lw = mesh->elem_info->get_dof_per_node();
    const size_t dof      = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    const size_t t = mesh->thickness;

    std::vector<size_t> lv_it(node_num);
    std::vector<const std::vector<size_t>*> lv_vec(node_num);

    const auto& lambda_list = mesh->lambda_elements;
    std::vector<double> k;
    std::vector<double> R(kw);
    std::vector<double> u_e(kw);
    std::vector<double> f_e(kw);
    size_t geom = 0;
    const Geometry* g = mesh->geometries[geom];
    bool use_D_cache = (topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous();
    size_t D_offset = 0;
    std::vector<double> D;
    if(!use_D_cache){
        D = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
    }
    std::fill(grad.begin(), grad.end(), 0);
    for(const auto& l_e:mesh->lambda_affected_elements){
        if(l_e.geom_id != geom){
            if(!use_D_cache){
                D_offset += g->mesh.size();
            }
            geom = l_e.geom_id;
            const Geometry* g = mesh->geometries[geom];
            use_D_cache = (topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous();
            if(!use_D_cache){
                D = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
            }
        }
        if(use_D_cache){
            k = l_e.e->get_k(D_cache[l_e.e->id - D_offset], t);
        } else {
            k = l_e.e->get_k(D, t);
        }
        const size_t ln = l_e.lambdas.size();
        std::vector<double> R(kw*ln, 0);
        std::vector<double> KR(kw*ln, 0);

        for(size_t i = 0; i < node_num; ++i){
            const auto& n = l_e.e->nodes[i];
            for(size_t j = 0; j < dof; ++j){
                const size_t p = n->id*dof + j;
                u_e[i*dof + j] = u_ext[p];
                f_e[i*dof + j] = f_ext[p];
            }
        }

        std::fill(lv_it.begin(), lv_it.end(), 0);
        for(size_t i = 0; i < node_num; ++i){
            const auto& vi = mesh->node_lambda_map.find(l_e.e->nodes[i]->id);
            if(vi != mesh->node_lambda_map.end()){
                lv_vec[i] = &vi->second;
            } else {
                lv_vec[i] = nullptr;
            }
        }
        for(size_t l_it = 0; l_it < ln; ++l_it){
            const size_t l_i = l_e.lambdas[l_it];
            const auto& ll_e = lambda_list[l_i];
            const double ld = lambda[2*l_num + l_i];
            const std::vector<double> R0 =
            {
                ll_e.n.X()*(2.0*ld),
                ll_e.n.Y()*(2.0*ld),
                ll_e.n.Z()*(2.0*ld)
            };

            for(size_t i = 0; i < node_num; ++i){
                if(lv_vec[i] != nullptr && lv_it[i] < lv_vec[i]->size() && lv_vec[i]->at(lv_it[i]) == l_i){
                    for(size_t ii = 0; ii < lw; ++ii){
                        // Transpose R0
                        R[(i*lw + ii)*ln + l_it] = R0[ii];
                    }
                    ++lv_it[i];
                }
            }
        }

        for(size_t i = 0; i < kw; ++i){
            for(size_t j = 0; j < ln; ++j){
                for(size_t l = 0; l < kw; ++l){
                    KR[i*ln + j] += k[i*kw + l]*R[l*ln + j];
                }
            }
        }

        for(size_t l_it = 0; l_it < ln; ++l_it){
            const size_t l_i = l_e.lambdas[l_it];
            double sum = 0;
            for(size_t i = 0; i < kw; ++i){
                sum += f_e[i]*R[i*ln + l_it] - u_e[i]*KR[i*ln + l_it];
            }
            grad[l_i] += sum;
        }
    }

    /* Parallelizable version (incomplete)
       (hopefully this is unnecessary)

    for(size_t i = 0; i < lambda_list.size(); ++i){
        const auto& l_e = lambda_list[i];
        const double li = lambda[2*l_num + i];
        if(l_e.parent->geom_id != geom){
            if(!use_D_cache){
                D_offset += g->mesh.size();
            }
            geom = l_e.parent->geom_id;
            const Geometry* g = mesh->geometries[geom];
            use_D_cache = (topopt && g->do_topopt) || !g->materials.get_materials()[0]->is_homogeneous();
            if(!use_D_cache){
                D = g->materials.get_D(g->mesh.front().get(), g->mesh.front()->get_centroid());
            }
        }
        if(use_D_cache){
            k = l_e.parent->parent->get_k(D_cache[l_e.parent->parent->id - D_offset], t);
        } else {
            k = l_e.parent->parent->get_k(D, t);
        }

        std::vector<double> R0 =
        {
            l_e.n.X()*(-2.0*li),
            l_e.p1.X()*(-2.0*li),
            l_e.p2.X()*(-2.0*li)
        };
        std::set<MeshElement*> affected_elements;


    }
    */
}

std::vector<double> FiniteElement::calculate_forces(const Meshing* const mesh, const std::vector<double>& displacements) const{
    logger::quick_log("Calculating forces...");
    const size_t dof       = mesh->elem_info->get_dof_per_node();
    const size_t node_num = mesh->elem_info->get_nodes_per_element();
    std::vector<double> results(mesh->node_list.size()*dof, 0);
  
   for(auto& g:mesh->geometries){
       const size_t num_mat = g->number_of_materials();
       if(num_mat == 1){
            for(auto& e:g->mesh){
                const gp_Pnt c = e->get_centroid();
                const auto D = g->materials.get_D(e.get(), c);
                auto f = e->get_internal_loads(D, mesh->thickness, displacements);
                for(size_t n = 0; n < node_num; ++n){
                    auto& node = e->nodes[n];
                    for(size_t i = 0; i < dof; ++i){
                        results[node->id*dof + i] += f[n*dof+i];
                    }
                }
            }
       } else {
           logger::quick_log(num_mat == 1, logger::ERROR, "calculate_forces() not implemented for geometries with more than one material");
       }
    }
    logger::quick_log("Done.");

    return results;
}

