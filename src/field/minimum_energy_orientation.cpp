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

#include <algorithm>
#include <set>
#include "logger.hpp"
#include "math/matrix.hpp"
#include "math/slicer.hpp"
#include "math/vector_view.hpp"
#include "project_data.hpp"
#include "field/minimum_energy_orientation.hpp"
#include "project_specification/registry.hpp"
#include "utils/basis_tensor.hpp"

namespace field{

MinimumEnergyOrientation::MinimumEnergyOrientation(const projspec::DataMap& data):
    proj_data(data.proj),
    elem_info(data.proj->topopt_element),
    geoms(),
    thickness(data.proj->thickness),
    show(data.get_bool("display", true)),
    DIM((elem_info->get_problem_type() == utils::PROBLEM_TYPE_2D) ? 2 : 3),
    NODES_PER_ELEM(elem_info->get_nodes_per_element()),
    max_it(data.get_int("max_it")),
    use_strain(data.get_bool("use_strain")),
    R0({0, 0, 1,
        0, 1, 0,
        -1, 0, 0}, 3, 3)
    {
   
    if(!data.get_bool("STUB")){ 
    }
}

void MinimumEnergyOrientation::generate(){
    this->geoms.clear();

    const size_t NODE_DOF = this->elem_info->get_dof_per_node();
    const size_t BNODE_NUM = this->elem_info->get_boundary_nodes_per_element();

    size_t elem_num = 0;
    {
        auto& proj_geoms = this->proj_data->geometries;
        std::vector<Geometry*> picked_geoms;
        picked_geoms.reserve(proj_geoms.size());
        for(auto& g:proj_geoms){
            auto mats = g->materials.get_materials();
            for(auto& m:mats){
                if(m->get_field() == this){
                    picked_geoms.push_back(g.get());
                    elem_num += g->mesh.size();
                    break;
                }
            }
        }
        this->geoms.insert(this->geoms.begin(), picked_geoms.begin(), picked_geoms.end());
    }
    this->theta.resize(elem_num, 0);
    this->omega.resize(elem_num, 0);
    this->phi.resize(elem_num, 0);

    this->matrices.resize(elem_num*DIM*DIM); // row major;

    if(this->geoms.size() == 0){
        return;
    }

    Meshing* mesh = this->proj_data->topopt_mesher.get();
    SolverManager* fem = this->proj_data->topopt_fea.get();

    std::vector<std::vector<utils::DelayedPointerView<Material>>> materials_backup(this->geoms.size());

    std::vector<std::vector<utils::DelayedPointerView<Material>>> ortho_material(this->geoms.size());
    {
        size_t elem_pos = 0;
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
            g->set_materials(ortho_material[i]);
            for(auto& e:g->mesh){
                this->elem_pos_map[e->id] = elem_pos;
                ++elem_pos;
            }
        }
    }

    std::vector<size_t> int_loads_id;
    std::vector<size_t> spring_id;
    if(this->proj_data->internal_loads.size() > 0){
        const size_t int_size = this->proj_data->internal_loads.size();
        int_loads_id.reserve(int_size);
        for(size_t i = 0; i < int_size; ++i){
            auto& inti = this->proj_data->internal_loads[i];
            if(inti.mat->get_field() == this){
                int_loads_id.push_back(i);
                inti.set_calculate_adjoint();
            }
        }
    }
    if(this->proj_data->springs.size() > 0){
        const size_t int_size = this->proj_data->springs.size();
        spring_id.reserve(int_size);
        for(size_t i = 0; i < int_size; ++i){
            auto& inti = this->proj_data->springs[i];
            if(inti.mat->get_field() == this){
                spring_id.push_back(i);
            }
        }
    }
    bool update_boundary_conditions = false;
    if(int_loads_id.size() > 0 || spring_id.size() > 0){
        update_boundary_conditions = true;
        mesh->apply_boundary_conditions(this->proj_data->forces, this->proj_data->supports, this->proj_data->springs, this->proj_data->internal_loads, this->proj_data->sub_problems);
    }

    std::vector<double> u(mesh->max_dofs, 0);
    auto u_old = u;
    double c_old = 0;

    const size_t var_size = elem_num;
    std::vector<double> dfdt(var_size, 0);
    std::vector<double> dfdo(var_size, 0);
    std::vector<double> dfdp(var_size, 0);

    const double s0 = 0.6;
    std::vector<double> st(var_size, s0);
    std::vector<double> so(var_size, s0);
    std::vector<double> sp(var_size, s0);

    std::vector<double> told1(var_size), told2(var_size);
    std::vector<double> pold1(var_size), pold2(var_size);
    std::vector<double> oold1(var_size), oold2(var_size);

    const double asymp_inc = 1.0;
    const double asymp_dec = 0.7;

    const double tmin = -std::numbers::pi_v<double>;
    const double tmax =  std::numbers::pi_v<double>;
    const double omin = -std::numbers::pi_v<double>;
    const double omax =  std::numbers::pi_v<double>;
    const double pmin = -std::numbers::pi_v<double>;
    const double pmax =  std::numbers::pi_v<double>;

    const auto sum_transpose = [](const math::Matrix& M){
        return M + M.T();
    };

    double du = 1.0;
    size_t it = 0;

    const auto run_mma = 
        [&](std::vector<double>& x, std::vector<double>& s, const std::vector<double>& dfdx,
           const double xmin, const double xmax, std::vector<double>& xold1, std::vector<double>& xold2){

        const double dxm = xmax - xmin;
        for(size_t i = 0; i < var_size; ++i){
            const double u = std::min(x[i] + s[i], xmax);
            const double l = std::max(x[i] - s[i], xmin);

            const double ux = u - x[i];
            const double lx = x[i] - l;
            const double dfp = std::max(0.0,  dfdx[i]);
            const double dfm = std::max(0.0, -dfdx[i]);
            const double p = ux*ux*( 1.001*dfp + 0.001*dfm + 1e-5/dxm );
            const double q = lx*lx*( 0.001*dfp + 1.001*dfm + 1e-5/dxm );

            const double sqp = std::sqrt(p);
            const double sqq = std::sqrt(q);

            x[i] = (u*sqq + l*sqp)/(sqq + sqp);
            du = std::max(du, std::abs(x[i] - xold1[i]));
            if(it > 1){
                if((xold2[i] - xold1[i])*(xold1[i] - x[i]) > 0){
                    s[i] *= asymp_inc;
                } else {
                    s[i] *= asymp_dec;
                }
                xold2[i] = xold1[i];
                xold1[i] = x[i];
            } else {
                xold2[i] = xold1[i];
                xold1[i] = x[i];
            }
        }
    };

    logger::quick_log("Generating minimum energy orientation field....");
    logger::quick_log("");
    while(du > 1e-1){
        fem->update_materials();
        fem->generate_matrix(mesh);

        fem->calculate_displacements_global(mesh, mesh->load_vector, u);


        this->lock = true;
        for(auto& g:this->geoms){
            for(auto& e:g->mesh){
                const size_t id = this->elem_pos_map[e->id];
                const auto pos = math::slicer::from_node_upos(e->nodes, NODES_PER_ELEM, NODE_DOF);
                const math::VectorSliceGeneralView us(u, pos);

                const double theta = this->theta[id];
                const double omega = this->omega[id];
                const double phi = this->phi[id];
                const auto Rt = Rtheta(theta);
                const auto Ro = Romega(omega);
                const auto Rp = Rphi(phi);
                const auto R = Rp*Ro*Rt*R0;
                const auto dRt = Rtheta_d1(theta);
                const auto dRo = Romega_d1(omega);
                const auto dRp = Rphi_d1(phi);
                const auto dRdt = Rp*Ro*dRt*R0;
                const auto dRdo = Rp*dRo*Rt*R0;
                const auto dRdp = dRp*Ro*Rt*R0;

                const auto D0 = g->materials.get_D(e.get(), e->get_centroid());
                
                const auto T(utils::basis_tensor_3D(R.T()));

                const auto dTdt(utils::basis_tensor_3D_d1(R.T(), dRdt.T()));
                const auto dTdo(utils::basis_tensor_3D_d1(R.T(), dRdo.T()));
                const auto dTdp(utils::basis_tensor_3D_d1(R.T(), dRdp.T()));

                const auto dTDTdt = sum_transpose(T*D0*dTdt.T());
                const auto dTDTdo = sum_transpose(T*D0*dTdo.T());
                const auto dTDTdp = sum_transpose(T*D0*dTdp.T());

                const auto dKdt = e->get_k(dTDTdt, this->thickness);
                const auto dKdo = e->get_k(dTDTdo, this->thickness);
                const auto dKdp = e->get_k(dTDTdp, this->thickness);

                dfdt[id] = -0.5*(dKdt*us).dot(us);
                dfdo[id] = -0.5*(dKdo*us).dot(us);
                dfdp[id] = -0.5*(dKdp*us).dot(us);
            }
        }
        for(auto& int_id:int_loads_id){
            auto& il = this->proj_data->internal_loads[int_id];
            il.curvature_adjoint(u);
            for(auto& be:il.boundary_mesh){
                const auto e = be->parent;
                const auto id = this->elem_pos_map[e->id];

                const double theta = this->theta[id];
                const double omega = this->omega[id];
                const double phi = this->phi[id];
                const auto Rt = Rtheta(theta);
                const auto Ro = Romega(omega);
                const auto Rp = Rphi(phi);
                const auto R = Rp*Ro*Rt*R0;
                const auto dRt = Rtheta_d1(theta);
                const auto dRo = Romega_d1(omega);
                const auto dRp = Rphi_d1(phi);
                const auto dRdt = Rp*Ro*dRt*R0;
                const auto dRdo = Rp*dRo*Rt*R0;
                const auto dRdp = dRp*Ro*Rt*R0;

                const auto S0 = il.mat->stiffness_inverse_3D(e, e->get_centroid());
                const auto D0 = S0.get_inverted_cholesky();
                
                const auto T(utils::basis_tensor_3D(R.T()));
                const auto Tinv(utils::basis_tensor_3D_inv_T(R.T()));

                const auto dTdt(utils::basis_tensor_3D_d1(R.T(), dRdt.T()));
                const auto dTinvdt(utils::basis_tensor_3D_inv_T_d1(R.T(), dRdt.T()));
                const auto dTdo(utils::basis_tensor_3D_d1(R.T(), dRdo.T()));
                const auto dTinvdo(utils::basis_tensor_3D_inv_T_d1(R.T(), dRdo.T()));
                const auto dTdp(utils::basis_tensor_3D_d1(R.T(), dRdp.T()));
                const auto dTinvdp(utils::basis_tensor_3D_inv_T_d1(R.T(), dRdp.T()));

                const auto TDT = T*D0*T.T();
                const auto TST = Tinv*S0*Tinv.T();

                const auto dTDTdt = sum_transpose(T*D0*dTdt.T());
                const auto dTSTdt = sum_transpose(Tinv*S0*dTinvdt.T());
                const auto dTDTdo = sum_transpose(T*D0*dTdo.T());
                const auto dTSTdo = sum_transpose(Tinv*S0*dTinvdo.T());
                const auto dTDTdp = sum_transpose(T*D0*dTdp.T());
                const auto dTSTdp = sum_transpose(Tinv*S0*dTinvdp.T());

                dfdt[id] += il.get_adjoint_gradient(u, TST, TDT, dTSTdt, dTDTdt, be.get());
                dfdo[id] += il.get_adjoint_gradient(u, TST, TDT, dTSTdo, dTDTdo, be.get());
                dfdp[id] += il.get_adjoint_gradient(u, TST, TDT, dTSTdp, dTDTdp, be.get());
            }
        }
        for(auto& spr_id:spring_id){
            const auto& spr = this->proj_data->springs[spr_id];
            std::vector<gp_Pnt> points(BNODE_NUM);
            const auto Rspr = utils::basis_tensor_3D_inv_T(spr.rot3D);
            for(auto& be:spr.submesh){
                const auto e = be->parent;
                const auto id = this->elem_pos_map[e->id];
                const auto pos = math::slicer::from_node_upos(e->nodes, NODES_PER_ELEM, NODE_DOF);
                const math::VectorSliceGeneralView us(u, pos);

                for(size_t i = 0; i < BNODE_NUM; ++i){
                    points[i] = be->nodes[i]->point;
                }

                const double theta = this->theta[id];
                const double omega = this->omega[id];
                const double phi = this->phi[id];
                const auto Rt = Rtheta(theta);
                const auto Ro = Romega(omega);
                const auto Rp = Rphi(phi);
                const auto R = Rp*Ro*Rt*R0;
                const auto dRt = Rtheta_d1(theta);
                const auto dRo = Romega_d1(omega);
                const auto dRp = Rphi_d1(phi);
                const auto dRdt = Rp*Ro*dRt*R0;
                const auto dRdo = Rp*dRo*Rt*R0;
                const auto dRdp = dRp*Ro*Rt*R0;

                const auto S0 = spr.mat->stiffness_inverse_3D(e, e->get_centroid());
                
                const auto Tinv(utils::basis_tensor_3D_inv_T(R.T()));

                const auto dTinvdt(utils::basis_tensor_3D_inv_T_d1(R.T(), dRdt.T()));
                const auto dTinvdo(utils::basis_tensor_3D_inv_T_d1(R.T(), dRdo.T()));
                const auto dTinvdp(utils::basis_tensor_3D_inv_T_d1(R.T(), dRdp.T()));

                const auto TST = Rspr*Tinv*S0*Tinv.T()*Rspr.T();

                const auto dTSTdt = sum_transpose(Rspr*Tinv*S0*dTinvdt.T()*Rspr.T());
                const auto dTSTdo = sum_transpose(Rspr*Tinv*S0*dTinvdo.T()*Rspr.T());
                const auto dTSTdp = sum_transpose(Rspr*Tinv*S0*dTinvdp.T()*Rspr.T());

                const auto dKdt = spr.get_dK(TST, dTSTdt);
                const auto dKdo = spr.get_dK(TST, dTSTdo);
                const auto dKdp = spr.get_dK(TST, dTSTdp);

                const auto dRedt = e->get_R(dKdt, this->thickness, points);
                const auto dRedo = e->get_R(dKdo, this->thickness, points);
                const auto dRedp = e->get_R(dKdp, this->thickness, points);

                dfdt[id] += -0.5*(dRedt*us).dot(us);
                dfdo[id] += -0.5*(dRedo*us).dot(us);
                dfdp[id] += -0.5*(dRedp*us).dot(us);
            }
        }
        du = 0;

        run_mma(theta, st, dfdt, tmin, tmax, told1, told2);
        run_mma(omega, so, dfdo, omin, omax, oold1, oold2);
        run_mma(phi, sp, dfdp, pmin, pmax, pold1, pold2);

        this->lock = false;
        //for(size_t i = 0; i < mesh->max_dofs; ++i){
        //    du = std::max(du, std::abs(u[i] - u_old[i]));
        //    u_old[i] = u[i];
        //}

        double c = 0;
        for(size_t i = 0; i < u.size(); ++i){
            c += u[i]*mesh->global_load_vector[i];
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
    for(auto& int_id:int_loads_id){
        auto& il = this->proj_data->internal_loads[int_id];
        il.clear_curvature_data();
    }

    for(size_t i = 0; i < this->geoms.size(); ++i){
        auto& g = this->geoms[i];
        g->set_materials(materials_backup[i]);
    }

    fem->update_materials();
    fem->generate_matrix(mesh);

    logger::quick_log("");
    logger::quick_log("Finished.");

}


void MinimumEnergyOrientation::initialize_views(Visualization* viz){
    if(this->show && this->geoms.size() != 0){
        this->L_view = viz->add_view("Longitudinal Vector", spview::defs::ViewType::VECTOR, spview::defs::DataType::OTHER);
        this->R_view = viz->add_view("Radial Vector", spview::defs::ViewType::VECTOR, spview::defs::DataType::OTHER);
        this->T_view = viz->add_view("Tangential Vector", spview::defs::ViewType::VECTOR, spview::defs::DataType::OTHER);
    }
}

void MinimumEnergyOrientation::display_views() const{
    if(this->geoms.size() == 0){
        return;
    }
    if(this->show){
        const size_t nodes_per_elem = elem_info->get_nodes_per_element();

        std::map<size_t, size_t> id_pos_map;
        std::set<Node*> unique_nodes;
        std::vector<size_t> ids(this->geoms.size(), 0);
        for(size_t i = 0; i < ids.size(); ++i){
            ids[i] = this->geoms[i]->id;
            unique_nodes.insert(this->geoms[i]->node_list.begin(), this->geoms[i]->node_list.end());
        }
        size_t new_id = 0;
        for(const auto& n:unique_nodes){
            id_pos_map[n->id] = new_id;
            ++new_id;
        }
        const size_t node_num = unique_nodes.size();
        unique_nodes.clear();
        std::vector<double> L_vecs(node_num*3, 0);
        std::vector<double> R_vecs(node_num*3, 0);
        std::vector<double> T_vecs(node_num*3, 0);
        // Assuming per-element rotation matrix storage
        for(const auto& g:geoms){
            for(const auto& e:g->mesh){
                const auto dirs = this->get_matrix(e.get(), e->get_centroid());
                for(size_t i = 0; i < nodes_per_elem; ++i){
                    const auto n = e->nodes[i];
                    const size_t id = id_pos_map.at(n->id);
                    for(size_t j = 0; j < 3; ++j){
                        T_vecs[id*3 + j] += dirs(j,0);
                        R_vecs[id*3 + j] += dirs(j,1);
                        L_vecs[id*3 + j] += dirs(j,2);
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
        this->L_view->update_view(L_vecs, ids);
        this->R_view->update_view(R_vecs, ids);
        this->T_view->update_view(T_vecs, ids);
    }
}
math::Matrix MinimumEnergyOrientation::get_matrix(const MeshElement* e, const gp_Pnt& p) const{
    (void)p;
    if(this->first_run || this->lock){
        return math::Matrix::identity(this->DIM);
    } else {
        const auto eid = this->elem_pos_map.at(e->id);
        const double theta = this->theta[eid];
        const double omega = this->omega[eid];
        const double phi = this->phi[eid];
        const auto Rt = Rtheta(theta);
        const auto Ro = Romega(omega);
        const auto Rp = Rphi(phi);
        const auto R = Rp*Ro*Rt*R0;

        return R;
    }
}

std::array<gp_Dir, 3> MinimumEnergyOrientation::get_array(const MeshElement* e, const gp_Pnt& p) const{
    auto M = this->get_matrix(e, p);
    if(this->DIM == 3){
        return {
            gp_Dir(M(0,0), M(1,0), M(2,0)),
            gp_Dir(M(0,1), M(1,1), M(2,1)),
            gp_Dir(M(0,2), M(1,2), M(2,2)),
        };
    } else {
        return {
            gp_Dir(M(0,0), M(1,0), 0),
            gp_Dir(M(0,1), M(1,1), 0),
            gp_Dir(     0,      0, 1),
        };
    }
}

using namespace projspec;
const bool MinimumEnergyOrientation::reg = Factory<Field>::add(
    [](const DataMap& data){
        return std::make_unique<MinimumEnergyOrientation>(data);
    },
    ObjectRequirements{
        "minimum_energy_orientation",
        {
            DataEntry{.name = "display", .type = TYPE_BOOL, .required = false},
        }
    }
);


}
