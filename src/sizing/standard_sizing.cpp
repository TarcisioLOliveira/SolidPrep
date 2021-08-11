/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   Foobar is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Foobar is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with SolidPrep.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "sizing/standard_sizing.hpp"
#include <TopoDS_Wire.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepOffsetAPI_MakePipe.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include "meshing/standard_beam_mesher.hpp"
#include "utils.hpp"
#include "lapacke.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <Precision.hxx>
#include <TopoDS.hxx>
#include <IntTools_EdgeEdge.hxx>
#include "utils/sparse_matrix.hpp"

namespace sizing{

StandardSizing::StandardSizing(ProjectData* data, FiniteElement* solver):
   Sizing(data), solver(solver){}

TopoDS_Shape StandardSizing::run(){
    double mesh_size = 4;

    meshing::StandardBeamMesher mesh(mesh_size, 1, utils::PROBLEM_TYPE_2D);
    TopoDS_Shape beams = this->build_initial_topology();
    auto m = mesh.mesh(beams);
    std::vector<MeshElement*> elems;
    std::vector<double> loads;
    mesh.prepare_for_FEM(m, MeshElementFactory::GT9, this->data);
    auto u = this->solver->calculate_displacements(this->data, &mesh);

    // Define fixed nodes
    std::vector<double> fixed_nodes(this->data->forces.size());
    for(size_t i = 0; i < mesh.node_list.size(); ++i){
        for(size_t j = 0; j < this->data->forces.size(); ++j){
            const gp_Pnt& p1 = mesh.node_list[i]->point;
            const gp_Pnt& p2 = mesh.node_list[fixed_nodes[j]]->point;
            const gp_Pnt& p3 = this->data->forces[j].S.get_centroid();
            if(p1.Distance(p3) < p2.Distance(p3)){
                fixed_nodes[j] = i;
            }
        }
    }

    TopoDS_Shape result = BRepBuilderAPI_Copy(data->ground_structure->shape);
    if(this->data->type == utils::PROBLEM_TYPE_2D){
        std::vector<long> uh(mesh.node_list.size()*2, 0);
        for(auto& n:fixed_nodes){
            uh[n*2] = -1;
            uh[n*2+1] = -1;
        }
        size_t id = 0;
        for(size_t i = 0; i < uh.size(); ++i){
            if(uh[i] > -1){
                uh[i] = id;
                ++id;
            }
        }
        double min_dim = std::numeric_limits<double>::max();
        for(auto& f:this->data->forces){
            if(f.S.get_dimension() < min_dim){
                min_dim = f.S.get_dimension();
            }
        }
        double div = min_dim/mesh_size;

        std::vector<double> h(id, 0);
        double t = this->data->thickness;
        std::vector<double> max_stresses = this->data->material->get_max_stresses(gp_Dir(1, 0, 0));
        double S_x = max_stresses[0];
        double S_y = max_stresses[1];
        double T_xy = max_stresses[2];
        for(size_t i = 0; i < mesh.element_list.size(); ++i){
            auto& e = mesh.element_list[i];
            for(size_t j = 0; j < e->nodes.size(); ++j){
                size_t k = (j+1) % e->nodes.size();
                auto& n1 = e->nodes[j];
                auto& n2 = e->nodes[k];
                gp_Pnt center = n1->point;
                center.BaryCenter(1, n2->point, 1);
                std::vector<double> loads = e->get_loads_at(center, u);
                if(uh[2*i] > -1){
                    double fn = std::abs(loads[1]/(t*S_y));
                    double fs = std::abs(3*loads[0]/(2*t*T_xy));
                    double mf = std::sqrt(6*std::abs(loads[2])/(t*S_y))/div;
                    h[uh[2*i]] = std::max({fn, fs, mf});
                }
                if(uh[2*i+1] > -1){
                    double fn = std::abs(loads[0]/(t*S_x));
                    double fs = std::abs(3*loads[1]/(2*t*T_xy));
                    double mf = std::sqrt(6*std::abs(loads[2])/(t*S_x))/div;
                    h[uh[2*i+1]] = std::max({fn, fs, mf});
                }
            }
        }

        std::vector<double> U = this->calculate_change(&mesh, uh, h, beams);

        for(auto& e:mesh.element_list){
            std::vector<gp_Vec> vecs;
            for(auto& n:e->nodes){
                gp_Vec v(0, 0, 0);
                if(uh[2*n->id] > -1){
                    v.SetX(U[uh[2*n->id]]);
                }
                if(uh[2*n->id+1] > -1){
                    v.SetY(U[uh[2*n->id+1]]);
                }
            }
            result = BRepAlgoAPI_Cut(result, e->get_shape(vecs));
        }
    }
    result = BRepAlgoAPI_Cut(this->data->ground_structure->shape, result);

    return result;
}

TopoDS_Shape StandardSizing::build_initial_topology() const{
    TopoDS_Shape geometry = BRepBuilderAPI_Copy(data->ground_structure->shape);
    for(auto& f:this->data->forces){
        for(auto& s:this->data->supports){
            std::vector<gp_Pnt> b(this->data->pathfinder->find_path(f.S, s.S));
            TopoDS_Wire w;
            for(auto p = b.crbegin()+1; p < b.crend(); ++p){
                TopoDS_Edge e = BRepBuilderAPI_MakeEdge(*(p-1), *p);
                w = BRepBuilderAPI_MakeWire(w, e);
            }
            TopoDS_Shape beam = BRepOffsetAPI_MakePipe(w, f.S.get_shape());
            geometry = BRepAlgoAPI_Cut(geometry, beam);
        }
    }
    geometry = BRepAlgoAPI_Cut(this->data->ground_structure->shape, geometry);

    return geometry;
}

std::vector<double> StandardSizing::calculate_change(BeamMeshing* mesh, const std::vector<long>& ids, std::vector<double> h, const TopoDS_Shape& beams) const{
    size_t dof = 0;
    size_t n = 0;
    if(this->data->type == utils::PROBLEM_TYPE_2D){
        dof = 2;
        n = 3;
    } else if(this->data->type == utils::PROBLEM_TYPE_3D){
        dof = 3;
        n = 4;
    }
    utils::SparseMatrix Km;
    for(size_t i = 0; i < mesh->element_list.size(); ++i){

        std::vector<double> k_mat(n*n*dof*dof, 0);
        std::vector<long> pos;
        for(size_t j = 0; j < n; ++j){
            size_t k = (j+1) % n;
            auto& n1 = mesh->element_list[i]->nodes[j];
            auto& n2 = mesh->element_list[i]->nodes[k];
            for(size_t I = 0; I < dof; ++I){
                k_mat[(I+j)*n+I+j] = 1;
                k_mat[(I+j)*n+I+k+dof] = -1;
            }
            for(size_t l = 0; l < dof; ++l){
                size_t id = n1->id*dof+l;
                if(ids[id] > -1){
                    pos.push_back(n*ids[id]);
                    h[n*ids[id]] = std::max(0.0, h[n*ids[id]] - std::abs(n1->point.Coord(l) - n2->point.Coord(l)));
                }
            }
        }
        Km.insert_matrix(k_mat, pos);
    }

    std::vector<double> h2(h.size(), 0);
    if(this->data->type == utils::PROBLEM_TYPE_2D){
        std::vector<TopoDS_Edge> edges;
        for (TopExp_Explorer exp(this->data->ground_structure->shape, TopAbs_EDGE); exp.More(); exp.Next()){
            edges.push_back(TopoDS::Edge(exp.Current()));
        }
        for (TopExp_Explorer exp(beams, TopAbs_EDGE); exp.More(); exp.Next()){
            edges.push_back(TopoDS::Edge(exp.Current()));
        }
        for(auto& n:mesh->boundary_nodes){
            gp_Lin line(n.node->point, n.normal);
            TopoDS_Edge line_edge = BRepBuilderAPI_MakeEdge(line, 0, Precision::Infinite());
            double min_dist = std::numeric_limits<double>::max();
            for(auto& e:edges){
                IntTools_EdgeEdge tool(line_edge, e);
                tool.Perform();
                auto common = tool.CommonParts();
                if(common.Size() > 0){
                    for(auto& c:common){
                        double dist = c.VertexParameter1(); // n.normal is a unit vector, and the line starts at 0
                        if(dist > 0 && dist < min_dist){
                            min_dist = dist;
                        }
                    }
                }
                if(min_dist < std::numeric_limits<double>::max()){
                    if(ids[n.node->id*2] > -1){
                        h2[ids[n.node->id*2]] = std::abs(min_dist*n.normal.X());
                    }
                    if(ids[n.node->id*2+1] > -1){
                        h2[ids[n.node->id*2+1]] = std::abs(min_dist*n.normal.Y());
                    }
                }
            }
        }
    } else if(this->data->type == utils::PROBLEM_TYPE_3D){
        // TODO
    }
    h2 = Km.multiply(h2);
    for(size_t i = 0; i < h.size(); ++i){
        h[i] = std::max(h[i], h2[i]);
    }
    size_t ku = 0;
    size_t kl = 0;
    std::vector<double> K = Km.to_general_band(h.size(), ku, kl);

    int info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, h.size(), kl, ku, 1, K.data(), h.size(), nullptr, h.data(), 1);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating topology expansions.", info);

    return h;
}

}
