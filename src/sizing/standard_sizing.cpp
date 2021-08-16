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
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepOffsetAPI_MakePipe.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include "logger.hpp"
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
#include <TopTools_ListOfShape.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Wireframe.hxx>
#include <string>

namespace sizing{

StandardSizing::StandardSizing(ProjectData* data, FiniteElement* solver):
   Sizing(data), solver(solver){}

TopoDS_Shape StandardSizing::run(){
    double mesh_size = 1;

    meshing::StandardBeamMesher mesh(mesh_size, 1, utils::PROBLEM_TYPE_2D);
    TopoDS_Shape beams = this->build_initial_topology();
    auto m = mesh.mesh(beams);
    mesh.prepare_for_FEM(m, MeshElementFactory::GT9, this->data);
    auto u = this->solver->calculate_displacements(this->data, &mesh);

    // Define fixed nodes
    std::vector<double> fixed_nodes(this->data->forces.size() + this->end_points.size());
    size_t offset = this->data->forces.size();
    for(size_t i = 0; i < mesh.node_list.size(); ++i){
        for(size_t j = 0; j < this->data->forces.size(); ++j){
            const gp_Pnt& p1 = mesh.node_list[i]->point;
            const gp_Pnt& p2 = mesh.node_list[fixed_nodes[j]]->point;
            const gp_Pnt& p3 = this->data->forces[j].S.get_centroid();
            if(p1.Distance(p3) < p2.Distance(p3)){
                fixed_nodes[j] = i;
            }
        }
        for(size_t j = 0; j < this->end_points.size(); ++j){
            const gp_Pnt& p1 = mesh.node_list[i]->point;
            const gp_Pnt& p2 = mesh.node_list[fixed_nodes[offset+j]]->point;
            const gp_Pnt& p3 = this->end_points[j];
            if(p1.Distance(p3) < p2.Distance(p3)){
                fixed_nodes[offset+j] = i;
            }
        }
    }

    logger::quick_log("Resizing geometry...");
    // TopoDS_Shape result = BRepBuilderAPI_Copy(data->ground_structure->shape);
    TopoDS_Shape result = BRepBuilderAPI_MakeSolid();
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

        std::vector<double> h(id, 0);
        double t = this->data->thickness;
        for(size_t i = 0; i < mesh.element_list.size(); ++i){
            auto& e = mesh.element_list[i];
            for(size_t j = 0; j < e->nodes.size(); ++j){
                size_t k = (j+1) % e->nodes.size();
                auto& n1 = e->nodes[j];
                auto& n2 = e->nodes[k];
                gp_Dir dir(gp_Vec(n2->point, n1->point));
                gp_Pnt center = n1->point;
                center.BaryCenter(1, n2->point, 1);
                gp_Dir perp = dir.Rotated(gp_Ax1(center, gp_Dir(1, 0, 0)), M_PI/2);
                std::vector<double> max_stresses = this->data->material->get_max_stresses(perp);
                double S_x = max_stresses[0];
                double T_xy = max_stresses[2];

                double h1 = n1->point.Distance(n2->point);
                std::vector<double> forces = e->get_average_loads(n2->point, n1->point, u);
                if(forces[0] == 0 && forces[1] == 0){
                    continue;
                }
                gp_Vec F(forces[0], forces[1], 0);
                // double Mz = std::abs(forces[2]);

                double hn = std::abs(perp.Dot(F))/(t*S_x);
                double hs = (F - perp.Dot(F)*perp).Magnitude()/(t*T_xy);
                // double hf = std::sqrt(6*Mz/(t*S_x));

                double hh = std::max({hn, hs});// hf});
                logger::quick_log(hh, h1);
                //hh = std::max(0.0, hh - h1);
                logger::quick_log(forces);
                if(uh[2*n1->id] > -1){
                    // int sign = dir.X()/std::abs(dir.X());
                    // h[uh[2*n1->id]] += sign*std::max(0.0, std::abs(F.X())/(t*S_x) - std::abs(n1->point.X() - n2->point.X()));
                    h[uh[2*n1->id]] += dir.X()*hh;
                }
                if(uh[2*n1->id+1] > -1){
                    // int sign = dir.Y()/std::abs(dir.Y());
                    // h[uh[2*n1->id+1]] += sign*std::max(0.0, std::abs(F.Y())/(t*S_x) - std::abs(n1->point.Y() - n2->point.Y()));
                    h[uh[2*n1->id+1]] += dir.Y()*hh;
                }
            }
        }

        std::vector<double> U = this->calculate_change(&mesh, uh, std::move(h), beams);
        logger::quick_log("Done.");

        logger::quick_log("Generating geometry...");
        std::cout << "0%";
        size_t i = 0;

        result = BRepBuilderAPI_Copy(this->data->ground_structure->shape);
        TopTools_ListOfShape shapes1;
        shapes1.Append(result);
        TopTools_ListOfShape shapes2;
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
                vecs.push_back(v);
            }
            TopoDS_Shape s = e->get_shape(vecs);
            // result = BRepAlgoAPI_Cut(result, s);
            //result = this->simplify_shape(result);
            
            shapes2.Clear();
            shapes2.Append(s);
            BRepAlgoAPI_Cut cut;
            cut.SetTools(shapes2);
            cut.SetArguments(shapes1);
            cut.SetNonDestructive(false);
            cut.SimplifyResult();
            cut.SetRunParallel(true);
            cut.Build();
            result = cut;
            shapes1.Clear();
            shapes1.Append(result);

            double pc = i/(double)(mesh.element_list.size()-1);
            std::cout << "\r" << pc*100 << "%         ";

            ++i;
        }
        std::cout << "\r" << 100 << "%         ";
        std::cout << std::endl;

        logger::quick_log("Done.");
    }
    result = BRepAlgoAPI_Cut(this->data->ground_structure->shape, result);
    result = this->simplify_shape(result);

    return result;
}
TopoDS_Shape StandardSizing::simplify_shape(TopoDS_Shape shape) const{
    double tol = 10;
    double prec = tol*1.5;

    Handle(ShapeFix_Shape) sfs = new ShapeFix_Shape(shape);
    sfs->SetPrecision(prec);
    sfs->SetMaxTolerance(2*tol);
    sfs->SetMinTolerance(tol);
    auto sfw = sfs->FixWireTool();
    sfw->ModifyGeometryMode() = true;
    sfw->ModifyTopologyMode() = true;
    sfw->FixSmallMode() = true;
    sfw->FixSmall(false, prec);
    sfs->Perform();
    shape = sfs->Shape();

    Handle(ShapeFix_Wireframe) SFWF = new ShapeFix_Wireframe(shape);
    SFWF->SetPrecision(prec);
    SFWF->SetMaxTolerance(2*tol);
    SFWF->SetMinTolerance(tol);
    SFWF->ModeDropSmallEdges() = Standard_True;
    SFWF->FixSmallEdges();
    SFWF->FixWireGaps();
    shape = SFWF->Shape();

    return shape;
}

TopoDS_Shape StandardSizing::build_initial_topology(){
    TopoDS_Shape geometry = BRepBuilderAPI_Copy(data->ground_structure->shape);
    this->end_points.reserve(this->data->supports.size());
    for(auto& f:this->data->forces){
        for(auto& s:this->data->supports){
            std::vector<gp_Pnt> b(this->data->pathfinder->find_path(f.S, s.S));
            this->end_points.push_back(b.front());
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

    return this->simplify_shape(geometry);
}

std::vector<double> StandardSizing::calculate_change(BeamMeshing* mesh, const std::vector<long>& ids, std::vector<double> h, const TopoDS_Shape& beams) const{
    double tol = 1e-13;
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
            Node* n1 = mesh->element_list[i]->nodes[j];
            Node* n2 = mesh->element_list[i]->nodes[k];
            for(size_t I = 0; I < dof; ++I){
                k_mat[(I+j*dof)*n*dof+I+j*dof] = 1;
                k_mat[(I+j*dof)*n*dof+I+k*dof] = -1;
            }
            for(size_t l = 0; l < dof; ++l){
                size_t id = n1->id*dof+l;
                pos.push_back(ids[id]);
                // if(ids[id] > -1){
                //     h[ids[id]] = std::max(0.0, h[ids[id]] - std::abs(n1->point.Coord(l+1) - n2->point.Coord(l+1)));
                // }
            }
        }
        Km.insert_matrix(k_mat, pos);
    }

    std::vector<double> h2(h.size(), 0);
    std::vector<size_t> boundary_ids;
    if(this->data->type == utils::PROBLEM_TYPE_2D){
        std::vector<TopoDS_Edge> edges_init;
        std::vector<TopoDS_Edge> edges_beam;
        for (TopExp_Explorer exp(this->data->ground_structure->shape, TopAbs_EDGE); exp.More(); exp.Next()){
            edges_init.push_back(TopoDS::Edge(exp.Current()));
        }
        for (TopExp_Explorer exp(beams, TopAbs_EDGE); exp.More(); exp.Next()){
            edges_beam.push_back(TopoDS::Edge(exp.Current()));
        }
        for(auto& n:mesh->boundary_nodes){
            gp_Lin line(n.node->point, n.normal);
            TopoDS_Edge line_edge = BRepBuilderAPI_MakeEdge(line, 0, Precision::Infinite());
            double min_dist = std::numeric_limits<double>::max();
            for(auto& e:edges_init){
                if(min_dist < tol){
                    min_dist = 0;
                    break;
                }
                IntTools_EdgeEdge tool(line_edge, e);
                tool.Perform();
                auto common = tool.CommonParts();
                if(common.Size() > 0){
                    for(auto& c:common){
                        double dist = std::abs(c.VertexParameter1()); // n.normal is a unit vector, and the line starts at 0
                        if(dist > 0 && dist < min_dist){
                            min_dist = dist;
                        }
                    }
                }
            }
            for(auto& e:edges_beam){
                if(min_dist < tol){
                    min_dist = 0;
                    break;
                }
                IntTools_EdgeEdge tool(line_edge, e);
                tool.Perform();
                auto common = tool.CommonParts();
                if(common.Size() > 0){
                    for(auto& c:common){
                        double dist = std::abs(c.VertexParameter1()); // n.normal is a unit vector, and the line starts at 0
                        if(dist > 0 && dist < min_dist){
                            min_dist = dist;
                        }
                    }
                }
            }
            if(min_dist < std::numeric_limits<double>::max()){
                if(ids[n.node->id*2] > -1){
                    boundary_ids.push_back(ids[n.node->id*2]);
                    h2[ids[n.node->id*2]] = min_dist*n.normal.X();
                }
                if(ids[n.node->id*2+1] > -1){
                    boundary_ids.push_back(ids[n.node->id*2+1]);
                    h2[ids[n.node->id*2+1]] = min_dist*n.normal.Y();
                }
            }
        }
    } else if(this->data->type == utils::PROBLEM_TYPE_3D){
        // TODO
    }
    // h2 = Km.multiply(h2);
    // std::vector<size_t> affected = Km.affected_ids(boundary_ids);
    // for(size_t i = 0; i < affected.size(); ++i){
    //     if(std::abs(h2[affected[i]]) < std::abs(h[affected[i]])){
    //         h[affected[i]] = h2[affected[i]];
    //     }
    // }
    size_t ku = 0;
    size_t kl = 0;
    std::vector<double> K = Km.to_general_band(h.size(), ku, kl);
    std::vector<int> pivot(h.size(), 0);

    //logger::quick_log(h);
    int info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, h.size(), kl, ku, 1, K.data(), h.size(), pivot.data(), h.data(), 1);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating topology expansions.", info);
    logger::quick_log(h);

    return h;
}

}
