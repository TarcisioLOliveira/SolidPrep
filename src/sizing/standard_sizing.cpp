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
#include <cstdlib>
#include <gp_Circ.hxx>
#include <NCollection_DataMap.hxx>
#include <NCollection_BaseMap.hxx>
#include <NCollection_StlIterator.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepOffsetAPI_MakePipe.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include "logger.hpp"
#include "meshing/standard_beam_mesher.hpp"
#include "utils.hpp"
#include "lapacke.h"
#include "cblas.h"
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
#include <vector>
#include <BOPAlgo_Splitter.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <GeomAPI_ExtremaCurveCurve.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <project_data.hpp>

namespace sizing{

StandardSizing::StandardSizing(ProjectData* data, FiniteElement* solver, double element_size, double multiplier):
   Sizing(data), solver(solver), element_size(element_size), multiplier(multiplier){}

TopoDS_Shape StandardSizing::run(){
    return this->boundary_expansion_approach();
    // return this->experimental_elemental_approach();
}

TopoDS_Shape StandardSizing::boundary_expansion_approach(){

    // Get beams and do FEA
    meshing::StandardBeamMesher mesh(this->element_size, 1, utils::PROBLEM_TYPE_2D);
    TopoDS_Shape beams = this->build_initial_topology();
    utils::shape_to_file("beams.step", beams);
    auto m = mesh.mesh(beams);
    mesh.prepare_for_FEM(m, MeshElementFactory::GT9, this->data);
    auto u = this->solver->calculate_displacements(this->data, &mesh);
    this->solver->calculate_forces(&mesh, u);

    if(this->data->type == utils::PROBLEM_TYPE_2D){
        return this->expansion_2D(mesh, u, beams);
    } else if(this->data->type == utils::PROBLEM_TYPE_3D){
        // TODO
    }

    return TopoDS_Shape();
}

TopoDS_Shape StandardSizing::expansion_2D(const meshing::StandardBeamMesher& mesh, const std::vector<double>& u, const TopoDS_Shape& beams){
    logger::quick_log("Calculating new geometry...");
    std::cout << "0%";

    // Get edges
    std::vector<TopoDS_Edge> edges_init;
    std::vector<TopoDS_Edge> edges_beam;
    for (TopExp_Explorer exp(this->data->ground_structure->shape, TopAbs_EDGE); exp.More(); exp.Next()){
        edges_init.push_back(TopoDS::Edge(exp.Current()));
    }
    for (TopExp_Explorer exp(beams, TopAbs_EDGE); exp.More(); exp.Next()){
        edges_beam.push_back(TopoDS::Edge(exp.Current()));
    }

    std::vector<ExternalForce> external_forces;
    // Get applied forces
    for(auto& f:this->data->forces){
        ExternalForce ef;
        ef.forces = f.vec;
        ef.position = f.S.get_centroid();
        ef.diameter = f.S.get_dimension();
        gp_Dir normal(f.S.get_normal());
        ef.line_dir = gp_Dir(normal.Rotated(gp_Ax1(ef.position, gp_Dir(0,0,1)), M_PI/2));
        external_forces.push_back(ef);
    }
    // Get reaction forces
    size_t Mn = 0;
    for(size_t i = 0; i < this->end_points.size(); ++i){
        // Get supports
        gp_Pnt cur_p = this->end_points[i];
        Support* sup = nullptr;
        for(auto& s:this->data->supports){
            if(s.S.is_inside(cur_p)){
                sup = &s;
                break;
            }
        }
        if(sup == nullptr){
            continue;
        }

        // Get edges
        gp_Dir normal(sup->S.get_normal());
        gp_Dir line_dir(normal.Y(), -normal.X(), 0);
        gp_Lin line(cur_p, line_dir);
        TopoDS_Edge line_edge = BRepBuilderAPI_MakeEdge(line, -Precision::Infinite(), Precision::Infinite());
        Handle(Geom_Line) geom_line(new Geom_Line(line));
        gp_Pnt p1, p2;
        double dist1 = Precision::Infinite();
        double dist2 = Precision::Infinite();
        if(dist1 - Precision::Infinite() >= -Precision::Confusion() || dist2 - Precision::Infinite() >= -Precision::Confusion()){
            for(TopExp_Explorer exp(beams, TopAbs_VERTEX); exp.More(); exp.Next()){
                gp_Pnt p = BRep_Tool::Pnt(TopoDS::Vertex(exp.Current()));
                if(line.Contains(p, Precision::Confusion())){
                    double dist = cur_p.Distance(p);
                    if(dist > Precision::Confusion()){
                        if(cur_p.Translated(dist*line_dir).IsEqual(p, Precision::Confusion()) && dist < dist1){
                            p1 = cur_p.Translated(dist*(line_dir));
                            dist1 = dist;
                        } else if(cur_p.Translated(-dist*line_dir).IsEqual(p, Precision::Confusion()) && dist < dist2){
                            p2 = cur_p.Translated(-dist*(line_dir));
                            dist2 = dist;
                        }
                    }
                }
            }
        }

        // Calculate reactions along cross-section
        gp_Pnt center(p1);
        center.BaryCenter(1, p2, 1);
        double distance = p1.Distance(p2);
        gp_Dir final_dir(gp_Vec(p1, p2));
        gp_Lin final_line(p1, final_dir);

        // Filter redundancy
        for(size_t j = i+1; j < this->end_points.size(); ++j){
            if(final_line.Contains(this->end_points[j], Precision::Confusion()) && center.Distance(this->end_points[j]) - distance/2 <= Precision::Confusion()){
                this->end_points.erase(this->end_points.begin()+j);
                --j;
            }
        }

        double Fx = 0;
        double Fy = 0;
        for(auto& n:mesh.node_list){
            if(final_line.Contains(n->point, Precision::Confusion())){
                double dist = p1.Distance(n->point);
                if(p1.Translated(dist*final_dir).IsEqual(n->point, Precision::Confusion())){
                    Fx -= n->results[0];
                    Fy -= n->results[1];
                }
            }
        }
        ExternalForce ef;
        ef.forces = gp_Vec(Fx, Fy, 0);
        ef.position = center;
        ef.Mx = sup->MX;
        ef.My = sup->MY;
        ef.Mz = sup->MZ;
        if(ef.Mz){
            ++Mn;
        }
        ef.diameter = distance;
        ef.line_dir = line_dir;
        external_forces.push_back(ef);
    }

    // Calculate reaction moments
    if(Mn > 1){
        // Fast inverted matrix because LAPACK was returning inf
        // Basically the inverse of ones() - I
        std::vector<double> Mat(Mn*Mn, 1.0/std::max(1.0, Mn-1.0));
        for(size_t i = 0; i < Mn; ++i){
            Mat[i*Mn + i] *= -std::max(0.0, Mn-2.0);
        }
        std::vector<double> Mvec(Mn);
        size_t idx = 0;
        for(size_t i = this->data->forces.size(); i < external_forces.size(); ++i){
            if(external_forces[i].Mz){
                auto& ef_i = external_forces[i];
                for(size_t j = 0; j < external_forces.size(); ++j){
                    auto& ef_j = external_forces[j];
                    if(i != j){
                        gp_Vec r(ef_i.position, ef_j.position);
                        Mvec[idx] -= r.Crossed(ef_j.forces).Z();
                    }
                }
                ++idx;
            }
        }
        auto Mtemp(Mvec);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, Mn, Mn, 1, Mat.data(), Mn, Mtemp.data(), 1, 0, Mvec.data(), 1);
        size_t Mi = 0;
        for(size_t i = this->data->forces.size(); i < external_forces.size(); ++i){
            if(external_forces[i].Mz){
                external_forces[i].moments.SetZ(Mvec[Mi]);
                ++Mi;
            }
        }
    } else if(Mn == 1){
        for(size_t i = this->data->forces.size(); i < external_forces.size(); ++i){
            if(external_forces[i].Mz){
                auto& ef_i = external_forces[i];
                for(size_t j = 0; j < external_forces.size(); ++j){
                    auto& ef_j = external_forces[j];
                    gp_Vec r(gp_Pnt(0,0,0), ef_j.position);
                    ef_i.moments -= r.Crossed(ef_j.forces);
                }
                break;
            }
        }
    }

    // Calculate expansions
    std::vector<ExpansionNode> exp_info;
    exp_info.reserve(mesh.boundary_nodes.size() + external_forces.size());
    size_t ef_count = 0;
    for(auto& ef:external_forces){
        // logger::quick_log(ef.forces.X(), ef.forces.Y(), ef.moments.Z(), ef.position.X(), ef.position.Y());
        if(ef_count >= this->data->forces.size()){
            exp_info.push_back(get_expansion_node_2D(-ef.line_dir, ef.position, ef.diameter, ef.forces.X(), ef.forces.Y(), ef.moments.Z(), edges_init));
        } else {
            exp_info.push_back(get_expansion_node_2D(-ef.line_dir, ef.position, ef.diameter, ef.forces.X(), ef.forces.Y(), ef.moments.Z()));
        }
        ++ef_count;
    }

    struct SeparateNodes{
        gp_Pnt center;
        double distance;
        gp_Dir dir;
        TopoDS_Shape splitter;
    };
    std::vector<SeparateNodes> separate_analysis;

    size_t count = 0;
    for(auto& n:mesh.boundary_nodes){
        if(!this->is_valid_boundary_point(n.node)){
            ++count;
            continue;
        }

        // Get cross-section
        gp_Dir line_dir = -n.normal;
        gp_Lin line(n.node->point, line_dir);
        Handle(Geom_Line) geom_line(new Geom_Line(line));
        gp_Pnt opposite(Precision::Infinite(), Precision::Infinite(), Precision::Infinite());
        gp_Pnt center;
        double distance = Precision::Infinite();
        for(auto& e:edges_beam){
            double a = 0;
            double b = 0;
            Handle(Geom_Curve) edge_curve = BRep_Tool::Curve(e, a, b);
            GeomAPI_ExtremaCurveCurve extrema(geom_line, edge_curve, 0, Precision::Infinite(), a, b);
            if(extrema.NbExtrema() > 0){
                if(extrema.TotalLowerDistance() < Precision::Confusion()){
                    gp_Pnt p1, p2;
                    extrema.TotalNearestPoints(p1, p2);
                    double dist = n.node->point.Distance(p1);
                    if(dist > Precision::Confusion() && n.node->point.Distance(p1) < distance){
                        center = p1;
                        center.BaryCenter(1, n.node->point,1);
                        if(this->is_inside_2D(center, beams)){
                            opposite = p1;
                            distance = dist;
                        }
                    }
                }
            }
        }
        if(distance - Precision::Infinite() >= -Precision::Confusion()){
            for(TopExp_Explorer exp(beams, TopAbs_VERTEX); exp.More(); exp.Next()){
                gp_Pnt p = BRep_Tool::Pnt(TopoDS::Vertex(exp.Current()));
                if(line.Contains(p, Precision::Confusion())){
                    double dist = n.node->point.Distance(p);
                    if(dist > Precision::Confusion()){
                        if(n.node->point.Translated(dist*(line_dir)).IsEqual(p, Precision::Confusion()) && dist < distance){
                            center = n.node->point.Translated(dist*(line_dir));
                            center.BaryCenter(1, n.node->point,1);
                            if(this->is_inside_2D(center, beams)){
                                opposite = n.node->point.Translated(dist*(line_dir));
                                distance = dist;
                            }
                        }
                    }
                }
            }
        }

        // As it's redundant, we can afford to lose points that aren't working
        // correctly.
        if(opposite.IsEqual(gp_Pnt(Precision::Infinite(), Precision::Infinite(), Precision::Infinite()), Precision::Confusion())){
            continue;
        }

        // Get loads within cross-section
        TopoDS_Edge crosssection = BRepBuilderAPI_MakeEdge(n.node->point, opposite);
        gp_Dir nn(-n.normal.Y(), n.normal.X(), 0);
        TopoDS_Shape cs_face = BRepPrimAPI_MakePrism(crosssection, 0.1*gp_Vec(nn));
        double Fx = 0;
        double Fy = 0;
        double Mz = 0;
        GProp_GProps props;
        BRepGProp::SurfaceProperties(beams, props);
        double A = props.Mass();
        GProp_GProps props_cs;
        BRepGProp::SurfaceProperties(cs_face, props_cs);
        double A_cs = props_cs.Mass();

        BOPAlgo_Splitter splitter;
        TopoDS_Shape beams_copy = BRepBuilderAPI_Copy(beams);
        splitter.AddArgument(beams_copy);
        splitter.AddTool(cs_face);
        splitter.Perform();
        TopoDS_Shape result = splitter.Shape();
        TopoDS_Shape part;
        for(TopExp_Explorer shell_getter(result, TopAbs_FACE); shell_getter.More(); shell_getter.Next()){
            TopoDS_Shape tmp = shell_getter.Current();
            GProp_GProps props2;
            BRepGProp::SurfaceProperties(part, props2);
            double A2 = props2.Mass();
            if((A - A_cs) - A2 > Precision::Confusion() ){
                part = tmp;
                break;
            }
        }
        if(part.IsNull()){
            SeparateNodes sn{center, distance, line_dir, cs_face};
            separate_analysis.push_back(sn);
        }

        for(auto& ef:external_forces){
            if(this->is_inside_2D(ef.position, part)){
                gp_Vec r(center, ef.position);
                Fx -= ef.forces.X();
                Fy -= ef.forces.Y();
                Mz -= ef.moments.Z() + r.Crossed(ef.forces).Z();
            }
        }
        //logger::quick_log(Fx, Fy, Mz, center.X(), center.Y());

        // Calculate new cross-section
        exp_info.push_back(get_expansion_node_2D(n.normal, center, distance, Fx, Fy, Mz, edges_init));

        double pc = count/(double)(mesh.boundary_nodes.size()-1);
        std::cout << "\r" << pc*100 << "%         ";

        ++count;
    }

    // In case the faster expansion method is not enough, go beam by beam:
    // if(separate_analysis.size()){
    //     for(auto& beam:this->separate_beams){
    //         mesh.clear_results();
    //     }
    // }

    std::cout << "\r" << 100 << "%         ";
    std::cout << std::endl;

    logger::quick_log("Done.");

    logger::quick_log("Building geometry...");
    // Generate topology
    TopoDS_Shape result = BRepBuilderAPI_Copy(this->data->ground_structure->shape);
    result = utils::cut_shape(result, beams);

    std::cout << "0%";
    count = 0;
    for(auto& n:exp_info){
        double pc = count/(double)(exp_info.size()-1);
        if(n.diameter < Precision::Confusion()){
            std::cout << "\r" << pc*100 << "%         ";

            ++count;
            continue;
        }
        gp_Ax2 axis(n.center, gp_Dir(0,0,1));
        gp_Circ circ(axis, n.diameter/2);
        TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(circ);
        TopoDS_Wire wire = BRepBuilderAPI_MakeWire(edge);
        TopoDS_Face face = BRepBuilderAPI_MakeFace(wire);

        // Fastest method I managed to find, cut from original shape then
        // inverse cut back into original shape.
        // Still a bit slow though.
        result = utils::cut_shape(result, face);

        std::cout << "\r" << pc*100 << "%         ";

        ++count;
    }
    std::cout << "\r" << 100 << "%         ";
    std::cout << std::endl;

    result = utils::cut_shape(this->data->ground_structure->shape, result);
    result = this->simplify_shape(result);

    logger::quick_log("Done.");

    return result;
}

StandardSizing::ExpansionNode StandardSizing::get_expansion_node_2D(const gp_Dir& line_dir, gp_Pnt center, double distance, double Fx, double Fy, double Mz, const std::vector<TopoDS_Edge>& edges_init) const{
    gp_Vec normal(line_dir);
    normal.Rotate(gp_Ax1(center, gp_Dir(0,0,1)), M_PI/2);
    gp_Vec F(Fx, Fy, 0);

    double t = this->data->thickness;

    std::vector<double> S(this->data->material->get_max_stresses(normal));

    double S_f = std::min(S[0], S[1]);
    double S_n = (normal.Dot(F) < 0) ? S[0] : S[1];
    double S_c = S[2];

    // Bending
    double h_f = std::sqrt(6*std::abs(Mz)/(t*S_f));
    // Normal
    double h_n = std::abs(normal.Dot(F))/(t*S_n);
    // Shear
    double h_c = (F - normal.Dot(F)*normal).Magnitude()*(3/(2*t*S_c));

    double h = this->multiplier*std::max({h_f, h_n, h_c});//, distance});
    if(h - distance < Precision::Confusion()){
        return {center, line_dir, h};
    }

    // Calculate new position for neutral axis, if desired
    if(edges_init.size() > 0){
        gp_Lin side1lin(center, line_dir);
        gp_Lin side2lin(center, -line_dir);
        Handle(Geom_Line) geom_line1(new Geom_Line(side1lin));
        Handle(Geom_Line) geom_line2(new Geom_Line(side2lin));
        
        double min_dist1 = Precision::Infinite();
        double min_dist2 = Precision::Infinite();
        for(auto& e:edges_init){
            BRepExtrema_DistShapeShape edge_checker(e, BRepBuilderAPI_MakeVertex(center));
            edge_checker.Perform();
            if(edge_checker.Value() < Precision::Confusion()){
                continue;
            }

            double a = 0;
            double b = 0;
            Handle(Geom_Curve) edge_curve = BRep_Tool::Curve(e, a, b);
            GeomAPI_ExtremaCurveCurve extrema1(geom_line1, edge_curve, 0, h/2, a, b);
            GeomAPI_ExtremaCurveCurve extrema2(geom_line2, edge_curve, 0, h/2, a, b);
            if(extrema1.NbExtrema() > 0){
                if(extrema1.TotalLowerDistance() < Precision::Confusion()){
                    gp_Pnt p1, p2;
                    extrema1.TotalNearestPoints(p1, p2);
                    double dist = center.Distance(p1);
                    if(dist > Precision::Confusion() && dist < min_dist1){
                        min_dist1 = dist;
                    }
                }
            }
            if(extrema2.NbExtrema() > 0){
                if(extrema2.TotalLowerDistance() < Precision::Confusion()){
                    gp_Pnt p1, p2;
                    extrema2.TotalNearestPoints(p1, p2);
                    double dist = center.Distance(p1);
                    if(dist > Precision::Confusion() && dist < min_dist2){
                        min_dist2 = dist;
                    }
                }
            }
        }
        if(min_dist1 < Precision::Infinite() - Precision::Confusion() && min_dist2 < Precision::Infinite() - Precision::Confusion()){
            gp_Pnt p1 = center.Translated(min_dist1*line_dir);
            gp_Pnt p2 = center.Translated(min_dist2*(-line_dir));
            center = p1;
            center.BaryCenter(1, p2, 1);
        } else if(min_dist1 < Precision::Infinite() - Precision::Confusion()){
            double dev = std::min(h/2-min_dist1, h/2-distance/2);
            center.Translate(-dev*line_dir);
        } else if(min_dist2 < Precision::Infinite() - Precision::Confusion()){
            double dev = std::min(h/2-min_dist2, h/2-distance/2);
            center.Translate(-dev*(-line_dir));
        }
    }

    return {center, line_dir, h};
}

bool StandardSizing::is_valid_boundary_point(MeshNode* n) const{
    for(auto& f:this->data->forces){
        if(f.S.is_inside(n->point)){
            return false;
        }
    }
    for(auto& s:this->data->supports){
        if(s.S.is_inside(n->point)){
            return false;
        }
    }
    return true;
}

TopoDS_Shape StandardSizing::simplify_shape(TopoDS_Shape shape) const{
    double tol = 5;
    double prec = tol*1.5;

    Handle(ShapeFix_Shape) sfs = new ShapeFix_Shape(shape);
    sfs->SetPrecision(prec);
    sfs->SetMaxTolerance(2*tol);
    sfs->SetMinTolerance(tol);
    auto sfw = sfs->FixWireTool();
    sfw->ModifyGeometryMode() = false;
    sfw->ModifyTopologyMode() = false;
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
            TopoDS_Shape csec = BRepBuilderAPI_Copy(f.S.get_shape());
            gp_Vec transl(f.S.get_centroid(), b.front());
            if(transl.Magnitude() > Precision::Confusion()){
                gp_Trsf trsf;
                trsf.SetTranslation(transl);
                csec = BRepBuilderAPI_Transform(csec, trsf);
            }
            gp_Dir fn = f.S.get_normal();
            gp_Dir sn(gp_Vec(b[0], b[1]));
            if(!fn.IsEqual(sn, Precision::Confusion()) && !fn.IsEqual(-sn, Precision::Confusion())){
                gp_Dir axis = fn.Crossed(sn);
                double ang = fn.AngleWithRef(sn, axis);
                gp_Ax1 ax(b.front(), axis);
                gp_Trsf trsf;
                trsf.SetRotation(ax, ang);
                csec = BRepBuilderAPI_Transform(csec, trsf);
            }

            TopoDS_Shape beam = utils::sweep_surface(b, csec, data->ground_structure->shape);
            this->separate_beams.push_back(beam);
            // TopoDS_Shape beam = utils::fast_make_2D_beam(b, f.S.get_dimension(), data->ground_structure->shape);
            geometry = utils::cut_shape(geometry, beam);
        }
    }

    return utils::cut_shape(this->data->ground_structure->shape, geometry);
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
            TopoDS_Edge line_edge = BRepBuilderAPI_MakeEdge(line, -Precision::Infinite(), Precision::Infinite());
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
                        if(dist >= 0 && dist < min_dist){
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
                        if(dist > Precision::Confusion() && dist < min_dist){
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
    h2 = Km.multiply(h2);
    std::vector<size_t> affected = Km.affected_ids(boundary_ids);
    for(size_t i = 0; i < affected.size(); ++i){
        if(std::abs(h2[affected[i]]) < std::abs(h[affected[i]])){
            h[affected[i]] = h2[affected[i]];
        }
    }
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

TopoDS_Shape StandardSizing::experimental_elemental_approach(){
    double mesh_size = 4;

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
                // std::vector<double> forces1 = e->get_loads_at(n1->point, u);
                // std::vector<double> forces2 = e->get_loads_at(n2->point, u);
                std::vector<double> forces = e->get_average_loads_and_torques(n2->point, n1->point, u);
                // std::vector<double> forces = e->get_average_loads_and_torques(n2->point, n1->point, u);
                // if((forces1[0]+forces2[0]) == 0 && (forces1[1]+forces2[1]) == 0){
                //     continue;
                // }
                // gp_Vec F(forces1[0]+forces2[0], forces1[1]+forces2[1], 0);
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
                logger::quick_log(F.X(), F.Y());
                if(uh[2*n1->id] > -1){
                    // int sign = dir.X()/std::abs(dir.X());
                    // h[uh[2*n1->id]] += sign*std::max(0.0, std::abs(F.Y())/(t*S_x) - std::abs(n1->point.X() - n2->point.X()));
                    h[uh[2*n1->id]] += dir.X()*hh;
                }
                if(uh[2*n1->id+1] > -1){
                    // int sign = dir.Y()/std::abs(dir.Y());
                    // h[uh[2*n1->id+1]] += sign*std::max(0.0, std::abs(F.X())/(t*S_x) - std::abs(n1->point.Y() - n2->point.Y()));
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
            
            result = utils::cut_shape(result, s);

            double pc = i/(double)(mesh.element_list.size()-1);
            std::cout << "\r" << pc*100 << "%         ";

            ++i;
        }
        std::cout << "\r" << 100 << "%         ";
        std::cout << std::endl;

        logger::quick_log("Done.");
    }
    result = utils::cut_shape(this->data->ground_structure->shape, result);

    return result;
}

bool StandardSizing::is_inside_2D(const gp_Pnt& p, const TopoDS_Shape& shape) const{
    BRepClass3d_SolidClassifier insider(shape, p, 0.1);
    return insider.State() == TopAbs_ON;
}

bool StandardSizing::is_inside_3D(const gp_Pnt& p, const TopoDS_Shape& shape) const{
    BRepClass3d_SolidClassifier insider(shape, p, 0.1);
    return insider.State() == TopAbs_IN;
}

}