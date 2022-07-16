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
#include <BRepAlgoAPI_Section.hxx>
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
#include <Geom_Circle.hxx>
#include <ShapeCustom.hxx>
#include <ShapeCustom_RestrictionParameters.hxx>
#include "element/GT9.hpp"

namespace sizing{

StandardSizing::StandardSizing(ProjectData* data, FiniteElement* solver, double element_size, double multiplier):
   Sizing(data), solver(solver), element_size(element_size), multiplier(multiplier){}

TopoDS_Shape StandardSizing::run(){
    return this->boundary_expansion_approach();
}

TopoDS_Shape StandardSizing::boundary_expansion_approach(){

    // Get beams and do FEA
    meshing::StandardBeamMesher mesh(this->element_size, 1, utils::PROBLEM_TYPE_2D, this->data);
    TopoDS_Shape beams = this->build_initial_topology();
    utils::shape_to_file("beams.step", beams);
    auto m = mesh.mesh(beams);
    std::unique_ptr<MeshElementFactory> gt9_maker(new MeshElementFactoryImpl<element::GT9>());
    mesh.prepare_for_FEM(m, gt9_maker, this->data);//, true);
    auto u = this->solver->calculate_displacements(this->data, &mesh);

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

        auto reactions = this->solver->calculate_forces(&mesh, u, this->data->topopt_element);
        auto dof = this->data->topopt_element->get_dof_per_node();
        double Fx = 0;
        double Fy = 0;
        double Mz = 0;
        for(auto& n:mesh.node_list){
            if(final_line.Contains(n->point, Precision::Confusion())){
                double dist = center.Distance(n->point);
                if(dist < distance/2 + Precision::Confusion()){
                    Fx += reactions[n->id*dof + 0];
                    Fy += reactions[n->id*dof + 1];
                    // This method is way too good to be true.
                    // But looks like it's true.
                    gp_Vec F(reactions[n->id*dof + 0], reactions[n->id*dof + 1], 0);
                    //logger::quick_log(F.X(), F.Y());
                    gp_Vec d(center, n->point);
                    //logger::quick_log(d.X(), d.Y());
                    gp_Vec M(d.Crossed(F));
                    Mz += M.Z();

                }
            }
        }
        logger::quick_log(Fx, Fy, Mz);
        ExternalForce ef;
        ef.forces = gp_Vec(Fx, Fy, 0);
        ef.moments = gp_Vec(0, 0, Mz);
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
    // this->calculate_reaction_moments(Mn, external_forces);

    // Calculate expansions
    std::vector<ExpansionNode> exp_info;
    std::vector<ExpansionNode> exp_tmp;
    exp_info.reserve(mesh.boundary_nodes.size() + external_forces.size());

    double min_diam = Precision::Infinite();
    double max_diam = 0;
    for(auto& f:this->data->forces){
        double d = f.S.get_dimension();
        if(d < min_diam){
            min_diam = d;
        }
        max_diam += d;
    }

    size_t count = 0;
    std::vector<BeamMeshing::BoundaryNode> bnodes = mesh.boundary_nodes;
    for(size_t i = 0; i < bnodes.size(); ++i){
        auto& n = bnodes[i];
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
        size_t id = 0;
        double distance = Precision::Infinite();
        for(size_t j = i+1; j < bnodes.size(); ++j){
            gp_Pnt p = bnodes[j].node->point;
            if(line.Contains(p, Precision::Confusion())){
                double dist = n.node->point.Distance(p);
                // using max_diam is a bit of a heuristic checking, as sometimes the algorithm
                // ends up picking a point on the other side of the geometry, for some reason.
                if(dist > Precision::Confusion() && dist < distance && dist < max_diam + Precision::Confusion()){
                    if(n.node->point.Translated(dist*(line_dir)).IsEqual(p, Precision::Confusion())){
                        gp_Pnt c = n.node->point.Translated(dist*(line_dir));
                        c.BaryCenter(1, n.node->point,1);
                        if(this->is_inside_2D(c, beams)){
                            opposite = p;//n.node->point.Translated(dist*(line_dir));
                            distance = dist;
                            center = c;
                            id = j;
                        }
                    }
                }
            }
        }
        if(id > 0){
            bnodes.erase(bnodes.begin()+id);
        } else {
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
                        if(dist > Precision::Confusion() && dist < distance){// && this->is_inside_2D(p1, beams)){
                            gp_Pnt c = p1;
                            c.BaryCenter(1, n.node->point, 1);
                            if(this->is_inside_2D(c, beams)){
                                opposite = p1;
                                distance = dist;
                                center = c;
                                continue;
                            }
                        }
                    }
                }
            }
        }

        // Get loads within cross-section
        TopoDS_Edge crosssection = BRepBuilderAPI_MakeEdge(n.node->point, opposite);
        gp_Vec nn(-n.normal.Y(), n.normal.X(), 0);
        nn.Normalize();

        std::vector<IntersectionNode> int_nodes;
        for(auto& e:mesh.element_list){
            // Heuristic filter
            TopoDS_Vertex v = BRepBuilderAPI_MakeVertex(e->get_centroid());
            auto extrema = BRepExtrema_DistShapeShape(crosssection, v, Extrema_ExtFlag_MIN);
            double dist = extrema.Value();
            if(dist > this->element_size){
                continue;
            }

            std::vector<gp_Pnt> points = e->get_intersection_points(crosssection);
            if(points.size() == 0){
                continue;
            }

            std::vector<double> force = e->get_internal_loads(u);
            size_t dof = e->get_element_info()->get_dof_per_node();
            for(size_t j = 0; j < e->nodes.size(); ++j){
                auto& ne = e->nodes[j];
                double line_pos_rel = nn.Dot(gp_Vec(ne->point, n.node->point));
                // Check for nodes that are on the same side
                if(line_pos_rel > Precision::Confusion()){
                    int_nodes.push_back({
                        ne->point,
                        gp_Vec(force[j*dof+0], force[j*dof+1], 0)
                    });
                }
            }
        }

        double Fx = 0;
        double Fy = 0;
        double Mz = 0;
        for(auto& in:int_nodes){
            Fx += in.force.X();
            Fy += in.force.Y();
            gp_Vec r(center, in.point);
            Mz += r.Crossed(in.force).Z();
        }

        auto node = get_expansion_node_2D(n.normal, center, distance, Fx, Fy, Mz, edges_init);

        if(node.diameter > distance){
            exp_info.push_back(node);
        }


        double pc = count/(double)(bnodes.size()-1);
        std::cout << "\r" << pc*100 << "%         ";

        ++count;
    }

    // Add external forces
    size_t ef_count = 0;
    for(auto& ef:external_forces){
        if(ef_count >= this->data->forces.size()){
            auto node = get_expansion_node_2D(-ef.line_dir, ef.position, ef.diameter, ef.forces.X(), ef.forces.Y(), ef.moments.Z(), edges_init);
            exp_info.push_back(node);
        } else {
            auto node = get_expansion_node_2D(-ef.line_dir, ef.position, ef.diameter, ef.forces.X(), ef.forces.Y(), ef.moments.Z());
            exp_info.push_back(node);
        }
        ++ef_count;
    }

    std::cout << "\r" << 100 << "%         ";
    std::cout << std::endl;

    logger::quick_log("Done.");

    logger::quick_log("Building geometry...");
    // Generate topology
    TopoDS_Shape result = BRepBuilderAPI_Copy(this->data->ground_structure->shape);
    result = utils::cut_shape(result, beams);

    // TODO: make this work
    //return this->bspline_simple2D(exp_info, result);
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


void StandardSizing::calculate_reaction_moments(size_t Mn, std::vector<ExternalForce>& external_forces) const{
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
}

StandardSizing::ExpansionNode StandardSizing::get_expansion_node_2D(const gp_Dir& line_dir, gp_Pnt center, double distance, double Fx, double Fy, double Mz, const std::vector<TopoDS_Edge>& edges_init) const{
    gp_Dir new_dir = line_dir;
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
        double min_dist1 = Precision::Infinite();
        double min_dist2 = Precision::Infinite();
        double norm1 = 0;
        double norm2 = 0;
        gp_Pnt p1(0,0,0);
        gp_Pnt p2(0,0,0);
        gp_Dir dir1(0,0,1);
        gp_Dir dir2(0,0,1);
        for(auto& e:edges_init){
            BRepExtrema_DistShapeShape edge_checker(BRepBuilderAPI_MakeVertex(center), e);
            edge_checker.Perform();
            double dist = edge_checker.Value();
            if(dist < Precision::Confusion() || dist > h/2){
                continue;
            }
            gp_Pnt p;
            for(int i = 1; i <= edge_checker.NbSolution(); ++i){
                p = edge_checker.PointOnShape2(i);
                if(center.Distance(p) <= dist + Precision::Confusion()){
                    break;
                }
            }
            gp_Vec vec(center, p);
            gp_Vec line_vec(line_dir);
            //gp_Vec translvec((h/2-dist)*(-vec.Normalized()));
            gp_Vec translvec(line_dir);

            if(vec.IsNormal(line_vec, Precision::Angular())){
                continue;
            }
            double b = dist;
            double c = h/2;
            double gamma = vec.Angle(line_dir);
            double gamma2 = vec.Angle(-line_dir);
            if(gamma2 > gamma){
                gamma = gamma2;
                translvec = -line_dir;
            }
            double delta = std::pow(2*b*std::cos(gamma), 2) - 4*(b*b-c*c);
            double a = (2*b*std::cos(gamma) + std::sqrt(delta))/2;
            a = std::max(a, (2*b*std::cos(gamma) - std::sqrt(delta))/2);
            double norm = std::abs(a);
            translvec.Multiply(a);

            if(line_dir.X()*translvec.X() >= -Precision::Confusion() && line_dir.Y()*translvec.Y() >= -Precision::Confusion() && norm > norm1){
                p1 = p;
                norm1 = norm;
            } else if(line_dir.X()*translvec.X() <= Precision::Confusion() && line_dir.Y()*translvec.Y() <= Precision::Confusion() && norm > norm2){
                p2 = p;
                norm2 = norm;
            }

        }
        if(norm1 > Precision::Confusion() && norm2 > Precision::Confusion()){
            gp_Pnt pp1 = center.Translated(norm1*line_dir);
            gp_Pnt pp2 = center.Translated(norm2*(-line_dir));
            center = pp1;
            center.BaryCenter(1, pp2, 1);
        } else if(norm1 > Precision::Confusion()){
            center.Translate(norm1*line_dir);
            new_dir = gp_Vec(center, p1);
        } else if(norm2 > Precision::Confusion()){
            center.Translate(norm2*(-line_dir));
            new_dir = -gp_Vec(center, p2);
        }
    }

    return {center, new_dir, h};
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
    double tol = 20;
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

    // Handle(ShapeCustom_RestrictionParameters) rp(new ShapeCustom_RestrictionParameters());
    // shape = ShapeCustom::BSplineRestriction(shape,
    //         tol,
    //         tol,
    //         200,
    //         400,
    //         GeomAbs_C2,
    //         GeomAbs_C2,
    //         true,
    //         true,
    //         rp);

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

bool StandardSizing::is_inside_2D(const gp_Pnt& p, const TopoDS_Shape& shape) const{
    BRepClass3d_SolidClassifier insider(shape, p, 0.1);
    return insider.State() == TopAbs_ON;
}

bool StandardSizing::is_inside_3D(const gp_Pnt& p, const TopoDS_Shape& shape) const{
    BRepClass3d_SolidClassifier insider(shape, p, 0.1);
    return insider.State() == TopAbs_IN;
}


TopoDS_Shape StandardSizing::bspline_simple2D(const std::vector<ExpansionNode>& exp_info, TopoDS_Shape base) const{
    TopoDS_Shape result = base;

    // gp_Pnt prev_p2 = exp_info[0].center.Translated( (exp_info[0].diameter/2)*exp_info[0].direction);
    // gp_Pnt prev_p3 = exp_info[0].center.Translated(-(exp_info[0].diameter/2)*exp_info[0].direction);
    // TopoDS_Vertex prev_v2 = BRepBuilderAPI_MakeVertex(prev_p2);
    // TopoDS_Vertex prev_v3 = BRepBuilderAPI_MakeVertex(prev_p3);
    std::cout << "0%";
    int count = 0;
    logger::quick_log(exp_info[0].center.X(), exp_info[0].center.Y(), exp_info[0].diameter);
    for(size_t i = 1; i < exp_info.size(); ++i){
        double pc = count/(double)(exp_info.size()-1);
        if(exp_info[i-1].diameter < Precision::Confusion() || exp_info[i].diameter < Precision::Confusion() 
            || exp_info[i-1].center.IsEqual(exp_info[i].center, Precision::Confusion())){
            std::cout << "\r" << pc*100 << "%         ";

            ++count;
            continue;
        }
        logger::quick_log(exp_info[i].center.X(), exp_info[i].center.Y(), exp_info[i].diameter);
        gp_Dir dir1 = exp_info[i-1].direction;
        gp_Dir dir2 = exp_info[i].direction;
        if(dir1.X()*dir2.X() < Precision::Confusion() && dir1.Y()*dir2.Y() < Precision::Confusion()){
            dir2 = -dir2;
        }
        gp_Pnt p1 = exp_info[i-1].center.Translated( (exp_info[i-1].diameter/2)*dir1);
        gp_Pnt p4 = exp_info[i-1].center.Translated(-(exp_info[i-1].diameter/2)*dir1);
        TopoDS_Vertex v1 = BRepBuilderAPI_MakeVertex(p1);
        TopoDS_Vertex v4 = BRepBuilderAPI_MakeVertex(p4);
        // TopoDS_Vertex v1(std::move(prev_v2));
        // TopoDS_Vertex v4(std::move(prev_v3));

        gp_Pnt p2 = exp_info[i].center.Translated( (exp_info[i].diameter/2)*dir2);
        gp_Pnt p3 = exp_info[i].center.Translated(-(exp_info[i].diameter/2)*dir2);
        TopoDS_Vertex v2 = BRepBuilderAPI_MakeVertex(p2);
        TopoDS_Vertex v3 = BRepBuilderAPI_MakeVertex(p3);

        TopoDS_Edge e1 = BRepBuilderAPI_MakeEdge(v1, v2);
        TopoDS_Edge e2 = BRepBuilderAPI_MakeEdge(v2, v3);
        TopoDS_Edge e3 = BRepBuilderAPI_MakeEdge(v3, v4);
        TopoDS_Edge e4 = BRepBuilderAPI_MakeEdge(v4, v1);
        TopoDS_Wire w = BRepBuilderAPI_MakeWire(e1, e2, e3, e4);
        TopoDS_Face f = BRepBuilderAPI_MakeFace(w);

        result = utils::cut_shape(result, f);

        // prev_v2 = std::move(v2);
        // prev_v3 = std::move(v3);

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

bool StandardSizing::insert_expansion_node(std::vector<ExpansionNode>& exp_info, ExpansionNode node) const{
    auto angle = [](const gp_Dir& d1, const gp_Dir& d2)->double{
        double dot = std::abs(d1.Dot(d2));

        return std::acos(dot);
    };

    bool inserted = false;
    if(exp_info.size() < 2){
        inserted = true;
        exp_info.push_back(node);
    } else {
        // Directional sort?
        gp_Pnt p(node.center);
        gp_Dir n(node.direction);
        for(auto i = exp_info.begin()+1; i < exp_info.end(); ++i){
            gp_Pnt p1((i-1)->center);
            gp_Pnt p2(i->center);
            gp_Vec dir1(p1, p);
            gp_Vec dir2(p2, p);
            gp_Dir n1((i-1)->direction);
            gp_Dir n2(i->direction);
            if(p2.IsEqual(p, Precision::Confusion())){
                if(std::abs(i->diameter - node.diameter) <= Precision::Confusion()
                  && (n2.IsEqual(n, Precision::Confusion()) || n2.IsEqual(-n, Precision::Confusion()))){
                    inserted = true;
                    break;
                } else {
                    continue;
                }
            }
            if(p1.IsEqual(p, Precision::Confusion())){
                if(std::abs((i-1)->diameter - node.diameter) <= Precision::Confusion()
                  && (n1.IsEqual(n, Precision::Confusion()) || n1.IsEqual(-n, Precision::Confusion()))){
                    inserted = true;
                    break;
                } else {
                    continue;
                }
            }
            bool is_between = dir1.X()*dir2.X() < Precision::Confusion() && dir1.Y()*dir2.Y() < Precision::Confusion();
            bool is_aligned = std::abs(angle(n, n1) + angle(n, n2) - angle(n1, n2)) < Precision::Confusion();//compare line dirs
            if(is_between && is_aligned){
                //if(node.diameter + Precision::Confusion() >= std::min(i->diameter, (i-1)->diameter)){
                    inserted = true;
                    exp_info.insert(i, node);
                //}
                break;
            }
        }
        if(!inserted){
            if(node.center.Distance(exp_info.front().center) < node.center.Distance(exp_info.back().center)){
                // gp_Vec dir1(exp_info[0].center, p);
                // gp_Vec dir2(exp_info[1].center, exp_info[0].center);
                // bool is_before = dir1.X()*dir2.X() > -Precision::Confusion() && dir1.Y()*dir2.Y() > -Precision::Confusion();
                // if(is_before){ // Otherwise, it's probably a miscalculated node.
                    exp_info.insert(exp_info.begin(), node);
                    inserted = true;
                // }
            } else {
                // gp_Vec dir1(exp_info[exp_info.size()-1].center, p);
                // gp_Vec dir2(exp_info[exp_info.size()-2].center, exp_info[exp_info.size()-1].center);
                // bool is_after = dir1.X()*dir2.X() > -Precision::Confusion() && dir1.Y()*dir2.Y() > -Precision::Confusion();
                // if(is_after){ // Otherwise, it's probably a miscalculated node.
                    exp_info.push_back(node);
                    inserted = true;
                // }
            }
        }
    }

    return inserted;
}

}
