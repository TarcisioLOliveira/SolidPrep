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

#include "optimizer.hpp"
#include "function.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include <BRepBuilderAPI_Transform.hxx>
#include <GProp_GProps.hxx>
#include <Interface_Static.hxx>
#include <STEPControl_Reader.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <TopoDS_Builder.hxx>

Constraint::Constraint(std::unique_ptr<DensityBasedFunction> fun, std::vector<Type> types, std::vector<double> bounds):
    fun(std::move(fun)), types(std::move(types)), bounds(std::move(bounds)){}

void Optimizer::initialize_optimizer(const Meshing* const mesh){
    this->number_of_elements = this->get_number_of_elements(mesh->geometries);
    this->volumes.resize(this->number_of_elements);
    this->stresses.resize(this->number_of_elements);

    this->get_volumes(mesh->geometries, mesh->thickness, this->volumes);
}

TopoDS_Shape Optimizer::STEP_workaround(const TopoDS_Shape& s) const{
    STEPControl_Writer writer;
    STEPControl_StepModelType mode = STEPControl_AsIs;
    Interface_Static::SetIVal("write.surfacecurve.mode",0);
    writer.PrintStatsTransfer(0, 0);
    writer.PrintStatsTransfer(1, 0);
    writer.PrintStatsTransfer(2, 0);
    writer.PrintStatsTransfer(3, 0);
    writer.PrintStatsTransfer(4, 0);
    writer.PrintStatsTransfer(5, 0);
    IFSelect_ReturnStatus stat = writer.Transfer(s,mode);
    (void) stat;
    auto ws = writer.WS();
    STEPControl_Reader reader(ws, false);

    Standard_Integer NbRoots = reader.NbRootsForTransfer();
    Standard_Integer num = reader.TransferRoots();
    (void) NbRoots;
    (void) num;
    return reader.OneShape();
}

TopoDS_Shape Optimizer::make_shape(const std::vector<double>& x, const std::vector<Geometry*>& geometries, const double result_threshold, const utils::ProblemType type) const{
    logger::quick_log(" "); 
    logger::quick_log("Saving resulting geometries...");
    // TODO: make it work with multiple geometries
    std::cout << "\r" << 0 << "%         ";

    if(type == utils::PROBLEM_TYPE_2D){
        TopoDS_Shape result = BRepBuilderAPI_Copy(geometries[0]->shape);
        for(size_t i = 0; i < x.size(); ++i){
            if(x[i] < result_threshold){
                TopoDS_Shape s(geometries[0]->mesh[i]->get_shape());
                result = utils::cut_shape(result, s);
            }
            const double pc = i/(double)(x.size()-1);
            std::cout << "\r" << pc*100 << "%         ";
        }
        std::cout << "\r" << 100 << "%         ";
        logger::quick_log(" "); 
        return result;
    } else if(type == utils::PROBLEM_TYPE_3D){
        TopoDS_Compound result;
        TopoDS_Builder builder;
        builder.MakeCompound(result);
        size_t g_num = 0;
        size_t count = 0;
        auto x_it = x.cbegin();
        for(auto& g:geometries){
            if(g->do_topopt){
                TopoDS_Compound result_geom;
                TopoDS_Builder builder_geom;
                builder_geom.MakeCompound(result_geom);
                for(size_t i = 0; i < g->mesh.size(); ++i){
                    if(*x_it >= result_threshold){
                        TopoDS_Shape s(g->mesh[i]->get_shape());
                        builder_geom.Add(result_geom, std::move(s));
                        // auto ss = STEP_workaround(s);

                        // // result = utils::cut_shape(result, ss);
                        // TopoDS_Shape cur_shape;
                        // double cur_mass = 0;

                        // TopExp_Explorer exp(result, TopAbs_SHAPE);
                        // if(exp.More()){
                        //     for(; exp.More(); exp.Next()){
                        //         auto cur = exp.Current();
                        //         GProp_GProps props;
                        //         BRepGProp::SurfaceProperties(cur, props);
                        //         if(props.Mass() > cur_mass){
                        //             cur_mass = props.Mass();
                        //           cur_shape = cur;
                        //         }
                        //         break;
                        //     }
                        //     auto test = utils::cut_shape(cur_shape, ss);
                        //     if(test.IsNull()){
                        //         result = utils::cut_shape(result, ss);
                        //     } else {
                        //         result = test;
                        //     }
                        // } else {
                        //     result = utils::cut_shape(result, ss);
                        // }
                    }
                    ++x_it;
                    const double pc = count/(double)(x.size()-1);
                    std::cout << "\r" << pc*100 << "%         ";
                    ++count;
                }
                std::cout << std::endl;
                utils::shape_to_file("result_"+std::to_string(g_num)+".step", result_geom);
                ++g_num;
                //builder.Add(result, result_geom);
            }
        }

        // // Another workaround.
        // // When using solids, BRepAlgoAPI_Cut sometimes only cuts the faces,
        // // leaving the solid part in place, but separate from the remaining
        // // solid. Therefore, at the end of the process, extract the result
        // // geometry by finding the solid with the greatest mass (assuming
        // // it remains completely united, of course).
        // TopoDS_Shape cur_shape;
        // double cur_mass = 0;

        // for(TopExp_Explorer exp(result, TopAbs_SOLID); exp.More(); exp.Next()){
        //     auto cur = exp.Current();
        //     GProp_GProps props;
        //     BRepGProp::VolumeProperties(cur, props);
        //     if(props.Mass() > cur_mass){
        //         cur_mass = props.Mass();
        //         cur_shape = cur;
        //     }
        // }
        // result = cur_shape;
        std::cout << "\r" << 100 << "%         ";
        logger::quick_log(" "); 
        return result;
    }
    return TopoDS_Shape();
}

void Optimizer::get_stresses(const std::vector<Geometry*> geometries, const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& stresses, double pc, double psi) const{
    (void)x;
    auto stress_it = stresses.begin();
    auto rho_it = x.begin();
    for(const auto& g:geometries){
        g->get_stresses(u, pc, psi, rho_it, stress_it);
    }
}

void Optimizer::get_volumes(const std::vector<Geometry*> geometries, const double thickness, std::vector<double>& volumes) const{
    auto V_it = volumes.begin();
    for(const auto& g:geometries){
        for(const auto& e:g->mesh){
            *V_it = e->get_volume(thickness);
            ++V_it;
        }
    }
}

size_t Optimizer::get_number_of_elements(const std::vector<Geometry*> geometries) const{
    size_t elem_num = 0;
    for(const auto& g:geometries){
        elem_num += g->mesh.size();
    }

    return elem_num;
}

void Optimizer::apply_densities(const std::vector<Geometry*> geometries, const std::vector<double>& x, std::vector<double>& vals, const double pc) const{
    if(pc == 0){
        return;
    } else if(pc == 1){
        auto x_it = x.cbegin();
        auto v_it = vals.begin();
        for(const auto& g:geometries){
            if(g->do_topopt){
                for(auto xi = x_it; xi < x_it + g->mesh.size(); ++xi, ++v_it){
                    *v_it *= *xi;
                }
                x_it += g->mesh.size();
            } else {
                v_it += g->mesh.size();
            }
        }
    } else {
        auto x_it = x.cbegin();
        auto v_it = vals.begin();
        for(const auto& g:geometries){
            if(g->do_topopt){
                for(auto xi = x_it; xi < x_it + g->mesh.size(); ++xi, ++v_it){
                    *v_it *= std::pow(*xi, pc);
                }
                x_it += g->mesh.size();
            } else {
                v_it += g->mesh.size();
            }
        }
    }
}
