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

#include "meshing/gmsh.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include <BRepBuilderAPI_Transform.hxx>
#include <BRep_Builder.hxx>
#include <gmsh.h>
#include <algorithm>
#include "project_data.hpp"
#include <iterator>
#include <unordered_map>
#include <BOPAlgo_Splitter.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <TopoDS.hxx>
#include <set>
#include <vector>

namespace meshing{

BoundaryRefinement::BoundaryRefinement(const projspec::DataMap* const data){
    if(data != nullptr){
        this->min_dist = data->get_double("min_dist");
        this->max_dist = data->get_double("max_dist");
        this->min_size = data->get_double("min_size");
        this->root_op = data->get_shape_op("region");
    }
}

std::vector<double> BoundaryRefinement::get_faces(const gmsh::vectorpair& geoms) const{
    auto tmp = this->get_faces_int(this->root_op.get(), geoms);
    std::vector<double> conv(tmp.size());
    std::copy(tmp.begin(), tmp.end(), conv.begin());

    return conv;
}

std::vector<int> BoundaryRefinement::get_faces_int(shape_op::ShapeOp* op, const gmsh::vectorpair& geoms) const{
    switch(op->get_type()){
        case shape_op::Code::GEOMETRY:{
            const size_t id = op->get_id();
            logger::log_assert(id < geoms.size(), logger::ERROR,
                    "unknown geometry id: {}", id);
            std::vector<int> up, s;
            gmsh::model::getAdjacencies(geoms[id].first, geoms[id].second, up, s);
            std::sort(s.begin(), s.end());

            return s;
        }
        case shape_op::Code::UNION:{
            auto s1 = this->get_faces_int(op->first(), geoms);
            auto s2 = this->get_faces_int(op->second(), geoms);
            s1.insert(s1.end(), s2.begin(), s2.end());
            s2.clear();
            std::sort(s1.begin(), s1.end());
            return s1;
        }
        case shape_op::Code::INTERSECTION:{
            auto s1 = this->get_faces_int(op->first(), geoms);
            auto s2 = this->get_faces_int(op->second(), geoms);
            std::vector<int> result;
            std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(result));
            s1.clear();
            s2.clear();
            return result;
        }
        case shape_op::Code::DIFFERENCE:{
            auto s1 = this->get_faces_int(op->first(), geoms);
            auto s2 = this->get_faces_int(op->second(), geoms);
            std::vector<int> result;
            std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(result));
            s1.clear();
            s2.clear();
            return result;
        }
        case shape_op::Code::SHELL:{
            logger::log_assert(false, logger::ERROR,
                    "SHELL shape operation is currently not implemented");
            return std::vector<int>();
        }
    }

    return std::vector<int>();
}

Gmsh::Gmsh(const projspec::DataMap& data):
    Meshing(data.proj->geometries, 
            data.proj->topopt_element,
            data.proj,
            data.proj->thickness),
    tmp_scale(data.get_double("tmp_scale", 1.0)),
    size(tmp_scale*data.get_double("element_size")),
    size_from_curvature(data.get_double("size_from_curvature", 0)),
    algorithm2D(data.get_int("algorithm2D", 6)),
    algorithm3D(data.get_int("algorithm3D", 4)),
    bound_ref(data.get_object("boundary_refinement"))    
{


}

void Gmsh::mesh(const std::vector<Force>& forces, 
                const std::vector<Support>& supports,
                std::vector<Spring>& springs){
    (void)springs;
    TopoDS_Shape shape = this->make_compound(this->geometries);

    bool has_condition_inside = false;
    const auto problem_type = this->elem_info->get_problem_type();
    TopoDS_Shape sh;
    if(this->tmp_scale == 1.0){
        sh = BRepBuilderAPI_Copy(shape);
    } else {
        sh = this->make_compound(this->geometries, this->tmp_scale);
    }
    if(problem_type == utils::PROBLEM_TYPE_2D){
        // I'll need to find out later how to deal with this.
        // Boundary conditions inside the domain don't make a lot of sense,
        // but they're useful for design automation.
        // I'm currently phasing this out because now I just duplicate
        // the boundary in these cases anyway, but I still don't know
        // if the behavior is exactly the same.
        has_condition_inside = this->adapt_for_boundary_condition_inside(shape, forces, supports);
    } else {
        // Maybe disable for cubic elements?
        for(auto& f:forces){
            if(f.cut_into_shape){
                //has_condition_inside = this->is_strictly_inside3D(f.S.get_centroid(), shape);
                BOPAlgo_Splitter splitter;
                splitter.SetNonDestructive(true);
                splitter.AddArgument(sh);
                splitter.AddTool(f.S.get_shape());
                splitter.Perform();
                sh = splitter.Shape();
            }
        }
        for(auto& s:supports){
            //has_condition_inside = this->is_strictly_inside3D(s.S.get_centroid(), shape);
            if(s.cut_into_shape && s.S.get_area() > 0){
                BOPAlgo_Splitter splitter;
                splitter.SetNonDestructive(true);
                splitter.AddArgument(sh);
                splitter.AddTool(s.S.get_shape());
                splitter.Perform();
                sh = splitter.Shape();
            }
        }
    }

    // if(has_condition_inside){
    //     // Workaround so that this does not break current (faster) method of
    //     // distributing elements to the different geometry instances, at least
    //     // considering global meshing for linear analysis.
    //     //
    //     // Its current use is mostly for simple topopt/beam sizing, so it's a good
    //     // workaround for now.
    //     //
    //     // Otherwise, the simpler idea would be to test for the center of mass
    //     // of each geometry, which would be a problem if the geometry has a 
    //     // hole in its center.
    //     //
    //     // Seemed better to sacrifice this gimmick than something more useful
    //     // such as support for complex geometries.
    //     logger::log_assert(this->geometries.size() == 1, logger::ERROR, "applying boundary conditions inside a geometry is currently not supported when the number of geometries is greater than 1.");
    // }

    gmsh::initialize();
    gmsh::option::setNumber("Geometry.AutoCoherence", 2);
    gmsh::option::setNumber("Mesh.MeshSizeMin", 0);
    gmsh::option::setNumber("Mesh.MeshSizeMax", this->size);

    gmsh::option::setNumber("Mesh.Algorithm", this->algorithm2D);
    gmsh::option::setNumber("Mesh.Algorithm3D", this->algorithm3D);

    const size_t order = this->elem_info->get_element_order();
    gmsh::option::setNumber("Mesh.ElementOrder", order);
    //if(order > 1){
    //    gmsh::option::setNumber("Mesh.HighOrderOptimize", 2);
    //}
    gmsh::option::setNumber("Mesh.Optimize", 1);
    //gmsh::option::setNumber("Mesh.OptimizeNetgen", 1);
    gmsh::option::setNumber("Mesh.OptimizeThreshold", 0.3);
    gmsh::option::setNumber("Geometry.Tolerance", 1e-4);
    gmsh::option::setNumber("Mesh.Smoothing", 200);
    gmsh::option::setNumber("Mesh.AnisoMax", 1e1);
    gmsh::option::setNumber("Mesh.AllowSwapAngle", 90);
    gmsh::option::setNumber("Mesh.RandomFactor", 1e-7);
    gmsh::option::setNumber("Mesh.RefineSteps", 100);

    gmsh::option::setNumber("Mesh.AngleToleranceFacetOverlap", 1);
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", this->size_from_curvature);
    gmsh::option::setNumber("Mesh.ToleranceInitialDelaunay", 1e-2);
    gmsh::option::setNumber("Mesh.LcIntegrationPrecision", 1e-7);

    // Quad/hex recombination
    if(this->elem_info->get_shape_type() == Element::Shape::QUAD){
        gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 2);
        gmsh::option::setNumber("Mesh.RecombineAll", 1);
        gmsh::option::setNumber("Mesh.Recombine3DAll", 1);
        if(problem_type == utils::PROBLEM_TYPE_2D){
            gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 1);
        } else if(problem_type == utils::PROBLEM_TYPE_3D){
            gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 2);
        }
    }

    if(this->geometries.size() > 1 && problem_type == utils::PROBLEM_TYPE_3D){
        BOPAlgo_Builder sec;
        sec.SetNonDestructive(true);
        sec.SetRunParallel(true);
        sec.SetUseOBB(false);
        for(size_t i = 0; i < this->geometries.size(); ++i){
            sec.AddArgument(this->geometries[i]->shape);
        }
        sec.Perform();
        sec.DumpErrors(std::cout);
        sh = sec.Shape();

    }

    std::vector<size_t> geom_elem_mapping, elem_node_tags, bound_elem_node_tags;
    std::unordered_map<size_t, size_t> duplicate_map;
    auto id_map = this->gmsh_meshing(has_condition_inside, sh, geom_elem_mapping, elem_node_tags, bound_elem_node_tags, this->elem_info, duplicate_map);

    bool deduplicate = false;
    if(geometries.size() > 1){
        deduplicate = true;
    }
    this->generate_elements(geom_elem_mapping, 
                            elem_node_tags, 
                            bound_elem_node_tags,
                            id_map,
                            duplicate_map,
                            deduplicate,
                            has_condition_inside);
}

std::unordered_map<size_t, MeshNode*> Gmsh::gmsh_meshing(bool has_condition_inside, TopoDS_Shape sh, std::vector<size_t>& geom_elem_mapping, std::vector<size_t>& elem_node_tags, std::vector<size_t>& bound_elem_node_tags, const MeshElementFactory* const elem_type, std::unordered_map<size_t, size_t>& duplicate_map){
    (void) has_condition_inside;

    const size_t type = elem_type->get_gmsh_element_type();
    const size_t bound_type = elem_type->get_boundary_gmsh_element_type();
    const auto problem_type = this->elem_info->get_problem_type();
    const size_t nodes_per_elem = this->elem_info->get_nodes_per_element();

    size_t dim = 0;
    if(problem_type == utils::PROBLEM_TYPE_2D){
        dim = 2;
    } else if(problem_type == utils::PROBLEM_TYPE_3D){
        dim = 3;
    }

    std::vector<size_t> node_tags, boundary_node_tags;
    std::vector<double> node_coords;
    geom_elem_mapping.resize(geometries.size(),0);

    if(geometries.size() == 1){
        gmsh::model::add("base");
        gmsh::vectorpair vec;
        gmsh::model::occ::importShapesNativePointer(static_cast<const void*>(&sh), vec);
        gmsh::model::occ::synchronize();

        gmsh::model::mesh::generate(dim);

        //gmsh::write("mesh.msh");

        std::vector<int> elem_types;
        gmsh::model::mesh::getElementTypes(elem_types);
        // Check if meshing went well
        logger::log_assert(std::find(elem_types.begin(), elem_types.end(), type) != elem_types.end(), logger::ERROR,
                            "element type not found in mesh's list of element types (this shouldn't happen).");
        logger::log_assert(std::find(elem_types.begin(), elem_types.end(), bound_type) != elem_types.end(), logger::ERROR,
                            "element type of boundary elements not found in mesh's list of element types (this shouldn't happen).");

        // Get nodes and their tags
        std::vector<double> node_params;
        gmsh::model::mesh::getNodes(node_tags, node_coords, node_params, dim, -1, true);
        std::vector<double> bnode_coords, bnode_params;
        gmsh::model::mesh::getNodes(boundary_node_tags, bnode_coords, bnode_params, dim-1, -1, true);

        // Would need to be changed to support multiple element types
        std::vector<size_t> elem_tags;
        gmsh::model::mesh::getElementsByType(type, elem_tags, elem_node_tags, -1);
        geom_elem_mapping[0] = elem_tags.size();
        elem_tags.clear();
        gmsh::model::mesh::getElementsByType(bound_type, elem_tags, bound_elem_node_tags);

        gmsh::clear();
    } else {
        gmsh::model::add("base");
        gmsh::vectorpair vec;
        gmsh::model::occ::importShapesNativePointer(static_cast<const void*>(&sh), vec);
        gmsh::model::occ::removeAllDuplicates();
        gmsh::model::occ::synchronize();

        // Refine mesh close to contact boundary (optional)
        if(this->geometries.size() > 1 && this->bound_ref.defined()){
            //const size_t N = vec.size()-1;
            //const size_t num_pairs = N*(N+1)/2;
            std::vector<double> result = this->bound_ref.get_faces(vec);

            gmsh::model::mesh::field::add("Distance", 1);
            gmsh::model::mesh::field::setNumber(1, "Sampling", 100);
            gmsh::model::mesh::field::setNumbers(1, "SurfacesList", result);

            //this->contact_size = 0.2; // temp
            gmsh::model::mesh::field::add("Threshold", 2);
            gmsh::model::mesh::field::setNumber(2, "InField", 1);
            gmsh::model::mesh::field::setNumber(2, "SizeMin", this->bound_ref.min_size);
            gmsh::model::mesh::field::setNumber(2, "SizeMax", this->size);
            gmsh::model::mesh::field::setNumber(2, "DistMin", this->bound_ref.min_dist);
            gmsh::model::mesh::field::setNumber(2, "DistMax", this->bound_ref.max_dist);

            //gmsh::model::mesh::field::add("AutomaticMeshSizeField", 3);
            //gmsh::model::mesh::field::setNumber(3, "hMax", this->size);
            //gmsh::model::mesh::field::setNumber(3, "hMin", 0);
            //gmsh::model::mesh::field::setNumber(3, "nPointsPerCircle", this->size_from_curvature);

            //gmsh::model::mesh::field::add("Min", 100);
            //gmsh::model::mesh::field::setNumbers(100, "FieldsList", {2, 3});

            gmsh::model::mesh::field::setAsBackgroundMesh(2);
            // https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_15_0/tutorials/c++/t10.cpp#L48
        }

        gmsh::model::mesh::generate(dim);

        std::vector<int> elem_types;
        gmsh::model::mesh::getElementTypes(elem_types);
        //logger::quick_log(elem_types);
        // Check if meshing went well
        logger::log_assert(std::find(elem_types.begin(), elem_types.end(), type) != elem_types.end(), logger::ERROR,
                            "element type not found in mesh's list of element types (this shouldn't happen).");
        logger::log_assert(std::find(elem_types.begin(), elem_types.end(), bound_type) != elem_types.end(), logger::ERROR,
                            "element type of boundary elements not found in mesh's list of element types (this shouldn't happen).");

        std::vector<double> node_params;
        gmsh::model::mesh::getNodes(node_tags, node_coords, node_params, dim, -1, true);
        std::vector<double> bnode_coords, bnode_params;
        gmsh::model::mesh::getNodes(boundary_node_tags, bnode_coords, bnode_params, dim-1, -1, true);

        std::vector<size_t> elem_tags;
        // Boundary elements
        gmsh::model::mesh::getElementsByType(bound_type, elem_tags, bound_elem_node_tags);
        elem_tags.clear();
        // I'm assuming here that the geometry's id and the internal model id
        // in Gmsh/OCCT are the same.
        size_t end = gmsh::model::occ::getMaxTag(dim);
        logger::log_assert(end == geometries.size(), logger::ERROR, "gmsh found more volumes than there are geometries. Please ensure that the loaded geometries are not compounds.");
        for(size_t i = 1; i <= end; ++i){
            std::vector<size_t> elem_tags_tmp;
            std::vector<size_t> elem_node_tags_tmp;
            gmsh::model::mesh::getElementsByType(type, elem_tags_tmp, elem_node_tags_tmp, i);
            elem_node_tags.insert(elem_node_tags.end(), elem_node_tags_tmp.begin(), elem_node_tags_tmp.end());
            geom_elem_mapping[i-1] = elem_tags_tmp.size();
        }
        for(size_t j = 1; j < geom_elem_mapping.size(); ++j){
            geom_elem_mapping[j] += geom_elem_mapping[j-1];
        }

        elem_tags.clear();

        size_t node_tag_new = 1+*std::max_element(node_tags.begin(), node_tags.end());

        struct AddedNode{
            size_t g1, g2;
            size_t id;
            gp_Pnt p;
        };
        std::unordered_map<size_t, AddedNode> new_node_tag_map;
        std::vector<std::vector<size_t>> new_bound_elems;
        {
            const size_t N = vec.size() - 1;
            new_bound_elems.reserve(N*(N+1)/2);
        }

        for(size_t i = 0; i < vec.size(); ++i){
            for(size_t j = i+1; j < vec.size(); ++j){
                gmsh::vectorpair tags{vec[i], vec[j]};

                gmsh::vectorpair b1, bc;

                gmsh::model::getBoundary(tags, bc, true, false);
                gmsh::model::getBoundary({vec[i]}, b1, true, false);

                std::set<int> vc_set;
                for(const auto& bi:bc){
                    vc_set.insert(bi.second);
                }
                for(const auto& bi:b1){
                    const auto res = vc_set.find(bi.second);
                    if(res == vc_set.end()){
                        std::vector<size_t> node_tags_tmp;
                        std::vector<double> node_coords_tmp, node_params_tmp;
                        gmsh::model::mesh::getNodes(node_tags_tmp, node_coords_tmp, node_params_tmp, dim-1, bi.second, true);

                        std::vector<size_t> elem_tags_tmp, bound_elem_node_tags_tmp;
                        gmsh::model::mesh::getElementsByType(bound_type, elem_tags_tmp, bound_elem_node_tags_tmp, bi.second);

                        new_bound_elems.push_back(std::move(bound_elem_node_tags_tmp));

                        for(size_t k = 0; k < node_tags_tmp.size(); ++k){
                            size_t old_id = node_tags_tmp[k];
                            auto pos = new_node_tag_map.find(old_id);
                            if(pos == new_node_tag_map.end()){
                                new_node_tag_map[old_id] = AddedNode{
                                    i, j, node_tag_new,
                                    gp_Pnt(
                                        node_coords_tmp[3*k + 0],
                                        node_coords_tmp[3*k + 1],
                                        node_coords_tmp[3*k + 2]
                                    )};
                                duplicate_map[old_id] = node_tag_new;

                                ++node_tag_new;
                            }
                        }
                    }
                }
            }
        }
        const size_t new_nodes_num = new_node_tag_map.size();
        const size_t old_nodes_size = node_tags.size();
        const size_t old_bound_nodes_size = boundary_node_tags.size();
        node_tags.resize(node_tags.size() + new_nodes_num);
        node_coords.resize(node_coords.size() + 3*new_nodes_num);
        boundary_node_tags.resize(boundary_node_tags.size() + new_nodes_num);
        auto new_node_it = new_node_tag_map.begin();
        for(size_t i = 0; i < new_nodes_num; ++i){
            node_tags[old_nodes_size + i] = new_node_it->second.id;
            boundary_node_tags[old_bound_nodes_size + i] = new_node_it->second.id;
            for(size_t j = 0; j < 3; ++j){
                node_coords[3*(old_nodes_size + i)+j] = new_node_it->second.p.Coord(1+j);
            }
            ++new_node_it;
        }

        size_t new_bound_elem_size = 0;
        for(auto& bl:new_bound_elems){
            new_bound_elem_size += bl.size();
            for(auto& bn:bl){
                bn = new_node_tag_map[bn].id;
            }
        }
        const size_t old_bound_elem_size = bound_elem_node_tags.size();
        bound_elem_node_tags.reserve(old_bound_elem_size + new_bound_elem_size);
        auto b_it = bound_elem_node_tags.begin() + old_bound_elem_size;
        for(auto& bl:new_bound_elems){
            b_it = bound_elem_node_tags.insert(b_it, bl.begin(), bl.end());
            b_it += bl.size();
        }

        for(size_t i = 1; i < vec.size(); ++i){
            #pragma omp parallel for
            for(size_t j = nodes_per_elem*geom_elem_mapping[i-1]; j < nodes_per_elem*geom_elem_mapping[i]; ++j){
                auto pos = new_node_tag_map.find(elem_node_tags[j]);
                if(pos != new_node_tag_map.end() && pos->second.g1 != i){
                    elem_node_tags[j] = pos->second.id;
                }
            }
        }

        gmsh::clear();
    }

    gmsh::finalize();

    const size_t dof = elem_type->get_dof_per_node();

    const auto node_comp = [](const MeshNode& n1, const MeshNode& n2) -> bool{
        return n1.id < n2.id;
    };
    std::set<MeshNode, decltype(node_comp)> node_set(node_comp);
    for(size_t i = 0; i < node_tags.size(); ++i){
        gp_Pnt p(node_coords[i*3]/this->tmp_scale, node_coords[i*3+1]/this->tmp_scale, node_coords[i*3+2]/this->tmp_scale);

        node_set.emplace(p, node_tags[i], dof);
    }

    std::unordered_map<size_t, MeshNode*> id_map;
    id_map.reserve(node_set.size());

    this->node_list.clear();
    this->node_list.reserve(node_set.size());
    for(auto it = node_set.cbegin(); it != node_set.cend(); ++it){

        this->node_list.emplace_back(std::make_unique<MeshNode>(it->point, it->id, dof));
        id_map.emplace(it->id, this->node_list.back().get());
    }
    node_set.clear();

    std::set<size_t> boundary_node_set(boundary_node_tags.begin(), boundary_node_tags.end());

    this->boundary_node_list.clear();
    this->boundary_node_list.reserve(boundary_node_set.size());
    for(auto it = boundary_node_set.cbegin(); it != boundary_node_set.cend(); ++it){
        this->boundary_node_list.push_back(id_map.at(*it));
    }

    return id_map;
}

using namespace projspec;
const bool Gmsh::reg = Factory<Meshing>::add(
    [](const DataMap& data){
        return std::make_unique<Gmsh>(data);
    },
    ObjectRequirements{
        "gmsh",
        {
            DataEntry{.name = "element_size", .type = TYPE_DOUBLE, .required = true},
            DataEntry{.name = "tmp_scale", .type = TYPE_DOUBLE, .required = false},
            DataEntry{.name = "size_from_curvature", .type = TYPE_DOUBLE, .required = false},
            DataEntry{.name = "algorithm2D", .type = TYPE_INT, .required = false},
            DataEntry{.name = "algorithm3D", .type = TYPE_INT, .required = false},
            DataEntry{.name = "boundary_refinement", .type = TYPE_OBJECT, .required = false,
                .object_data = {
                    DataEntry{.name = "min_dist", .type = TYPE_DOUBLE, .required = true},
                    DataEntry{.name = "max_dist", .type = TYPE_DOUBLE, .required = true},
                    DataEntry{.name = "min_size", .type = TYPE_DOUBLE, .required = true},
                    DataEntry{.name = "region", .type = TYPE_SHAPE_OP, .required = true},
                }
            }
        }
    }
);

}
