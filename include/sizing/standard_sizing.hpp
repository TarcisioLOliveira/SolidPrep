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

#ifndef STANDARD_SIZING_HPP
#define STANDARD_SIZING_HPP

#include "project_specification/data_map.hpp"
#include "sizing.hpp"
#include "finite_element.hpp"
#include "meshing/standard_beam_mesher.hpp"
#include <vector>

namespace sizing{

class StandardSizing : public Sizing{
    public:
    struct ExpansionNode{
        gp_Pnt center;
        gp_Dir direction;
        double diameter;
        std::vector<size_t> used_ef = std::vector<size_t>();
    };

    struct ExternalForce{
        gp_Vec forces = gp_Vec(0,0,0);
        gp_Vec moments = gp_Vec(0,0,0);
        gp_Pnt position = gp_Pnt(0,0,0);
        gp_Dir line_dir = gp_Dir(0,0,1);
        bool Mx = false;
        bool My = false;
        bool Mz = false;
        double diameter = 0;
    };

    struct IntersectionNode{
        gp_Pnt point;
        gp_Vec force;
    };

    StandardSizing(const projspec::DataMap& data);

    virtual TopoDS_Shape run() override;

    private:
    static const bool reg;
    FiniteElement* solver;
    double element_size;
    double multiplier;
    std::vector<gp_Pnt> end_points;

    TopoDS_Shape build_initial_topology();
    // Not recommended to use it with beams. Other shapes tend work.
    // Some kind of simplification is necessary for the resulting resized beam
    // though, as the small edges can make Gmsh generate way too many nodes.
    TopoDS_Shape simplify_shape(TopoDS_Shape shape) const;
    bool is_inside_2D(const gp_Pnt& p, const TopoDS_Shape& shape) const;
    bool is_inside_3D(const gp_Pnt& p, const TopoDS_Shape& shape) const;

    std::vector<TopoDS_Shape> separate_beams;
    TopoDS_Shape boundary_expansion_approach();
    bool is_valid_boundary_point(MeshNode* n) const;
    TopoDS_Shape expansion_2D(const meshing::StandardBeamMesher& mesh, const std::vector<double>& u, const TopoDS_Shape& beams);
    ExpansionNode get_expansion_node_2D(const gp_Dir& line_dir, gp_Pnt center, double distance, double Fx, double Fy, double Mz, const std::vector<TopoDS_Edge>& edges_init = std::vector<TopoDS_Edge>()) const;
    void calculate_reaction_moments(size_t Mn, std::vector<ExternalForce>& external_forces) const;

    /**
     * Currently unused.
     * Previously used to test for a way of ordering points along the beams' central
     * lines having only the points and their normals.
     * Does not work correctly, may be removed later.
     */
    TopoDS_Shape bspline_simple2D(const std::vector<ExpansionNode>& exp_info, TopoDS_Shape base) const;
    bool insert_expansion_node(std::vector<ExpansionNode>& exp_info, ExpansionNode node) const;
};

}

#endif
