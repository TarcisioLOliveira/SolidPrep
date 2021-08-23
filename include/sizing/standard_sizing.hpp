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

#ifndef STANDARD_SIZING_HPP
#define STANDARD_SIZING_HPP

#include "sizing.hpp"
#include "finite_element.hpp"
#include "beam_meshing.hpp"
#include "meshing/standard_beam_mesher.hpp"

namespace sizing{

class StandardSizing : public Sizing{
    public:
    struct ExpansionNode{
        gp_Pnt center;
        gp_Dir direction;
        double diameter;
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

    StandardSizing(ProjectData* data, FiniteElement* solver);

    virtual TopoDS_Shape run() override;

    private:
    FiniteElement* solver;

    TopoDS_Shape build_initial_topology();
    // Doesn't work very well. Can distort some beams.
    TopoDS_Shape simplify_shape(TopoDS_Shape shape) const;
    bool is_inside_2D(const gp_Pnt& p, const TopoDS_Shape& shape) const;
    bool is_inside_3D(const gp_Pnt& p, const TopoDS_Shape& shape) const;

    TopoDS_Shape boundary_expansion_approach();
    bool is_valid_boundary_point(MeshNode* n) const;
    TopoDS_Shape expansion_2D(const meshing::StandardBeamMesher& mesh, const std::vector<double>& u, const TopoDS_Shape& beams);
    ExpansionNode get_expansion_node_2D(const gp_Dir& line_dir, gp_Pnt center, double distance, double Fx, double Fy, double Mz, const std::vector<TopoDS_Edge>& edges_init = std::vector<TopoDS_Edge>()) const;
   
    /**
     * Related to elemental approach. 
     * Fully functional, but does not give the expected results. 
     */
    std::vector<gp_Pnt> end_points;
    TopoDS_Shape experimental_elemental_approach();
    std::vector<double> calculate_change(BeamMeshing* mesh, const std::vector<long>& ids, std::vector<double> h, const TopoDS_Shape& beams) const;
};

}

#endif
