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

namespace sizing{

class StandardSizing : public Sizing{
    public:
    StandardSizing(ProjectData* data, FiniteElement* solver);

    virtual TopoDS_Shape run() override;

    private:
    FiniteElement* solver;
    std::vector<gp_Pnt> end_points;

    TopoDS_Shape build_initial_topology();
    std::vector<double> calculate_change(BeamMeshing* mesh, const std::vector<long>& ids, std::vector<double> h, const TopoDS_Shape& beams) const;
    TopoDS_Shape simplify_shape(TopoDS_Shape shape) const;
};

}

#endif
