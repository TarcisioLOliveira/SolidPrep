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

#ifndef BEAM_SIZING_HPP
#define BEAM_SIZING_HPP

#include "sizing.hpp"
#include <vector>
#include "element_factory.hpp"

namespace sizing{

class BeamSizing : public Sizing{
    public:
    BeamSizing(ProjectData* data, BeamElementFactory::BeamElementType t);

    virtual TopoDS_Shape run() override;

    BeamElementFactory::BeamElementType type;
};

}

#endif
