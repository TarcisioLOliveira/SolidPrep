/*
 *   Copyright (C) 2025 Tarcísio Ladeia de Oliveira.
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

#ifndef TRABECULAR_CUTOFF_HPP
#define TRABECULAR_CUTOFF_HPP

#include <map>
#include <unordered_map>
#include <vector>
#include "field.hpp"
#include "project_specification/data_map.hpp"

namespace field{

class TrabecularCutoff : public ScalarField {
    public:
    virtual ~TrabecularCutoff() = default;
    TrabecularCutoff(const projspec::DataMap& data);

    virtual void generate() override;
    virtual void initialize_views(Visualization* viz) override;
    virtual void display_views() const override;
    virtual double get(const MeshElement* e, const gp_Pnt& p) const override{
        (void) p;
        return this->trabecular_map.at(e->id);
    }

    inline virtual SubType get_sub_type() const override{
        return SubType::DOMAIN;
    }

    private:
    static const bool reg;

    std::vector<Geometry*> geoms;

    ScalarField* density;
    std::unordered_map<size_t, double> trabecular_map; // 1 == trabecular; 0 == cortical

    bool show;
    const double cutoff;

    TopoDS_Shape shell;
    const double bone_thickness;


    ViewHandler* distr = nullptr;
};

}

#endif
