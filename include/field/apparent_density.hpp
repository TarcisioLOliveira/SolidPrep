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

#ifndef APPARENT_DENSITY_HPP
#define APPARENT_DENSITY_HPP

#include <map>
#include <vector>
#include "field.hpp"
#include "math/matrix.hpp"
#include "project_specification/data_map.hpp"
#include "teem/nrrd.h"

namespace field{

class NRRDReader {
    public:
    NRRDReader();
    ~NRRDReader() = default;

    void load(const std::string& file_path);
    double get_averaged(const math::Vector& point) const;

    private:
    std::unique_ptr<Nrrd, std::function<void(Nrrd*)>> file;

    double get_data(const std::vector<size_t>& p) const;

    math::Vector origin;
    math::Matrix R;
    math::Matrix Rinv;
    double DIM;
    std::vector<size_t> len;

    double (*lup)(const void *, size_t I);
};

class ApparentDensity : public ScalarField {
    public:
    virtual ~ApparentDensity() = default;
    ApparentDensity(const projspec::DataMap& data);

    virtual void generate() override;
    virtual void initialize_views(Visualization* viz) override;
    virtual void display_views() const override;
    virtual double get(const MeshElement* e, const gp_Pnt& p) const override;

    inline virtual SubType get_sub_type() const override{
        return SubType::DOMAIN;
    }

    private:
    static const bool reg;

    const MeshElementFactory* elem_info;
    std::vector<Geometry*> geoms;
    double thickness;
    NRRDReader nrrd;
    math::Matrix R;
    math::Vector t;
    double radius;
    std::array<double, 2> rho_func;

    std::map<size_t, size_t> id_pos_map;
    std::vector<double> nodal_densities;

    Visualization* viz;

    bool show;
    size_t DIM;
    size_t NODES_PER_ELEM;

    ViewHandler* density = nullptr;
};

}

#endif
