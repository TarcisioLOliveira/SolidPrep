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

#ifndef ORTHOTROPIC_FLOW_HPP
#define ORTHOTROPIC_FLOW_HPP

#include <map>
#include <vector>
#include "field.hpp"

namespace field{

class OrthotropicFlow : public CoordinateField {
    public:
    OrthotropicFlow(const MeshElementFactory* elem_info, std::vector<Geometry*> geoms, std::vector<CrossSection> entries, std::vector<double> coeffs, double alpha, double thickness, bool show);

    virtual void generate() override;
    virtual void initialize_views(Visualization* viz) override;
    virtual void display_views() const override;

    virtual std::array<gp_Dir, 3> get_array(const MeshElement* e, const gp_Pnt& p) const override;
    virtual Eigen::Matrix<double, 3, 3> get_matrix(const MeshElement* e, const gp_Pnt& p) const override;

    inline virtual Class get_class() const override{
        return Class::ORTHOTROPIC_FLOW;
    }
    inline virtual SubType get_sub_type() const override{
        return SubType::DOMAIN;
    }

    private:
    class UniqueNodeForView{
        public:
        const Node* node;
        const MeshElement* element;

        bool operator<(const UniqueNodeForView& n) const{
            return this->node->id < n.node->id;
        }
    };
    const MeshElementFactory* elem_info;
    std::vector<Geometry*> geoms;
    std::vector<CrossSection> entries;
    std::vector<double> coeffs;
    std::vector<double> elem_mult_long;
    std::vector<double> elem_mult_rad;
    double ALPHA = 1;
    double thickness;
    bool show;
    size_t DIM;
    size_t NODES_PER_ELEM;

    std::map<size_t, size_t> id_pos_map;
    std::vector<UniqueNodeForView> node_view_list;

    std::vector<double> longitudinal;
    std::vector<double> radial;

    Visualization* viz;

    ViewHandler* L_view = nullptr;
    ViewHandler* R_view = nullptr;
    ViewHandler* T_view = nullptr;
};

}

#endif
