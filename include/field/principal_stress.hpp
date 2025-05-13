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

#ifndef PRINCIPAL_STRESS_HPP
#define PRINCIPAL_STRESS_HPP

#include <map>
#include <unordered_map>
#include <vector>
#include "field.hpp"
#include "math/matrix.hpp"
#include "project_specification/data_map.hpp"
#include "utils/delayed_pointer.hpp"

namespace field{

class PrincipalStress : public CoordinateField {
    public:
    virtual ~PrincipalStress() = default;
    PrincipalStress(const projspec::DataMap& data);

    virtual void generate() override;
    virtual void initialize_views(Visualization* viz) override;
    virtual void display_views() const override;

    virtual std::array<gp_Dir, 3> get_array(const MeshElement* e, const gp_Pnt& p) const override;
    virtual math::Matrix get_matrix(const MeshElement* e, const gp_Pnt& p) const override;

    inline virtual Class get_class() const override{
        return Class::PRINCIPAL_STRESS;
    }
    inline virtual SubType get_sub_type() const override{
        return SubType::DOMAIN;
    }
    inline virtual bool is_fea_dependent() const override{ return true; }

    private:
    static const bool reg;
    class UniqueNodeForView{
        public:
        const Node* node;
        const MeshElement* element;

        bool operator<(const UniqueNodeForView& n) const{
            return this->node->id < n.node->id;
        }
    };
    ProjectData* proj_data;
    utils::DelayedPointerView<Material> initial_material;
    const MeshElementFactory* elem_info;
    std::vector<Geometry*> geoms;
    std::vector<math::Matrix> base_D;
    std::unordered_map<const MeshElement*, math::Matrix*> elem_mapping;
    double thickness;
    bool show;
    size_t DIM;
    size_t NODES_PER_ELEM;
    bool first_run = true;
    std::vector<double> u;
    const size_t max_it;

    std::map<size_t, size_t> id_pos_map;
    std::vector<UniqueNodeForView> node_view_list;

    Visualization* viz;

    ViewHandler* L_view = nullptr;
    ViewHandler* R_view = nullptr;
    ViewHandler* T_view = nullptr;
};

}

#endif
