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
    virtual math::Matrix get_dsh(const MeshElement* e, const gp_Pnt& p, const size_t n, const size_t dof) const override;

    inline void set_max_it(size_t m){
        this->max_it = m;
    }

    inline virtual Class get_class() const override{
        return Class::PRINCIPAL_STRESS;
    }
    inline virtual SubType get_sub_type() const override{
        return SubType::DOMAIN;
    }
    inline virtual bool is_fea_dependent() const override{ return true; }

    protected:
    virtual void regenerate_internal() override{
        this->max_it = 1;
        this->generate();
    }

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
    inline math::Matrix get_eigen(const MeshElement* e, const gp_Pnt& p) const{
        math::Matrix tensor;
        if(this->use_strain){
            tensor = e->get_strain_tensor(p, this->u);
        } else {
            if(this->DIM == 3){
                const auto D = this->initial_material->stiffness_3D(e, p);
                tensor = e->get_stress_tensor(D, p, this->u);
            } else {
                const auto D = this->initial_material->stiffness_2D(e, p);
                tensor = e->get_stress_tensor(D, p, this->u);
            }
        }
        math::Eigen eigen(tensor, 1e-14);
        std::vector<double> abs_sort(DIM);
        for(size_t i = 0; i < DIM; ++i){
            abs_sort[i] = std::abs(eigen.eigenvalues[i]);
        }
        std::sort(abs_sort.begin(), abs_sort.end(), std::less<double>());
        std::vector<size_t> pos(DIM);
        for(size_t i = 0; i < DIM; ++i){
            for(size_t j = 0; j < DIM; ++j){
                if(abs_sort[i] == std::abs(eigen.eigenvalues[j])){
                    pos[i] = j;
                    break;
                }
            }
        }
        math::Matrix dirs(DIM, DIM);
        for(size_t i = 0; i < DIM; ++i){
            for(size_t j = 0; j < DIM; ++j){
                dirs(i, j) = eigen.eigenvectors(i, pos[j]);
            }
        }
        return dirs;
    }
    ProjectData* proj_data;
    utils::DelayedPointerView<Material> initial_material;
    const MeshElementFactory* elem_info;
    std::vector<Geometry*> geoms;
    double thickness;
    bool show;
    size_t DIM;
    size_t NODES_PER_ELEM;
    size_t DOF;
    bool first_run = true;
    std::vector<double> u;
    size_t max_it;
    const bool use_strain;

    std::map<size_t, size_t> id_pos_map;
    std::vector<UniqueNodeForView> node_view_list;

    Visualization* viz;

    ViewHandler* L_view = nullptr;
    ViewHandler* R_view = nullptr;
    ViewHandler* T_view = nullptr;
};

}

#endif
