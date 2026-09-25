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
#include "logger.hpp"
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
    virtual std::vector<math::Matrix> get_dsh_full(const MeshElement* e, const gp_Pnt& p, const size_t n, const size_t dof) const override;

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
    void generate_vectors();
    inline math::Matrix get_eigen(const MeshElement* e, const gp_Pnt& p) const{
        math::Matrix tensor;
        if(!this->use_smooth){
            tensor = this->get_tensor_base(e, p);
        } else {
            tensor = this->get_tensor_interp(e, p);
        }
        math::Eigen eigen(tensor, 1e-14);
        std::vector<size_t> pos(DIM);
        std::iota(pos.begin(), pos.end(), 0);
        std::vector<double> abs_sort(DIM);
        for(size_t i = 0; i < DIM; ++i){
            abs_sort[i] = std::abs(eigen.eigenvalues[i]);
        }
        std::sort(pos.begin(), pos.end(), [&](size_t i, size_t j){ return abs_sort[i] < abs_sort[j]; });

        math::Matrix dirs(DIM, DIM);
        for(size_t i = 0; i < DIM; ++i){
            for(size_t j = 0; j < DIM; ++j){
                dirs(i, j) = eigen.eigenvectors(i, pos[j]);
            }
        }
        return dirs;
    }
    inline math::Matrix get_tensor_base(const MeshElement* e, const gp_Pnt& p) const{
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
        return tensor;
    }
    inline math::Matrix get_tensor_interp(const MeshElement* e, const gp_Pnt& p) const{
        math::Matrix tensor(3, 3, 0);
        const auto N = e->get_Ni_1dof(p);

        for(size_t i = 0; i < N.get_N(); ++i){
            const size_t id = this->id_pos_map.at(e->nodes[i]->id);
            size_t count = 0;
            for(size_t j = 0; j < 3; ++j){
                for(size_t k = 0; k <= j; ++k){
                    if(j == k){
                        tensor(j,k) += N[i]*this->nodal_strains[6*id + count];
                    } else {
                        tensor(j,k) += N[i]*this->nodal_strains[6*id + count];
                        tensor(k,j) += N[i]*this->nodal_strains[6*id + count];
                    }
                    ++count;
                }
            }
        }
        return tensor;
    }
    math::Matrix get_tensor_dsh(const MeshElement* e, const gp_Pnt& p, const size_t n, const size_t dof) const;

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
    bool use_smooth = false;
    std::vector<double> u;
    size_t max_it;
    const bool use_strain;
    std::unordered_multimap<size_t, MeshElement*> inverse_mesh;
    const double ALPHA;

    std::map<size_t, size_t> id_pos_map;
    std::vector<UniqueNodeForView> node_view_list;

    Visualization* viz;

    std::vector<double> nodal_weights;
    std::vector<double> nodal_strains;

    ViewHandler* L_view = nullptr;
    ViewHandler* R_view = nullptr;
    ViewHandler* T_view = nullptr;
};

}

#endif
