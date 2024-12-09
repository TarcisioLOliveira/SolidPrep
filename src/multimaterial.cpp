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

#include <string>
#include <array>
#include <cmath>
#include "math/matrix.hpp"
#include "utils/D_operations.hpp"
#include "logger.hpp"
#include "multimaterial.hpp"
#include "utils.hpp"

MultiMaterial::MultiMaterial(std::vector<Material*> materials, utils::ProblemType type, bool has_void):
    materials(std::move(materials)), problem_type(type), has_void(has_void){

    for(const auto& M:this->materials){    
       if(!M->is_homogeneous()){
           this->is_homogeneous = false;
           break;
       } 
    } 
}

void MultiMaterial::get_D(std::vector<double>::const_iterator rho, const double mix, const MeshElement* e, const gp_Pnt& p, math::Matrix& D) const{
    if(has_void){
        ++rho;
    }

    if(this->problem_type == utils::PROBLEM_TYPE_2D){
        if(this->materials.size() > 1){
            this->get_D_internal(rho, mix, e, p, D);
        } else {
            const auto Dtmp = this->materials[0]->stiffness_2D(e, p);
            std::copy(Dtmp.data(), Dtmp.data() + 3*3, D.data());
        }
    } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
        if(this->materials.size() > 1){
            this->get_D_internal(rho, mix, e, p, D);
        } else {
            const auto Dtmp = this->materials[0]->stiffness_3D(e, p);
            std::copy(Dtmp.data(), Dtmp.data() + 6*6, D.data());
        }
    }
}

void MultiMaterial::get_gradD(std::vector<double>::const_iterator rho, const double mix, const MeshElement* e, const gp_Pnt& p, std::vector<math::Matrix>& gradD) const{
    if(has_void){
        ++rho;
    }

    if(this->problem_type == utils::PROBLEM_TYPE_2D){
        if(this->materials.size() > 1){
            this->get_gradD_internal(rho, mix, e, p, gradD);
        }
        if(has_void){
            const auto Dtmp = this->materials[0]->stiffness_2D(e, p);
            std::copy(Dtmp.data(), Dtmp.data() + 3*3, gradD[0].data());
        }
    } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
        if(this->materials.size() > 1){
            this->get_gradD_internal(rho, mix, e, p, gradD);
        }
        if(has_void){
            const auto Dtmp = this->materials[0]->stiffness_3D(e, p);
            std::copy(Dtmp.data(), Dtmp.data() + 6*6, gradD[0].data());
        }
    }
}

void MultiMaterial::get_D_internal(std::vector<double>::const_iterator& rho, const double mix, const MeshElement* e, const gp_Pnt& p, math::Matrix& D) const{
    std::vector<double> v(this->materials.size(), 0);
    double total = 0;
    for(size_t i = 0; i < v.size()-1; ++i){
        total += *(rho+i);
        v[i] = *(rho+i);
    }
    v[v.size()-1] = 1 - total;

    const double psi_v = std::abs(mix);
    const double psi_r = 1 - psi_v;

    size_t N = 0;
    math::Matrix vD;
    math::Matrix rD;
    if(this->problem_type == utils::PROBLEM_TYPE_2D){
        N = 3;
    } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
        N = 6;
    }
    vD = math::Matrix(N, N);
    rD = math::Matrix(N, N);

    if(mix >= 0){
        for(size_t j = 0; j < v.size(); ++j){
            if(this->problem_type == utils::PROBLEM_TYPE_2D){
                const auto MD = this->materials[j]->stiffness_2D(e, p);
                const auto MS = this->materials[j]->stiffness_inverse_2D(e, p);
                vD += v[j]*MD;
                rD += v[j]*MS;
            } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                const auto MD = this->materials[j]->stiffness_3D(e, p);
                const auto MS = this->materials[j]->stiffness_inverse_3D(e, p);
                vD += v[j]*MD;
                rD += v[j]*MS;
            }
        }
        rD.invert_cholesky();
        rD *= psi_r;
        vD *= psi_v;
        D = rD + vD;
    } else {
        for(size_t j = 0; j < v.size(); ++j){
            if(this->problem_type == utils::PROBLEM_TYPE_2D){
                const auto MD = this->materials[j]->stiffness_2D(e, p);
                const auto MS = MD*MD;
                vD += v[j]*MD;
                rD += v[j]*MS;
            } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                const auto MD = this->materials[j]->stiffness_3D(e, p);
                const auto MS = MD*MD;
                vD += v[j]*MD;
                rD += v[j]*MS;
            }
        }
        math::Eigen e(rD);
        rD = e.square_root();
        rD *= psi_r;
        vD *= psi_v;
        D = rD + vD;
    }
}

void MultiMaterial::get_gradD_internal(std::vector<double>::const_iterator& rho, const double mix, const MeshElement* e, const gp_Pnt& p, std::vector<math::Matrix>& gradD) const{
    std::vector<double> v(this->materials.size(), 0);
    double total = 0;
    for(size_t i = 0; i < v.size()-1; ++i){
        total += *(rho+i);
        v[i] = *(rho+i);
    }
    v[v.size()-1] = 1 - total;

    const double psi_v = std::abs(mix);
    const double psi_r = 1 - psi_v;
    const size_t offset = (has_void) ? 1 : 0;

    size_t N = 0;
    math::Matrix vD;
    math::Matrix rD;
    if(this->problem_type == utils::PROBLEM_TYPE_2D){
        N = 3;
    } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
        N = 6;
    }
    vD = math::Matrix(N, N);
    rD = math::Matrix(N, N);

    if(mix >= 0){
        math::Matrix rD_last;
        math::Matrix rD_mult;
        // Reuss matrix for last material
        math::Matrix MDl;
        if(this->problem_type == utils::PROBLEM_TYPE_2D){
            MDl = this->materials.back()->stiffness_2D(e, p);
            rD_last = this->materials.back()->stiffness_inverse_2D(e, p);
        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
            MDl = this->materials.back()->stiffness_3D(e, p);
            rD_last = this->materials.back()->stiffness_inverse_3D(e, p);
        }
        // Full Reuss matrix, necessary due to derivative chain rule
        std::vector<math::Matrix> S_cache(v.size());
        {
            for(size_t j = 0; j < v.size(); ++j){
                if(this->problem_type == utils::PROBLEM_TYPE_2D){
                    const auto MS = this->materials[j]->stiffness_inverse_2D(e, p);
                    S_cache[j] = MS;
                    rD += v[j]*MS;
                } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                    const auto MS = this->materials[j]->stiffness_inverse_3D(e, p);
                    S_cache[j] = MS;
                    rD += v[j]*MS;
                }
            }
            rD_mult = rD.get_inverted_cholesky();
        }
        // Generate gradients
        math::Matrix D_tmp(N, N);
        for(size_t j = 0; j < v.size()-1; ++j){
            if(this->problem_type == utils::PROBLEM_TYPE_2D){
                const auto MD = this->materials[j]->stiffness_2D(e, p);
                const auto MS = S_cache[j];
                D_tmp = MS - rD_last;
                gradD[j+offset] = psi_v*(MD - MDl) - psi_r*rD_mult*(MS - rD_last)*rD_mult;
            } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                const auto MD = this->materials[j]->stiffness_3D(e, p);
                const auto MS = S_cache[j];
                D_tmp = MS - rD_last;
                gradD[j+offset] = psi_v*(MD - MDl) - psi_r*rD_mult*(MS - rD_last)*rD_mult;
            }
        }
    } else {
        math::Matrix rD_last;
        math::Matrix rD_mult;
        // Reuss matrix for last material
        math::Matrix MDl;
        if(this->problem_type == utils::PROBLEM_TYPE_2D){
            MDl = this->materials.back()->stiffness_2D(e, p);
            rD_last = MDl*MDl;
        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
            MDl = this->materials.back()->stiffness_3D(e, p);
            rD_last = MDl*MDl;
        }
        // Full Reuss matrix, necessary due to derivative chain rule
        std::vector<math::Matrix> D_cache(v.size());
        {
            for(size_t j = 0; j < v.size(); ++j){
                if(this->problem_type == utils::PROBLEM_TYPE_2D){
                    D_cache[j] = this->materials[j]->stiffness_2D(e, p);
                } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                    D_cache[j] = this->materials[j]->stiffness_3D(e, p);
                }
                const auto MS = D_cache[j]*D_cache[j];
                rD += v[j]*MS;
            }
            math::Eigen e(rD);
            rD_mult = e.square_root();
            rD_mult.invert_cholesky();
            rD_mult *= 0.5;

            rD_last = rD_mult*rD_last;
        }
        // Generate gradients
        for(size_t j = 0; j < v.size()-1; ++j){
            const auto& MD = D_cache[j];
            const auto MS = MD*MD;
            // The last part is actually rD_mult*(MS - rD_last), but
            // rD_mult*rD_last is precalculated, as it is constant
            gradD[j+offset] = psi_v*(MD - MDl) + psi_r*(rD_mult*MS - rD_last);
        }
    }
}

double MultiMaterial::get_density(std::vector<double>::const_iterator rho, const MeshElement* e, const gp_Pnt& p) const{
    double void_rho = 1.0;
    if(has_void){
        void_rho = *rho;
        ++rho;
    }

    std::vector<double> v(this->materials.size(), 0);
    double total = 0;
    for(size_t i = 0; i < v.size()-1; ++i){
        total += *(rho+i);
        v[i] = *(rho+i);
    }
    v[v.size()-1] = 1 - total;

    rho += v.size()-1;

    double d = 0;
    for(size_t i = 0; i < v.size(); ++i){
        d += v[i]*this->materials[i]->get_density(e, p);
    }

    d *= void_rho;

    return d;
}

double MultiMaterial::get_density_deriv(std::vector<double>::const_iterator rho, const MeshElement* e, const gp_Pnt& p, std::vector<double>::iterator& grad) const{
    double void_rho = 1.0;
    if(has_void){
        void_rho = *rho;
        ++rho;
    }

    std::vector<double> v(this->materials.size(), 0);
    double total = 0;
    for(size_t i = 0; i < v.size()-1; ++i){
        total += *(rho+i);
        v[i] = *(rho+i);
    }
    v[v.size()-1] = 1 - total;

    rho += v.size()-1;

    double d = 0;
    for(size_t i = 0; i < v.size(); ++i){
        d += v[i]*this->materials[i]->get_density(e, p);
    }

    if(has_void){
        *grad = d;
        ++grad;
    }
    for(size_t i = 0; i < v.size()-1; ++i){
        *grad = void_rho*(this->materials[i]->get_density(e, p) - this->materials.back()->get_density(e, p));
        ++grad;
    }

    return void_rho*d;
}
