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

void MultiMaterial::get_D(std::vector<double>::const_iterator rho, const double mix, const MeshElement* e, const gp_Pnt& p, std::vector<double>& D) const{
    if(has_void){
        ++rho;
    }

    if(this->problem_type == utils::PROBLEM_TYPE_2D){
        if(this->materials.size() > 1){
            this->get_D_internal(rho, mix, e, p, D);
        } else {
            const auto Dtmp = this->materials[0]->stiffness_2D(e, p);
            std::copy(Dtmp.begin(), Dtmp.end(), D.begin());
        }
    } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
        if(this->materials.size() > 1){
            this->get_D_internal(rho, mix, e, p, D);
        } else {
            const auto Dtmp = this->materials[0]->stiffness_3D(e, p);
            std::copy(Dtmp.begin(), Dtmp.end(), D.begin());
        }
    }
}

void MultiMaterial::get_gradD(std::vector<double>::const_iterator rho, const double mix, const MeshElement* e, const gp_Pnt& p, std::vector<std::vector<double>>& gradD) const{
    if(has_void){
        ++rho;
    }

    if(this->problem_type == utils::PROBLEM_TYPE_2D){
        if(this->materials.size() > 1){
            this->get_gradD_internal(rho, mix, e, p, gradD);
        }
        if(has_void){
            const auto Dtmp = this->materials[0]->stiffness_2D(e, p);
            std::copy(Dtmp.begin(), Dtmp.end(), gradD[0].begin());
        }
    } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
        if(this->materials.size() > 1){
            this->get_gradD_internal(rho, mix, e, p, gradD);
        }
        if(has_void){
            const auto Dtmp = this->materials[0]->stiffness_3D(e, p);
            std::copy(Dtmp.begin(), Dtmp.end(), gradD[0].begin());
        }
    }
}

void MultiMaterial::get_D_internal(std::vector<double>::const_iterator& rho, const double mix, const MeshElement* e, const gp_Pnt& p, std::vector<double>& D) const{
    std::vector<double> v(this->materials.size(), 0);
    double total = 0;
    for(size_t i = 0; i < v.size()-1; ++i){
        total += *(rho+i);
        v[i] = *(rho+i);
    }
    v[v.size()-1] = 1 - total;

    const double psi_v = std::abs(mix);
    const double psi_r = 1 - psi_v;

    size_t D_size = 0;
    std::vector<double> vD;
    std::vector<double> rD;
    if(this->problem_type == utils::PROBLEM_TYPE_2D){
        D_size = 9;
        vD.resize(D_size,0);
        rD.resize(D_size,0);
    } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
        D_size = 36;
        vD.resize(D_size,0);
        rD.resize(D_size,0);
    }

    if(mix >= 0){
        for(size_t j = 0; j < v.size(); ++j){
            if(this->problem_type == utils::PROBLEM_TYPE_2D){
                const auto MD = this->materials[j]->stiffness_2D(e, p);
                const auto MS = this->materials[j]->stiffness_inverse_2D(e, p);
                for(size_t i = 0; i < 9; ++i){
                    vD[i] += v[j]*MD[i];
                    rD[i] += v[j]*MS[i];
                }
            } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                const auto MD = this->materials[j]->stiffness_3D(e, p);
                const auto MS = this->materials[j]->stiffness_inverse_3D(e, p);
                for(size_t i = 0; i < 36; ++i){
                    vD[i] += v[j]*MD[i];
                    rD[i] += v[j]*MS[i];
                }
            }
        }
        if(this->problem_type == utils::PROBLEM_TYPE_2D){
            rD = utils::D_op::invert_2D(rD);
            for(size_t i = 0; i < 9; ++i){
                D[i] = psi_v*vD[i] + psi_r*rD[i];
            }
        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
            rD = utils::D_op::invert_3D(rD);
            for(size_t i = 0; i < 36; ++i){
                D[i] = psi_v*vD[i] + psi_r*rD[i];
            }
        }
    } else {
        for(size_t j = 0; j < v.size(); ++j){
            if(this->problem_type == utils::PROBLEM_TYPE_2D){
                const auto MD = this->materials[j]->stiffness_2D(e, p);
                const auto MS = utils::D_op::square_2D(MD);
                for(size_t i = 0; i < 9; ++i){
                    vD[i] += v[j]*MD[i];
                    rD[i] += v[j]*MS[i];
                }
            } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                const auto MD = this->materials[j]->stiffness_3D(e, p);
                const auto MS = utils::D_op::square_3D(MD);
                for(size_t i = 0; i < 36; ++i){
                    vD[i] += v[j]*MD[i];
                    rD[i] += v[j]*MS[i];
                }
            }
        }
        if(this->problem_type == utils::PROBLEM_TYPE_2D){
            rD = utils::D_op::square_root_2D(rD);
            for(size_t i = 0; i < 9; ++i){
                D[i] = psi_v*vD[i] + psi_r*rD[i];
            }
        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
            rD = utils::D_op::square_root_3D(rD);
            for(size_t i = 0; i < 36; ++i){
                D[i] = psi_v*vD[i] + psi_r*rD[i];
            }
        }
    }
}

void MultiMaterial::get_gradD_internal(std::vector<double>::const_iterator& rho, const double mix, const MeshElement* e, const gp_Pnt& p, std::vector<std::vector<double>>& gradD) const{
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

    size_t D_size = 0;
    std::vector<double> vD;
    std::vector<double> rD;
    if(this->problem_type == utils::PROBLEM_TYPE_2D){
        D_size = 9;
        vD.resize(D_size,0);
        rD.resize(D_size,0);
    } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
        D_size = 36;
        vD.resize(D_size,0);
        rD.resize(D_size,0);
    }

    if(mix >= 0){
        std::vector<double> rD_last;
        std::vector<double> rD_mult;
        // Reuss matrix for last material
        std::vector<double> MDl;
        if(this->problem_type == utils::PROBLEM_TYPE_2D){
            MDl = this->materials.back()->stiffness_2D(e, p);
            rD_last = this->materials.back()->stiffness_inverse_2D(e, p);
        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
            MDl = this->materials.back()->stiffness_3D(e, p);
            rD_last = this->materials.back()->stiffness_inverse_3D(e, p);
        }
        // Full Reuss matrix, necessary due to derivative chain rule
        {
            std::fill(rD.begin(), rD.end(), 0);
            for(size_t j = 0; j < v.size(); ++j){
                if(this->problem_type == utils::PROBLEM_TYPE_2D){
                    const auto MS = this->materials[j]->stiffness_inverse_2D(e, p);
                    for(size_t i = 0; i < 9; ++i){
                        rD[i] += v[j]*MS[i];
                    }
                } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                    const auto MS = this->materials[j]->stiffness_inverse_3D(e, p);
                    for(size_t i = 0; i < 36; ++i){
                        rD[i] += v[j]*MS[i];
                    }
                }
            }
            if(this->problem_type == utils::PROBLEM_TYPE_2D){
                rD_mult = utils::D_op::invert_2D(rD);
            } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                rD_mult = utils::D_op::invert_3D(rD);
            }
        }
        // Generate gradients
        std::vector<double> D_tmp(D_size, 0);
        for(size_t j = 0; j < v.size()-1; ++j){
            if(this->problem_type == utils::PROBLEM_TYPE_2D){
                const auto MD = this->materials[j]->stiffness_2D(e, p);
                const auto MS = this->materials[j]->stiffness_inverse_2D(e, p);
                for(size_t i = 0; i < 9; ++i){
                    D_tmp[i] = -(MS[i] - rD_last[i]);
                }
                rD = utils::D_op::mult_2D(rD_mult, D_tmp);
                rD = utils::D_op::mult_2D(rD, rD_mult);
                for(size_t i = 0; i < 9; ++i){
                    gradD[j+offset][i] = psi_v*(MD[i] - MDl[i]) + psi_r*rD[i];
                }
            } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                const auto MD = this->materials[j]->stiffness_3D(e, p);
                const auto MS = this->materials[j]->stiffness_inverse_3D(e, p);
                for(size_t i = 0; i < 36; ++i){
                    D_tmp[i] = -(MS[i] - rD_last[i]);
                }
                rD = utils::D_op::mult_3D(rD_mult, D_tmp);
                rD = utils::D_op::mult_3D(rD, rD_mult);
                for(size_t i = 0; i < 36; ++i){
                    gradD[j+offset][i] = psi_v*(MD[i] - MDl[i]) + psi_r*rD[i];
                }
            }
        }
    } else {
        std::vector<double> rD_last;
        std::vector<double> rD_mult;
        // Reuss matrix for last material
        std::vector<double> MDl;
        if(this->problem_type == utils::PROBLEM_TYPE_2D){
            MDl = this->materials.back()->stiffness_2D(e, p);
            rD_last = utils::D_op::square_2D(MDl);
        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
            MDl = this->materials.back()->stiffness_3D(e, p);
            rD_last = utils::D_op::square_3D(MDl);
        }
        // Full Reuss matrix, necessary due to derivative chain rule
        std::vector<std::vector<double>> D_cache(v.size());
        {
            std::fill(rD.begin(), rD.end(), 0);
            for(size_t j = 0; j < v.size(); ++j){
                if(this->problem_type == utils::PROBLEM_TYPE_2D){
                    D_cache[j] = this->materials[j]->stiffness_2D(e, p);
                    const auto MS = utils::D_op::square_2D(D_cache[j]);
                    for(size_t i = 0; i < 9; ++i){
                        rD[i] += v[j]*MS[i];
                    }
                } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                    D_cache[j] = this->materials[j]->stiffness_3D(e, p);
                    const auto MS = utils::D_op::square_3D(D_cache[j]);
                    for(size_t i = 0; i < 36; ++i){
                        rD[i] += v[j]*MS[i];
                    }
                }
            }
            if(this->problem_type == utils::PROBLEM_TYPE_2D){
                rD_mult = utils::D_op::invert_2D(utils::D_op::square_root_2D(rD));
                for(size_t i = 0; i < 9; ++i){
                    rD_mult[i] *= 0.5;
                }
                rD_last = utils::D_op::mult_2D(rD_mult, rD_last);
            } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                rD_mult = utils::D_op::invert_3D(utils::D_op::square_root_3D(rD));
                for(size_t i = 0; i < 36; ++i){
                    rD_mult[i] *= 0.5;
                }
                rD_last = utils::D_op::mult_3D(rD_mult, rD_last);
            }
        }
        // Generate gradients
        for(size_t j = 0; j < v.size()-1; ++j){
            if(this->problem_type == utils::PROBLEM_TYPE_2D){
                const auto& MD = D_cache[j];
                const auto MS = utils::D_op::square_2D(MD);
                rD = utils::D_op::mult_2D(rD_mult, MS);
                for(size_t i = 0; i < 9; ++i){
                    gradD[j+offset][i] = psi_v*(MD[i] - MDl[i]) + psi_r*(rD[i] - rD_last[i]);
                }
            } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                const auto& MD = D_cache[j];
                const auto MS = utils::D_op::square_3D(MD);
                rD = utils::D_op::mult_3D(rD_mult, MS);
                for(size_t i = 0; i < 36; ++i){
                    gradD[j+offset][i] = psi_v*(MD[i] - MDl[i]) + psi_r*(rD[i] - rD_last[i]);
                }
            }
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
