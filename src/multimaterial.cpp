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

#include <array>
#include <cmath>
#include "multimaterial.hpp"
#include "utils.hpp"

MultiMaterial::MultiMaterial(std::vector<Material*> materials, utils::ProblemType type):
    materials(std::move(materials)), problem_type(type){}

void MultiMaterial::get_D(std::vector<double>::const_iterator& rho, const bool has_void, const double p, const double p_min, const double mix, std::vector<double>& D) const{
    double void_rho = 0;
    if(has_void){
        void_rho = p_min + (1-p_min)*std::pow(*rho, p);
        ++rho;
    }

    if(this->problem_type == utils::PROBLEM_TYPE_2D){
        std::array<double, 5> pos{0, 1, 3, 4, 8};
        if(this->materials.size() > 0){
            this->get_D_internal(rho, pos.data(), pos.size(), mix, D);
        } else {
            std::copy(this->materials[0]->stiffness_2D().begin(), this->materials[0]->stiffness_2D().end(), D.begin());
        }
        if(has_void){
            this->apply_void(void_rho, pos.data(), pos.size(), D);
        }
    } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
        std::array<double, 12> pos{0, 1, 2, 6, 7, 8, 12, 13, 14, 21, 28, 35};
        if(this->materials.size() > 0){
            this->get_D_internal(rho, pos.data(), pos.size(), mix, D);
        } else {
            std::copy(this->materials[0]->stiffness_3D().begin(), this->materials[0]->stiffness_3D().end(), D.begin());
        }
        if(has_void){
            this->apply_void(void_rho, pos.data(), pos.size(), D);
        }
    }
    rho += this->materials.size()-1;
}

void MultiMaterial::get_gradD(std::vector<double>::const_iterator& rho, const bool has_void, const double p, const double p_min, const double mix, std::vector<std::vector<double>>& gradD) const{
    double void_rho = 0;
    double void_drho = 0;
    if(has_void){
        void_rho = p_min + (1-p_min)*std::pow(*rho, p);
        void_drho = (1-p_min)*p*std::pow(*rho, p-1);
        ++rho;
    }

    if(this->problem_type == utils::PROBLEM_TYPE_2D){
        std::array<double, 5> pos{0, 1, 3, 4, 8};
        if(this->materials.size() > 0){
            this->get_gradD_internal(rho, has_void, pos.data(), pos.size(), mix, gradD);
        } else {
            std::copy(this->materials[0]->stiffness_2D().begin(), this->materials[0]->stiffness_2D().end(), gradD[0].begin());
        }
        if(has_void){
            this->get_D(rho, false, p, p_min, mix, gradD[0]);
            this->apply_void(void_drho, pos.data(), pos.size(), gradD[0]);
            for(size_t i = 1; i < gradD.size(); ++i){
                this->apply_void(void_rho, pos.data(), pos.size(), gradD[i]);
            }
        }
    } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
        std::array<double, 12> pos{0, 1, 2, 6, 7, 8, 12, 13, 14, 21, 28, 35};
        if(this->materials.size() > 0){
            this->get_gradD_internal(rho, has_void, pos.data(), pos.size(), mix, gradD);
        } else {
            std::copy(this->materials[0]->stiffness_3D().begin(), this->materials[0]->stiffness_3D().end(), gradD[0].begin());
        }
        if(has_void){
            this->get_D(rho, false, p, p_min, mix, gradD[0]);
            this->apply_void(void_drho, pos.data(), pos.size(), gradD[0]);
            for(size_t i = 1; i < gradD.size(); ++i){
                this->apply_void(void_rho, pos.data(), pos.size(), gradD[i]);
            }
        }
    }
    rho += this->materials.size()-1;
}

void MultiMaterial::get_D_internal(std::vector<double>::const_iterator& rho, const double* pos, const size_t posN, const double mix, std::vector<double>& D) const{
    std::vector<double> v(this->materials.size(), 0);
    double total = 0;
    for(size_t i = 0; i < v.size()-1; ++i){
        total += *(rho+i);
        v[i] = *(rho+i);
    }
    v[v.size()-1] = 1 - total;

    const double psi_v = std::abs(mix);
    const double psi_r = 1 - psi_v;

    if(mix >= 0){
        if(this->problem_type == utils::PROBLEM_TYPE_2D){
            for(size_t i = 0; i < posN; ++i){
                double voigt = 0;
                double reuss = 0;
                for(size_t j = 0; j < v.size(); ++j){
                    voigt += v[j]*this->materials[j]->stiffness_2D()[pos[i]];
                    reuss += v[j]/this->materials[j]->stiffness_2D()[pos[i]];
                }
                D[pos[i]] = psi_v*voigt + psi_r/reuss;
            }
        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
            for(size_t i = 0; i < posN; ++i){
                double voigt = 0;
                double reuss = 0;
                for(size_t j = 0; j < v.size(); ++j){
                    voigt += v[j]*this->materials[j]->stiffness_3D()[pos[i]];
                    reuss += v[j]/this->materials[j]->stiffness_3D()[pos[i]];
                }
                D[pos[i]] = psi_v*voigt + psi_r/reuss;
            }
        }
    } else {
        if(this->problem_type == utils::PROBLEM_TYPE_2D){
            for(size_t i = 0; i < posN; ++i){
                double voigt = 0;
                double reuss = 0;
                for(size_t j = 0; j < v.size(); ++j){
                    const double Ci = this->materials[j]->stiffness_2D()[pos[i]];
                    voigt += v[j]*Ci;
                    reuss += v[j]*Ci*Ci;
                }
                D[pos[i]] = psi_v*voigt + psi_r*std::sqrt(reuss);
            }
        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
            for(size_t i = 0; i < posN; ++i){
                double voigt = 0;
                double reuss = 0;
                for(size_t j = 0; j < v.size(); ++j){
                    const double Ci = this->materials[j]->stiffness_3D()[pos[i]];
                    voigt += v[j]*Ci;
                    reuss += v[j]*Ci*Ci;
                }
                D[pos[i]] = psi_v*voigt + psi_r*std::sqrt(reuss);
            }
        }
    }
}

void MultiMaterial::get_gradD_internal(std::vector<double>::const_iterator& rho, const bool has_void, const double* pos, const size_t posN, const double mix, std::vector<std::vector<double>>& gradD) const{
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

    if(mix >= 0){
        if(this->problem_type == utils::PROBLEM_TYPE_2D){
            for(size_t i = 0; i < posN; ++i){
                double reuss_tmp = 0;
                for(size_t j = 0; j < v.size(); ++j){
                    reuss_tmp += v[j]/this->materials[j]->stiffness_2D()[pos[i]];
                }
                for(size_t k = 0; k < v.size()-1; ++k){
                    double voigt = this->materials[k]->stiffness_2D()[pos[i]] - this->materials.back()->stiffness_2D()[pos[i]];
                    double reuss = -(1.0/this->materials[k]->stiffness_2D()[pos[i]] - 1.0/this->materials.back()->stiffness_2D()[pos[i]])/(reuss_tmp*reuss_tmp);

                    gradD[k+offset][pos[i]] = psi_v*voigt + psi_r*reuss;
                }
            }
        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
            for(size_t i = 0; i < posN; ++i){
                double reuss_tmp = 0;
                for(size_t j = 0; j < v.size(); ++j){
                    reuss_tmp += v[j]/this->materials[j]->stiffness_3D()[pos[i]];
                }
                for(size_t k = 0; k < v.size()-1; ++k){
                    const double voigt = this->materials[k]->stiffness_3D()[pos[i]] - this->materials.back()->stiffness_3D()[pos[i]];
                    const double reuss = -(1.0/this->materials[k]->stiffness_3D()[pos[i]] - 1.0/this->materials.back()->stiffness_3D()[pos[i]])/(reuss_tmp*reuss_tmp);

                    gradD[k+offset][pos[i]] = psi_v*voigt + psi_r*reuss;
                }
            }
        }
    } else {
        // TODO
        if(this->problem_type == utils::PROBLEM_TYPE_2D){
            for(size_t i = 0; i < posN; ++i){
                double reuss_tmp = 0;
                for(size_t j = 0; j < v.size(); ++j){
                    const double Ci = this->materials[j]->stiffness_2D()[pos[i]];
                    reuss_tmp += v[j]*Ci*Ci;
                }
                for(size_t k = 0; k < v.size()-1; ++k){
                    const double Ci = this->materials[k]->stiffness_2D()[pos[i]];
                    const double Cb = this->materials.back()->stiffness_2D()[pos[i]];
                    const double voigt = Ci - Cb;
                    const double reuss = (Ci*Ci - Cb*Cb)/(2*std::sqrt(reuss_tmp*reuss_tmp));;

                    gradD[k+offset][pos[i]] = psi_v*voigt + psi_r*reuss;
                }
            }
        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
            for(size_t i = 0; i < posN; ++i){
                double reuss_tmp = 0;
                for(size_t j = 0; j < v.size(); ++j){
                    const double Ci = this->materials[j]->stiffness_3D()[pos[i]];
                    reuss_tmp += v[j]*Ci*Ci;
                }
                for(size_t k = 0; k < v.size()-1; ++k){
                    const double Ci = this->materials[k]->stiffness_3D()[pos[i]];
                    const double Cb = this->materials.back()->stiffness_3D()[pos[i]];
                    const double voigt = Ci - Cb;
                    const double reuss = (Ci*Ci - Cb*Cb)/(2*std::sqrt(reuss_tmp*reuss_tmp));;

                    gradD[k+offset][pos[i]] = psi_v*voigt + psi_r*reuss;
                }
            }
        }
    }
}

double MultiMaterial::get_density(std::vector<double>::const_iterator& rho, const bool has_void) const{
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

    rho += v.size();

    double d = 0;
    for(size_t i = 0; i < v.size(); ++i){
        d += v[i]*this->materials[i]->density;
    }

    d *= void_rho;

    return d;
}

double MultiMaterial::get_density_deriv(std::vector<double>::const_iterator& rho, const bool has_void, std::vector<double>::iterator& grad) const{
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

    rho += v.size();

    double d = 0;
    for(size_t i = 0; i < v.size(); ++i){
        d += v[i]*this->materials[i]->density;
    }

    if(has_void){
        *grad = d;
        ++grad;
    }
    for(size_t i = 0; i < v.size()-1; ++i){
        *grad = void_rho*(this->materials[i]->density - this->materials.back()->density);
        ++grad;
    }

    return void_rho*d;
}
