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
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "logger.hpp"
#include "multimaterial.hpp"
#include "utils.hpp"

MultiMaterial::MultiMaterial(std::vector<Material*> materials, utils::ProblemType type, bool has_void):
    materials(std::move(materials)), problem_type(type), has_void(has_void){}

void MultiMaterial::get_D(std::vector<double>::const_iterator rho, const double mix, std::vector<double>& D) const{
    if(has_void){
        ++rho;
    }

    if(this->problem_type == utils::PROBLEM_TYPE_2D){
        std::array<double, 5> pos{0, 1, 3, 4, 8};
        if(this->materials.size() > 1){
            this->get_D_internal(rho, pos.data(), pos.size(), mix, D);
        } else {
            std::copy(this->materials[0]->stiffness_2D().begin(), this->materials[0]->stiffness_2D().end(), D.begin());
        }
    } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
        std::array<double, 12> pos{0, 1, 2, 6, 7, 8, 12, 13, 14, 21, 28, 35};
        if(this->materials.size() > 1){
            this->get_D_internal(rho, pos.data(), pos.size(), mix, D);
        } else {
            std::copy(this->materials[0]->stiffness_3D().begin(), this->materials[0]->stiffness_3D().end(), D.begin());
        }
    }
}

void MultiMaterial::get_gradD(std::vector<double>::const_iterator rho, const double mix, std::vector<std::vector<double>>& gradD) const{
    if(has_void){
        ++rho;
    }

    if(this->problem_type == utils::PROBLEM_TYPE_2D){
        std::array<double, 5> pos{0, 1, 3, 4, 8};
        if(this->materials.size() > 1){
            this->get_gradD_internal(rho, pos.data(), pos.size(), mix, gradD);
        }
        if(has_void){
            std::copy(this->materials[0]->stiffness_2D().begin(), this->materials[0]->stiffness_2D().end(), gradD[0].begin());
        }
    } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
        std::array<double, 12> pos{0, 1, 2, 6, 7, 8, 12, 13, 14, 21, 28, 35};
        if(this->materials.size() > 1){
            this->get_gradD_internal(rho, pos.data(), pos.size(), mix, gradD);
        }
        if(has_void){
            std::copy(this->materials[0]->stiffness_3D().begin(), this->materials[0]->stiffness_3D().end(), gradD[0].begin());
        }
    }
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
                auto& MD = this->materials[j]->stiffness_2D();
                auto& MS = this->materials[j]->stiffness_inverse_2D();
                for(size_t i = 0; i < posN; ++i){
                    vD[pos[i]] += v[j]*MD[pos[i]];
                    rD[pos[i]] += v[j]*MS[pos[i]];
                }
            } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                auto& MD = this->materials[j]->stiffness_3D();
                auto& MS = this->materials[j]->stiffness_inverse_3D();
                for(size_t i = 0; i < posN; ++i){
                    vD[pos[i]] += v[j]*MD[pos[i]];
                    rD[pos[i]] += v[j]*MS[pos[i]];
                }
            }
        }
        if(this->problem_type == utils::PROBLEM_TYPE_2D){
            rD = this->invert_2D(rD);
        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
            rD = this->invert_3D(rD);
        }
        for(size_t i = 0; i < posN; ++i){
            D[pos[i]] = psi_v*vD[pos[i]] + psi_r*rD[pos[i]];
        }
    } else {
        for(size_t j = 0; j < v.size(); ++j){
            if(this->problem_type == utils::PROBLEM_TYPE_2D){
                auto& MD = this->materials[j]->stiffness_2D();
                auto MS = this->square_2D(this->materials[j]->stiffness_2D());
                for(size_t i = 0; i < posN; ++i){
                    vD[pos[i]] += v[j]*MD[pos[i]];
                    rD[pos[i]] += v[j]*MS[pos[i]];
                }
            } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                auto& MD = this->materials[j]->stiffness_3D();
                auto MS = this->square_3D(this->materials[j]->stiffness_3D());
                for(size_t i = 0; i < posN; ++i){
                    vD[pos[i]] += v[j]*MD[pos[i]];
                    rD[pos[i]] += v[j]*MS[pos[i]];
                }
            }
        }
        if(this->problem_type == utils::PROBLEM_TYPE_2D){
            rD = this->square_root_2D(rD);
        } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
            rD = this->square_root_3D(rD);
        }
        for(size_t i = 0; i < posN; ++i){
            D[pos[i]] = psi_v*vD[pos[i]] + psi_r*rD[pos[i]];
        }
    }
}

void MultiMaterial::get_gradD_internal(std::vector<double>::const_iterator& rho, const double* pos, const size_t posN, const double mix, std::vector<std::vector<double>>& gradD) const{
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
        {
            const size_t j = v.size() - 1;
            if(this->problem_type == utils::PROBLEM_TYPE_2D){
                rD_last = this->materials[j]->stiffness_inverse_2D();
            } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                rD_last = this->materials[j]->stiffness_inverse_3D();
            }
        }
        // Full Reuss matrix, necessary due to derivative chain rule
        {
            std::fill(rD.begin(), rD.end(), 0);
            for(size_t j = 0; j < v.size(); ++j){
                if(this->problem_type == utils::PROBLEM_TYPE_2D){
                    auto& MS = this->materials[j]->stiffness_inverse_2D();
                    for(size_t i = 0; i < posN; ++i){
                        rD[pos[i]] += v[j]*MS[pos[i]];
                    }
                } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                    auto& MS = this->materials[j]->stiffness_inverse_3D();
                    for(size_t i = 0; i < posN; ++i){
                        rD[pos[i]] += v[j]*MS[pos[i]];
                    }
                }
            }
            if(this->problem_type == utils::PROBLEM_TYPE_2D){
                rD_mult = this->square_2D(this->invert_2D(rD));
                for(size_t i = 0; i < posN; ++i){
                    rD_mult[pos[i]] *= -1.0;
                }
                rD_last = this->mult_2D(rD_mult, rD_last);
            } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                rD_mult = this->square_3D(this->invert_3D(rD));
                for(size_t i = 0; i < posN; ++i){
                    rD_mult[pos[i]] *= -1.0;
                }
                rD_last = this->mult_3D(rD_mult, rD_last);
            }
        }
        // Generate gradients
        for(size_t j = 0; j < v.size()-1; ++j){
            if(this->problem_type == utils::PROBLEM_TYPE_2D){
                auto& MD = this->materials[j]->stiffness_2D();
                auto& MDl = this->materials.back()->stiffness_2D();
                auto& MS = this->materials[j]->stiffness_inverse_2D();
                rD = this->mult_2D(rD_mult, MS);
                for(size_t i = 0; i < posN; ++i){
                    gradD[j+offset][pos[i]] = psi_v*(MD[pos[i]] - MDl[pos[i]]) + psi_r*(rD[pos[i]] - rD_last[pos[i]]);
                }
            } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                auto& MD = this->materials[j]->stiffness_3D();
                auto& MDl = this->materials.back()->stiffness_3D();
                auto& MS = this->materials[j]->stiffness_inverse_3D();
                rD = this->mult_3D(rD_mult, MS);
                for(size_t i = 0; i < posN; ++i){
                    gradD[j+offset][pos[i]] = psi_v*(MD[pos[i]] - MDl[pos[i]]) + psi_r*(rD[pos[i]] - rD_last[pos[i]]);
                }
            }
        }
    } else {
        std::vector<double> rD_last;
        std::vector<double> rD_mult;
        // Reuss matrix for last material
        {
            const size_t j = v.size() - 1;
            if(this->problem_type == utils::PROBLEM_TYPE_2D){
                rD_last = this->square_2D(this->materials[j]->stiffness_2D());
            } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                rD_last = this->square_3D(this->materials[j]->stiffness_3D());
            }
        }
        // Full Reuss matrix, necessary due to derivative chain rule
        {
            std::fill(rD.begin(), rD.end(), 0);
            for(size_t j = 0; j < v.size(); ++j){
                if(this->problem_type == utils::PROBLEM_TYPE_2D){
                    auto MS = this->square_2D(this->materials[j]->stiffness_2D());
                    for(size_t i = 0; i < posN; ++i){
                        rD[pos[i]] += v[j]*MS[pos[i]];
                    }
                } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                    auto MS = this->square_3D(this->materials[j]->stiffness_3D());
                    for(size_t i = 0; i < posN; ++i){
                        rD[pos[i]] += v[j]*MS[pos[i]];
                    }
                }
            }
            if(this->problem_type == utils::PROBLEM_TYPE_2D){
                rD_mult = this->invert_2D(this->square_root_2D(rD));
                for(size_t i = 0; i < posN; ++i){
                    rD_mult[pos[i]] *= 0.5;
                }
                rD_last = this->mult_2D(rD_mult, rD_last);
            } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                rD_mult = this->invert_3D(this->square_root_3D(rD));
                for(size_t i = 0; i < posN; ++i){
                    rD_mult[pos[i]] *= 0.5;
                }
                rD_last = this->mult_3D(rD_mult, rD_last);
            }
        }
        // Generate gradients
        for(size_t j = 0; j < v.size()-1; ++j){
            if(this->problem_type == utils::PROBLEM_TYPE_2D){
                auto& MD = this->materials[j]->stiffness_2D();
                auto& MDl = this->materials.back()->stiffness_2D();
                auto MS = this->square_2D(this->materials[j]->stiffness_2D());
                rD = this->mult_2D(rD_mult, MS);
                for(size_t i = 0; i < posN; ++i){
                    gradD[j+offset][pos[i]] = psi_v*(MD[pos[i]] - MDl[pos[i]]) + psi_r*(rD[pos[i]] - rD_last[pos[i]]);
                }
            } else if(this->problem_type == utils::PROBLEM_TYPE_3D){
                auto& MD = this->materials[j]->stiffness_3D();
                auto& MDl = this->materials.back()->stiffness_3D();
                auto MS = this->square_3D(this->materials[j]->stiffness_3D());
                rD = this->mult_3D(rD_mult, MS);
                for(size_t i = 0; i < posN; ++i){
                    gradD[j+offset][pos[i]] = psi_v*(MD[pos[i]] - MDl[pos[i]]) + psi_r*(rD[pos[i]] - rD_last[pos[i]]);
                }
            }
        }
    }
}

double MultiMaterial::get_density(std::vector<double>::const_iterator& rho) const{
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
        d += v[i]*this->materials[i]->density;
    }

    d *= void_rho;

    return d;
}

double MultiMaterial::get_density_deriv(std::vector<double>::const_iterator& rho, std::vector<double>::iterator& grad) const{
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

std::vector<double> MultiMaterial::invert_2D(const std::vector<double>& d) const{
    std::vector<double> D = {
    d[4]/(d[0]*d[4] - d[1]*d[3])
    ,
    -d[1]/(d[0]*d[4] - d[1]*d[3])
    ,
    0
    ,
    -d[3]/(d[0]*d[4] - d[1]*d[3])
    ,
    d[0]/(d[0]*d[4] - d[1]*d[3])
    ,
    0
    ,
    0
    ,
    0
    ,
    1/d[8]
    };
    return D;
}
std::vector<double> MultiMaterial::invert_3D(const std::vector<double>& d) const{
    std::vector<double> D = {
    (d[13]*d[8] - d[14]*d[7])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    (d[1]*d[14] - d[13]*d[2])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    (-d[1]*d[8] + d[2]*d[7])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    0
    ,
    0
    ,
    0
    ,
    (-d[12]*d[8] + d[14]*d[6])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    (-d[0]*d[14] + d[12]*d[2])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    (d[0]*d[8] - d[2]*d[6])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    0
    ,
    0
    ,
    0
    ,
    (d[12]*d[7] - d[13]*d[6])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    (d[0]*d[13] - d[1]*d[12])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    (-d[0]*d[7] + d[1]*d[6])/(d[0]*d[13]*d[8] - d[0]*d[14]*d[7] - d[1]*d[12]*d[8] + d[1]*d[14]*d[6] + d[12]*d[2]*d[7] - d[13]*d[2]*d[6])
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    1/d[21]
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    1/d[28]
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    1/d[35]
    };

    return D;
}
std::vector<double> MultiMaterial::square_2D(const std::vector<double>& d) const{
    std::vector<double> D = {
    d[0]*d[0] + d[1]*d[3]
    ,
    d[1]*(d[0] + d[4])
    ,
    0
    ,
    d[3]*(d[0] + d[4])
    ,
    d[1]*d[3] + d[4]*d[4]
    ,
    0
    ,
    0
    ,
    0
    ,
    d[8]*d[8]
    };
    return D;
}
std::vector<double> MultiMaterial::square_3D(const std::vector<double>& d) const{
    std::vector<double> D = {
    d[0]*d[0] + d[1]*d[6] + d[12]*d[2]
    ,
    d[0]*d[1] + d[1]*d[7] + d[13]*d[2]
    ,
    d[0]*d[2] + d[1]*d[8] + d[14]*d[2]
    ,
    0
    ,
    0
    ,
    0
    ,
    d[0]*d[6] + d[12]*d[8] + d[6]*d[7]
    ,
    d[1]*d[6] + d[13]*d[8] + d[7]*d[7]
    ,
    d[14]*d[8] + d[2]*d[6] + d[7]*d[8]
    ,
    0
    ,
    0
    ,
    0
    ,
    d[0]*d[12] + d[12]*d[14] + d[13]*d[6]
    ,
    d[1]*d[12] + d[13]*d[14] + d[13]*d[7]
    ,
    d[12]*d[2] + d[13]*d[8] + d[14]*d[14]
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    d[21]*d[21]
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    d[28]*d[28]
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    d[35]*d[35]
    };

    return D;
}
std::vector<double> MultiMaterial::mult_2D(const std::vector<double>& d, const std::vector<double>& e) const{
    std::vector<double> D = {
    d[0]*e[0] + d[1]*e[3]
    ,
    d[0]*e[1] + d[1]*e[4]
    ,
    0
    ,
    d[3]*e[0] + d[4]*e[3]
    ,
    d[3]*e[1] + d[4]*e[4]
    ,
    0
    ,
    0
    ,
    0
    ,
    d[8]*e[8]
    };

    return D;
}
std::vector<double> MultiMaterial::mult_3D(const std::vector<double>& d, const std::vector<double>& e) const{
    std::vector<double> D = {
    d[0]*e[0] + d[1]*e[6] + d[2]*e[12]
    ,
    d[0]*e[1] + d[1]*e[7] + d[2]*e[13]
    ,
    d[0]*e[2] + d[1]*e[8] + d[2]*e[14]
    ,
    0
    ,
    0
    ,
    0
    ,
    d[6]*e[0] + d[7]*e[6] + d[8]*e[12]
    ,
    d[6]*e[1] + d[7]*e[7] + d[8]*e[13]
    ,
    d[6]*e[2] + d[7]*e[8] + d[8]*e[14]
    ,
    0
    ,
    0
    ,
    0
    ,
    d[12]*e[0] + d[13]*e[6] + d[14]*e[12]
    ,
    d[12]*e[1] + d[13]*e[7] + d[14]*e[13]
    ,
    d[12]*e[2] + d[13]*e[8] + d[14]*e[14]
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    d[21]*e[21]
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    d[28]*e[28]
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    0
    ,
    d[35]*e[35]
    };

    return D;
}
std::vector<double> MultiMaterial::square_root_2D(const std::vector<double>& d) const{
    std::vector<double> D = {
    std::sqrt(2)*d[1]*d[3]*((-d[0] + d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*std::sqrt(d[0] + d[4] - std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4])) + (d[0] - d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*std::sqrt(d[0] + d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4])))/((-d[0] + d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*(d[0] - d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))
    ,
    std::sqrt(2)*d[1]*(2*d[1]*d[3]*(d[0] - d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*std::sqrt(d[0] + d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4])) - (-d[0] + d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*std::sqrt(d[0] + d[4] - std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*(d[0]*d[0] - 2*d[0]*d[4] + d[0]*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]) + 2*d[1]*d[3] + d[4]*d[4] - d[4]*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4])))/((-d[0] + d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*(d[0] - d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*(d[0]*d[0] - 2*d[0]*d[4] + d[0]*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]) + 4*d[1]*d[3] + d[4]*d[4] - d[4]*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4])))
    ,
    0
    ,
    std::sqrt(2)*d[3]*(-std::sqrt(d[0] + d[4] - std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4])) + std::sqrt(d[0] + d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4])))/(2*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))
    ,
    std::sqrt(2)*(d[1]*d[3]*std::sqrt(d[0] + d[4] + std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4])) + std::sqrt(d[0] + d[4] - std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))*(d[0]*d[0] - 2*d[0]*d[4] + d[0]*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]) + 2*d[1]*d[3] + d[4]*d[4] - d[4]*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]))/2)/(d[0]*d[0] - 2*d[0]*d[4] + d[0]*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] + d[4]*d[4]) + 4*d[1]*d[3] + d[4]*d[4] - d[4]*std::sqrt(d[0]*d[0] - 2*d[0]*d[4] + 4*d[1]*d[3] +d[4]*d[4]))
    ,
    0
    ,
    0
    ,
    0
    ,
    std::sqrt(d[8])
    };

    return D;
}
std::vector<double> MultiMaterial::square_root_3D(const std::vector<double>& d) const{
    typedef Eigen::Matrix<double, 6, 6> MatType;
    MatType M = Eigen::Map<const MatType>(d.data(), 6, 6);

    Eigen::SelfAdjointEigenSolver<MatType> eigen(M);

    M = eigen.operatorSqrt();

    std::vector<double> result(36);
    std::copy(M.data(), M.data()+36, result.begin());

    return result;
}
