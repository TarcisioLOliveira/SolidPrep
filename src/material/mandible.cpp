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

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <fstream>
#include <algorithm>
#include <gp_Ax1.hxx>
#include "logger.hpp"
#include "utils/D_operations.hpp"
#include "material/mandible.hpp"
#include "utils/basis_tensor.hpp"

namespace material{

Mandible::Mandible(const std::string& name, Material* outer, Material* inner, const std::string& path_points1, const std::string& path_points2, double C, bool with_implant, ImplantRegion imp):
    Material(name, std::vector<double>(1), std::vector<double>(1)),
    outer(outer), inner(inner), C(C), ring1(this), ring2(this), with_implant(with_implant), implant(imp){

    {
        auto p1 = this->load_points(path_points1);
        auto p2 = this->load_points(path_points2);
        this->ring1.initialize_center(p1);
        this->ring2.initialize_center(p2);
        this->center_normal = gp_Vec(this->ring1.center, this->ring2.center);
        gp_Vec dist(this->ring1.center, p1[0]);
        gp_Vec proj = dist.Dot(this->center_normal)*center_normal;
        this->center_ref = dist - proj;
        this->ring1.initialize_points(p1);
        this->ring2.initialize_points(p2);
    }

}
std::vector<double> Mandible::stiffness_2D(const MeshElement* const e, const gp_Pnt& p) const{
    auto Do = this->outer->stiffness_2D(e, p);
    auto Di = this->inner->stiffness_2D(e, p);
    auto coeff = this->get_multiplier(p);

    if(!this->with_implant){
        for(size_t i = 0; i < 9; ++i){
            Do[i] = (1 - coeff)*Do[i] + coeff*Di[i];
        }
    } else {
        const auto m = this->implant.get_implant_multiplier(p);
        for(size_t i = 0; i < 9; ++i){
            Do[i] = m*((1 - coeff)*Do[i] + coeff*Di[i]);
        }
    }
    return Do;
}
std::vector<double> Mandible::stiffness_3D(const MeshElement* const e, const gp_Pnt& p) const{
    auto Do = this->outer->stiffness_3D(e, p);
    auto Di = this->inner->stiffness_3D(e, p);
    auto coeff = this->get_multiplier(p);

    if(!this->with_implant){
        for(size_t i = 0; i < 36; ++i){
            Do[i] = (1 - coeff)*Do[i] + coeff*Di[i];
        }
    } else {
        const auto m = this->implant.get_implant_multiplier(p);
        for(size_t i = 0; i < 36; ++i){
            Do[i] = m*((1 - coeff)*Do[i] + coeff*Di[i]);
        }
    }
    return Do;
}
std::vector<double> Mandible::stiffness_inverse_2D(const MeshElement* const e, const gp_Pnt& p) const{
    auto coeff = this->get_multiplier(p);
    constexpr double TOL = 1e-7;
    if(coeff < TOL){
        if(!this->with_implant){
            return this->outer->stiffness_inverse_2D(e, p);
        } else {
            auto S = this->outer->stiffness_inverse_2D(e, p);
            const auto m = this->implant.get_implant_multiplier(p);
            for(size_t i = 0; i < 9; ++i){
                S[i] /= m;
            }
            return S;
        }
    } else if(coeff > 1 - TOL){
        if(!this->with_implant){
            return this->inner->stiffness_inverse_2D(e, p);
        } else {
            auto S = this->inner->stiffness_inverse_2D(e, p);
            const auto m = this->implant.get_implant_multiplier(p);
            for(size_t i = 0; i < 9; ++i){
                S[i] /= m;
            }
            return S;
        }
    } else {
        return utils::D_op::invert_2D(this->stiffness_2D(e, p));
    }
}
std::vector<double> Mandible::stiffness_inverse_3D(const MeshElement* const e, const gp_Pnt& p) const{
    auto coeff = this->get_multiplier(p);
    constexpr double TOL = 1e-7;
    if(coeff < TOL){
        if(!this->with_implant){
            return this->outer->stiffness_inverse_3D(e, p);
        } else {
            auto S = this->outer->stiffness_inverse_3D(e, p);
            const auto m = this->implant.get_implant_multiplier(p);
            for(size_t i = 0; i < 36; ++i){
                S[i] /= m;
            }
            return S;
        }
    } else if(coeff > 1 - TOL){
        if(!this->with_implant){
            return this->inner->stiffness_inverse_3D(e, p);
        } else {
            auto S = this->inner->stiffness_inverse_3D(e, p);
            const auto m = this->implant.get_implant_multiplier(p);
            for(size_t i = 0; i < 36; ++i){
                S[i] /= m;
            }
            return S;
        }
    } else {
        return utils::D_op::invert_3D(this->stiffness_3D(e, p));
    }
}
double Mandible::get_density(const MeshElement* const e, const gp_Pnt& p) const{
    double d_o = this->outer->get_density(e, p);
    double d_i = this->inner->get_density(e, p);
    auto coeff = this->get_multiplier(p);

    if(!this->with_implant){
        return (1 - coeff)*d_o + coeff*d_i;
    } else {
        const auto m = this->implant.get_implant_multiplier(p);
        return m*((1 - coeff)*d_o + coeff*d_i);
    }
}
double Mandible::beam_E_2D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 2, 2>& R) const{
    auto S = this->stiffness_inverse_2D(e, p);
    const auto Rt = utils::basis_tensor_2D_inv_T(R);
    this->rotate_S_2D(S, Rt);

    return 1/S[0];
}
double Mandible::beam_E_3D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 3, 3>& R) const{
    auto S = this->stiffness_inverse_3D(e, p);
    const auto Rt = utils::basis_tensor_3D_inv_T(R);
    this->rotate_S_3D(S, Rt);

    return 1.0/S[0];
}
std::array<double, 2> Mandible::beam_EG_2D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 2, 2>& R) const{
    auto S = this->stiffness_inverse_2D(e, p);
    const auto Rt = utils::basis_tensor_2D_inv_T(R);
    this->rotate_S_2D(S, Rt);

    return {1.0/S[0], 1.0/S[8]};
}
std::array<double, 4> Mandible::beam_EG_3D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 3, 3>& R) const{
    auto S = this->stiffness_inverse_3D(e, p);
    const auto Rt = utils::basis_tensor_3D_inv_T(R);
    this->rotate_S_3D(S, Rt);

    return {1.0/S[0], 1.0/S[21], 1.0/S[28], 1.0/S[35]};
}
double Mandible::S12_2D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 2, 2>& R) const{
    auto EG_o = this->outer->S12_2D(e, p, R);
    auto EG_i = this->inner->S12_2D(e, p, R);
    auto coeff = this->get_multiplier(p);

    if(!this->with_implant){
        return (1 - coeff)*EG_o + coeff*EG_i;
    } else {
        const auto m = this->implant.get_implant_multiplier(p);
        return m*((1 - coeff)*EG_o + coeff*EG_i);
    }
}
std::array<double, 2> Mandible::S12_S13_3D(const MeshElement* const e, const gp_Pnt& p, const Eigen::Matrix<double, 3, 3>& R) const{
    auto EG_o = this->outer->S12_S13_3D(e, p, R);
    auto EG_i = this->inner->S12_S13_3D(e, p, R);
    auto coeff = this->get_multiplier(p);

    if(!this->with_implant){
        for(size_t i = 0; i < 2; ++i){
            EG_o[i] = (1 - coeff)*EG_o[i] + coeff*EG_i[i];
        }
    } else {
        const auto m = this->implant.get_implant_multiplier(p);
        for(size_t i = 0; i < 2; ++i){
            EG_o[i] = m*((1 - coeff)*EG_o[i] + coeff*EG_i[i]);
        }
    }
    return EG_o;
}
    
void Mandible::Ring::initialize_center(const std::vector<gp_Pnt>& init){
    double X = 0, Y = 0, Z = 0;
    for(const auto& p:init){
        X += p.X();
        Y += p.Y();
        Z += p.Z();
    }
    X /= init.size();
    Y /= init.size();
    Z /= init.size();
    this->center.SetCoord(X, Y, Z);
    gp_Vec v1(center, init[0]);
    gp_Vec v2(center, init[1]);
    this->normal = v1.Crossed(v2);
}
void Mandible::Ring::initialize_points(const std::vector<gp_Pnt>& init){
    this->points.clear();
    this->points.reserve(init.size());
    for(const auto& p:init){
        this->points.push_back(M->to_ring_point_cartesian(p));
    }
    std::sort(this->points.begin(), this->points.end());
}

Mandible::RingPoint Mandible::Ring::get_r_max(double theta) const{
    auto pos = std::find_if(this->points.begin(), this->points.end(), 
            [theta](const RingPointCartesian& p){return p.theta > theta;});

    size_t top = (pos - this->points.begin()) % this->points.size();
    size_t bot = (top + this->points.size() - 1) % this->points.size();

    auto p1 = this->points[bot];
    auto p2 = this->points[top];
    double theta1 = p1.theta;
    double theta2 = p2.theta;
    if(theta1 > theta2){
        theta1 -= 2*M_PI;
    }
    if(theta > theta2){
        theta -= 2*M_PI;
    }

    double coeff = (theta - theta1)/(theta2 - theta1);
    gp_Pnt pn = p1.p;
    pn.Translate(coeff*gp_Vec(p1.p, p2.p));
    //auto p = this->M->to_ring_point(pn);
    RingPoint p;
    p.theta = theta;
    p.lambda = p1.lambda + (p2.lambda - p1.lambda)*coeff;
    gp_Pnt c2 = this->M->ring1.center.Translated(p.lambda*this->M->center_normal);
    p.r = gp_Vec(c2, pn).Magnitude();

    return p;
}

double Mandible::get_multiplier(const gp_Pnt& p) const{
    auto rp = this->to_ring_point(p);
    auto r1 = this->ring1.get_r_max(rp.theta);
    auto r2 = this->ring2.get_r_max(rp.theta);
    double coeff = (rp.lambda - r1.lambda)/(r2.lambda - r1.lambda);

    double r_max = r1.r + (r2.r - r1.r)*coeff;

    return this->heaviside(rp.r, r_max);
}

std::vector<gp_Pnt> Mandible::load_points(const std::string& path) const{
    std::vector<gp_Pnt> points;

    std::ifstream file;
    file.open(path);
    logger::log_assert(file.good(), logger::ERROR, "File not found: {}", path.c_str());
    std::string val;
    while(std::getline(file, val, ' ')){
        double X = std::stod(val);
        std::getline(file, val, ' ');
        double Y = std::stod(val);
        std::getline(file, val, '\n');
        double Z = std::stod(val);

        points.emplace_back(X, Y, Z);
    }
    points.shrink_to_fit();

    return points;
}

Mandible::RingPoint Mandible::to_ring_point(const gp_Pnt& p) const{
    RingPoint rp;
    gp_Vec v1(this->ring1.center, p);

    rp.lambda = v1.Dot(this->center_normal);
    gp_Pnt p2 = this->ring1.center.Translated(rp.lambda*this->center_normal);
    rp.r = p.Distance(p2);
    gp_Vec distv(p, p2);
    rp.theta = this->angle_with_ref(distv);

    return rp;
}
Mandible::RingPointCartesian Mandible::to_ring_point_cartesian(const gp_Pnt& p) const{
    RingPointCartesian rp;
    gp_Vec v1(this->ring1.center, p);

    rp.lambda = v1.Dot(this->center_normal);
    gp_Pnt p2 = this->ring1.center.Translated(rp.lambda*this->center_normal);
    rp.p = p;
    gp_Vec distv(p, p2);
    rp.theta = this->angle_with_ref(distv);

    return rp;
}

double Mandible::angle_with_ref(const gp_Vec& dist) const{
    const gp_Vec v(dist.Normalized());
    const gp_Vec n = (v - (v.Dot(this->center_normal)*this->center_normal)).Normalized();
    const gp_Vec cv = n.Crossed(this->center_ref);
    const double mult = cv.Dot(this->center_normal);

    const int sgn = (0 <= mult) - (mult < 0);

    return std::acos(n.Dot(this->center_ref))*sgn + M_PI;
}

Mandible::ImplantRegion::ImplantRegion(const gp_Pnt& center_1, const gp_Pnt& center_2, double r1, double r2, const std::vector<double>& a, double dl):
    center_1(center_1), center_2(center_2), r1(r1), r2(r2),
    normal(gp_Vec(center_1, center_2)),
    decay_distance(dl), min_str(a[0]),
    a(a), a_len(a.size()),
    max_l(center_1.Distance(center_2)){
    
}

double Mandible::ImplantRegion::get_implant_multiplier(const gp_Pnt& p) const{
    const gp_Vec v1(this->center_1, p);

    const double l = v1.Dot(this->normal);
    const gp_Vec v2 = v1 - l*this->normal;
    const double dist = v2.Magnitude();
    if(l >= 0 && l <= max_l){
        const double r = (r1 - r2)*l/max_l + r1;
        if(dist < r){
            return min_str;
        } else if(dist > r + decay_distance){
            return 1;
        } else {
            const double x = dist - r;
            return this->f(x);
        }
    } else if(l < 0){
        if(dist <= r1){
            if(-l > decay_distance){
                return 1;
            } else {
                return this->f(-l);
            }
        } else {
            const double dr = dist - r1;
            const double dc = std::sqrt(l*l + dr*dr);
            if(dc > decay_distance){
                return 1;
            } else {
                return this->f(dc);
            }
        }
    } else if(l > max_l){
        const double dl = l - max_l;
        if(dist <= r2){
            if(dl > decay_distance){
                return 1;
            } else {
                return this->f(dl);
            }
        } else {
            const double dr = dist - r2;
            const double dc = std::sqrt(dl*dl + dr*dr);
            if(dc > decay_distance){
                return 1;
            } else {
                return this->f(dc);
            }
        }
    }
    return 1;
}

}
