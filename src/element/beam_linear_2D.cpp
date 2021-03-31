/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   Foobar is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Foobar is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with SolidPrep.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "element/beam_linear_2D.hpp"
#include "logger.hpp"

namespace element{

BeamLinear2D::BeamLinear2D(BeamNode* n1, BeamNode* n2, std::vector<long> u_pos, double I, double A, double E):
    BeamElement(n1, n2, std::move(u_pos)), I(I), A(A), E(E){}

std::vector<double> BeamLinear2D::get_k() const{
    // A First Course in the Finite Element Method, Daryl L. Logan, 6 ed. (2016).
    
    double L = this->nodes[0]->point.Distance(this->nodes[1]->point)*1e-3;
    gp_Vec v(this->nodes[0]->point, this->nodes[1]->point);
    gp_Dir n(v);
    double theta = n.AngleWithRef(gp_Dir(1,0,0), gp_Dir(0,0,1));

    double C = std::cos(theta);
    double S = std::sin(theta);
    std::vector<double> k(21);
    k[0]  = (E/L)*(A*C*C + (12*I/(L*L))*S*S);
    k[1]  = (E/L)*(A - (12*I/(L*L)))*C*S;
    k[2]  = (E/L)*(A*S*S + (12*I/(L*L))*C*C);
    k[3]  = (E/L)*(-6*I/L)*S;
    k[4]  = (E/L)*(6*I/L)*C;
    k[5]  = (E/L)*4*I;
    k[6]  = (E/L)*(-(A*C*C + (12*I/(L*L))*S*S));
    k[7]  = (E/L)*(-(A - (12*I/(L*L)))*C*S);
    k[8]  = (E/L)*(6*I/L)*S;
    k[9]  = (E/L)*(A*C*C + (12*I/(L*L))*S*S);
    k[10] = (E/L)*(-(A - (12*I/(L*L)))*C*S);
    k[11] = (E/L)*(-(A*S*S + (12*I/(L*L))*C*C));
    k[12] = (E/L)*(-6*I/L)*C;
    k[13] = (E/L)*(A - (12*I/(L*L)))*C*S;
    k[14] = (E/L)*(A*S*S + (12*I/(L*L))*C*C);
    k[15] = (E/L)*(-6*I/L)*S;
    k[16] = (E/L)*(6*I/L)*C;
    k[17] = (E/L)*2*I;
    k[18] = (E/L)*(6*I/L)*S;
    k[19] = (E/L)*(-6*I/L)*C;
    k[20] = (E/L)*4*I;

    return k;
}

BeamNode* BeamLinear2D::get_internal_loads(size_t node, const std::vector<double>& u) const{
    logger::log_assert(node == 0 || node == 1, logger::ERROR, "wrong value for BeamLinear2D node, must be either 0 or 1.");

    std::vector<double> k = this->get_k();

    BeamNode2D* n = static_cast<BeamNode2D*>(this->nodes[node]);
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 6; ++j){
            n->results[i] += k[utils::to_triangular(node*3+i, j)]*u[this->u_pos[j]];
        }
        n->results[i] = std::abs(n->results[i]);
    }

    return this->get_node(node);
}

}
