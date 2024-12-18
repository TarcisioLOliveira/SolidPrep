/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
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

#include "beam_graph.hpp"
#include "element_factory.hpp"
#include "logger.hpp"
#include "math/matrix.hpp"
#include "utils.hpp"

std::vector<double> BeamGraph::run(){
    // Uses UPLO=L lower symmetric band matrices

    std::vector<std::vector<gp_Pnt>> beams;
    size_t k_dim = BeamElementFactory::get_k_dimension(this->type);

    // K dimensions
    int W = 0;
    int N = k_dim;

    // Get beam paths
    int node_qnt = 0;
    for(auto& f:this->data->forces){
        for(auto& s:this->data->supports){
            std::vector<gp_Pnt> b(this->data->pathfinder->find_path(f.S, s.S));
            W += (k_dim/2)*b.size() - s.DOF();
            node_qnt += b.size();
            beams.push_back(std::move(b));
        }
    }

    // Initialize matrices and vectors
    std::vector<double> K(W*N);
    std::vector<double> F(W);
    std::vector<BeamElement*> elems(node_qnt-1);
    this->clear_nodes();
    this->nodes.resize(node_qnt);

    BeamNodeFactory::BeamNodeType node_type = BeamElementFactory::get_node_type(this->type);

    size_t cur_id = 0;
    size_t pos_offset = 0;

    // Fill K
    for(size_t i = 0; i < beams.size(); ++i){
        const std::vector<gp_Pnt>& beam = beams[i];
        const Force& f = this->data->forces[i/this->data->supports.size()];
        double A = f.S.get_area();
        double I = f.S.get_moment_of_inertia(1, 1);
        double dim = f.S.get_dimension();
        size_t sup_id = i%this->data->supports.size();

        // Get element attached to support
        {
            // Negative id indicates support, indexed at -(n+1) in this->data->supports
            gp_Vec v1 = gp_Vec(beam[0], beam[1]);
            this->nodes[0] = BeamNodeFactory::make_node(beam[0], -(1+sup_id), dim, v1, node_type);
            gp_Vec v2 = gp_Vec(beam[0], beam[1]) + gp_Vec(beam[1], beam[2]);
            this->nodes[1] = BeamNodeFactory::make_node(beam[1], cur_id, dim, v2, node_type);
            std::vector<long> pos(k_dim);

            // Calculate needed offsets to ensure K will not be singular
            if(k_dim >= 2 && this->data->supports[sup_id].X){
                pos[0] = -1;
            } else {
                pos[0] = cur_id + pos_offset;
                ++pos_offset;
            }
            if(k_dim >= 4 && this->data->supports[sup_id].Y){
                pos[1] = -1;
            } else {
                pos[1] = cur_id + pos_offset + 1;
                ++pos_offset;
            }
            if(k_dim >= 6 && this->data->supports[sup_id].MZ){
                pos[2] = -1;
            } else {
                pos[2] = cur_id + pos_offset + 2;
                ++pos_offset;
            }
            for(size_t j = k_dim/2; j < k_dim; ++j){
                pos[j] = cur_id + pos_offset + (j-k_dim/2);
            }
            double E = 0;
            if(data->type == utils::PROBLEM_TYPE_2D){
                math::Matrix R(
                    {v1.X(), -v1.Y(),
                     v1.Y(), v1.X()}, 2, 2);
                E = this->data->materials[0]->beam_E_2D(nullptr, gp_Pnt(), R);
            } else {
                logger::log_assert(false, logger::ERROR, "3D beam graph not fully implemented");
                // TODO
                //E = this->data->materials[0]->beam_E_3D(nullptr, gp_Pnt(), R);
            }
            for(size_t i = 0; i < 3; ++i){
                this->nodes[0]->u_pos[i] = pos[i];
                this->nodes[1]->u_pos[i] = pos[i+3];
            }
            elems[0] = BeamElementFactory::make_element(BeamElementFactory::BEAM_LINEAR_2D, this->nodes[0], this->nodes[1], I, A, E);
            this->insert_element_matrix(K, elems[0]->get_k(), pos, W, N);
        }

        // Get the rest of the elements
        for(size_t j = 1; j < beam.size()-1; ++j){
            size_t true_pos = cur_id+sup_id+1;
            gp_Vec v(beam[true_pos], beam[true_pos+1]);
            if(true_pos+2 != beam.size()){
                v += gp_Vec(beam[true_pos+1], beam[true_pos+2]);
            }
            this->nodes[true_pos+1] = BeamNodeFactory::make_node(beam[j+1], cur_id+1, dim, v, node_type);
            std::vector<long> pos(k_dim);
            size_t id0 = this->nodes[true_pos]->id;
            size_t id1 = this->nodes[true_pos+1]->id;
            for(size_t l = 0; l < k_dim/2; ++l){
                pos[l] = id0*k_dim/2+pos_offset+l;
            }
            for(size_t l = 0; l < k_dim/2; ++l){
                pos[l+k_dim/2] = id1*k_dim/2+pos_offset+l;
            }
            double E = 0;
            if(data->type == utils::PROBLEM_TYPE_2D){
                math::Matrix R(
                    {v.X(), -v.Y(),
                     v.Y(), v.X()}, 2, 2);
                E = this->data->materials[0]->beam_E_2D(nullptr, gp_Pnt(), R);
            } else {
                logger::log_assert(false, logger::ERROR, "3D beam graph not fully implemented");
                // TODO
                //E = this->data->materials[0]->beam_E_3D(nullptr, gp_Pnt(), R);
            }
            for(size_t i = 0; i < 3; ++i){
                this->nodes[true_pos+1]->u_pos[i] = pos[i+3];
            }
            elems[true_pos] = BeamElementFactory::make_element(BeamElementFactory::BEAM_LINEAR_2D, this->nodes[true_pos], this->nodes[true_pos+1], I, A, E);
            this->insert_element_matrix(K, elems[true_pos]->get_k(), pos, W, N);
            ++cur_id;
        }
        // Adds forces to force vector
        {
            if(k_dim >= 2){
                F[cur_id*k_dim/2 + pos_offset] = -f.vec.X();
            }
            if(k_dim >= 4){
                F[cur_id*k_dim/2 + pos_offset + 1] = -f.vec.Y();
            }
        }
    }

    // TODO NODE MERGING/CONNECTING ON INTERSECTION

    //  While LAPACK is naturally column major, and using row major may
    //  make this function slower due to the transposing, it is still better
    //  than what would be necessary to do in order to both accomodate
    //  column major and the possibility of expanding down the matrix band.

    int info = LAPACKE_dpbsv(LAPACK_ROW_MAJOR, 'L', W, N-1, 1, K.data(), W, F.data(), 1);
    logger::log_assert(info == 0, logger::ERROR, "LAPACKE returned {} while calculating displacements during sizing step.", info);
    K.clear();

    // Calculate reactions
    std::vector<double> reactions;
    reactions.reserve((elems.size()+1)*3);

    auto res = elems[0]->get_internal_loads(0, F);
    reactions.insert(reactions.end(), res.begin(), res.end());
    for(size_t i = 1; i < elems.size(); ++i){
        auto res = elems[i]->get_internal_loads(1, F);
        reactions.insert(reactions.end(), res.begin(), res.end());
    }
    // Delete elements
    while(!elems.empty()){
        delete elems.back();
        elems.pop_back();
    }

    return reactions;
}

void BeamGraph::insert_element_matrix(std::vector<double>& K, const std::vector<double>& k, const std::vector<long>& pos, int w, int& n) const{
    // Check if K needs to be resized
    long first = -1;
    long last = 0;
    for(size_t i = 0; i < pos.size(); ++i){
        if(first == -1 && pos[i] > -1){
            first = i;
        }
        if(pos[i] > -1){
            last = i;
        }
    }
    // Resize K if needed
    if(pos[last] - pos[first] > n){
        n = pos[last] - pos[first];
        K.resize(w*n);
    }
    
    // Add k to K
    // If pos[i] = -1, row and column are not to be inserted 
    for(long i = first; i < last+1; ++i){
        for(long j = i; j < last+1; ++j){
            if(pos[i] > -1 && pos[j] > -1){
                K[utils::to_band(pos[j], pos[i], w)] += k[utils::to_upper_triangular(i,j)];
            }
        }
    }
}

BeamGraph::~BeamGraph(){
    this->clear_nodes();
}

void BeamGraph::clear_nodes(){
    while(!this->nodes.empty()){
        delete this->nodes.back();
        this->nodes.pop_back();
    }
}

