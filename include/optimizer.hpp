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

#ifndef OPTIMIZER_HPP
#define OPTIMIZER_HPP

#include <BRepBuilderAPI_Copy.hxx>
#include <TopoDS_Shape.hxx>
#include "finite_element.hpp"
#include "visualization.hpp"

class Optimizer{
    public:
    virtual ~Optimizer() = default;

    virtual void initialize_views(Visualization* viz) = 0;
    virtual TopoDS_Shape optimize(FiniteElement* fem, Meshing* mesh) = 0;

    protected:

    inline TopoDS_Shape make_shape(const std::vector<double>& x, const std::vector<Geometry*> geometries, const double result_threshold) const{
        logger::quick_log(" "); 
        logger::quick_log("Saving resulting geometries...");
        // TODO: make it work with multiple geometries
        std::cout << "\r" << 0 << "%         ";
        TopoDS_Shape result = BRepBuilderAPI_Copy(geometries[0]->shape);
        for(size_t i = 0; i < x.size(); ++i){
            if(x[i] >= result_threshold){
                result = utils::cut_shape(result, geometries[0]->mesh[i]->get_shape());
            }
            double pc = i/(double)(x.size()-1);
            std::cout << "\r" << pc*100 << "%         ";
        }
        result = utils::cut_shape(geometries[0]->shape, result);
        std::cout << "\r" << 100 << "%         ";
        logger::quick_log(" "); 
        return result;
    }
};


#endif
