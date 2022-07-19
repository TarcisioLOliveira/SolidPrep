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

#include "geometry.hpp"
#include "STEPCAFControl_Reader.hxx"
#include "BRepClass3d_SolidClassifier.hxx"
#include "BRepBuilderAPI_Transform.hxx"

Geometry::Geometry(const std::string& path, double scale, utils::ProblemType type,
        MeshElementFactory* elem_type, bool do_topopt, Material* material, std::vector<Material*> alt_materials ):
    shape(this->load_shape(path, scale)), scale(scale), material(material), alternate_materials(alt_materials.begin(), alt_materials.end()),
    element_type(elem_type), do_topopt(do_topopt), mesh(), 
    type(type){}

TopoDS_Shape Geometry::load_shape(const std::string& path, double scale) const{
    TopoDS_Shape s;
    
    STEPControl_Reader reader;
    IFSelect_ReturnStatus stat = reader.ReadFile(path.c_str());
    if(stat != IFSelect_RetDone){
        reader.PrintCheckLoad(false, IFSelect_ItemsByEntity);
        exit(EXIT_FAILURE);
    }

    Standard_Integer NbRoots = reader.NbRootsForTransfer();
    Standard_Integer num = reader.TransferRoots();
    (void) NbRoots;
    (void) num;
    s = reader.OneShape();
    if(s.IsNull()){
        reader.PrintCheckTransfer(true, IFSelect_ItemsByEntity);
        exit(EXIT_FAILURE);
    }

    if(scale != 1){
        gp_Trsf t;
        t.SetScale(gp_Pnt(0, 0, 0), scale);
        BRepBuilderAPI_Transform transf(s, t, true);
        s = transf.Shape();
    }

    return s;
}

bool Geometry::is_inside(const gp_Pnt& p) const{
    if(this->type == utils::PROBLEM_TYPE_2D){
        return this->is_inside_2D(p);
    } else if(this->type == utils::PROBLEM_TYPE_3D){
        return this->is_inside_3D(p);
    }
    return false;
}

bool Geometry::is_inside_2D(const gp_Pnt& p) const{
    BRepClass3d_SolidClassifier insider(this->shape, p, 0.01);
    return insider.State() == TopAbs_ON;
}

bool Geometry::is_inside_3D(const gp_Pnt& p) const{
    BRepClass3d_SolidClassifier insider(this->shape, p, 0.01);
    return insider.State() == TopAbs_IN;
}
