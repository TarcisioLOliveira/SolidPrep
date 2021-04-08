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

#include "utils.hpp"
#include "logger.hpp"
#include "Interface_Static.hxx"

namespace utils{

void shape_to_file(const std::string& s, const TopoDS_Shape& t){
    STEPControl_Writer writer;
    STEPControl_StepModelType mode = STEPControl_AsIs;
    Interface_Static::SetIVal("write.surfacecurve.mode",0);
    IFSelect_ReturnStatus stat = writer.Transfer(t,mode);
    if(stat != IFSelect_RetDone){
        writer.PrintStatsTransfer(4, 0);
        std::cout << "ERROR: failed to translate model to STEP. Filename: " << s << std::endl;
        exit(EXIT_FAILURE);
    }
    IFSelect_ReturnStatus stat2 = writer.Write(s.c_str());
    if(stat2 != IFSelect_RetDone){
        writer.PrintStatsTransfer(4, 0);
        std::cout << "ERROR: failed to write translated model. Filename: " << s << std::endl;
        exit(EXIT_FAILURE);
    }
    (void) stat;
    (void) stat2;
}

}
