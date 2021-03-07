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

#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string>

namespace logger{

    enum Type{
        WARNING,
        ERROR
    };

    template<typename ... Args>
    static void log_assert(bool expr, Type t, std::string message, Args ... args){
        if(!expr){
            if(t == WARNING){
                message = "WARNING: "+message;
            } else if(t == ERROR){
                message = "ERROR: "+message;
            }
            fprintf(stderr, message.c_str(), args ...);
            if(t == ERROR){
                exit(EXIT_FAILURE);
            }
        }
    }
}

#endif
