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

#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <utility>

namespace utils{

    /**
     * Formats string, replacing each "{}" with one of the arguments. Fails
     * silently if number of "{}" is different from the number of arguments.
     * Overload for lack of arguments, just returns the string.
     *
     * @param str String to be formatted.
     *
     * @return Formatted string.
     */
    static std::string format(std::string str){
        return str;
    }

    /**
     * Formats string, replacing each "{}" with one of the arguments. Fails
     * silently if number of "{}" is different from the number of arguments.
     * Overload for single argument.
     *
     * @param str String to be formatted.
     * @param a Argument to be added to string.
     *
     * @return Formatted string.
     */
    template <typename Arg1>
    static std::string format(std::string str, Arg1&& a){
        size_t pos = str.find("{}");

        if(pos == std::string::npos){
            return str;
        } else {
            return str.substr(0, pos) + std::to_string(a) + str.substr(pos+2);
        }
    }

    template<>
    std::string format<std::string>(std::string str, std::string&& a){
        size_t pos = str.find("{}");

        if(pos == std::string::npos){
            return str;
        } else {
            return str.substr(0, pos) + a + str.substr(pos+2);
        }
    }

    template<>
    std::string format<std::string&>(std::string str, std::string& a){
        size_t pos = str.find("{}");

        if(pos == std::string::npos){
            return str;
        } else {
            return str.substr(0, pos) + a + str.substr(pos+2);
        }
    }

    template<>
    std::string format<const char*>(std::string str, const char*&& a){
        size_t pos = str.find("{}");

        if(pos == std::string::npos){
            return str;
        } else {
            return str.substr(0, pos) + std::string(a) + str.substr(pos+2);
        }
    }

    /**
     * Formats string, replacing each "{}" with one of the arguments. Fails
     * silently if number of "{}" is different from the number of arguments.
     *
     * @param str String to be formatted.
     * @param a1 Argument to be added to string.
     * @param args Other arguments to be added to string.
     *
     * @return Formatted string.
     */
    template <typename Arg1, typename ... Args>
    static std::string format(std::string str, Arg1&& a1, Args&& ... args){
        return format(format(str, std::forward<Arg1>(a1)), std::forward<Args>(args) ...);
    }

}

#endif