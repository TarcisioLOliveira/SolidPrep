/*
 *   Copyright (C) 2024 Tarcísio Ladeia de Oliveira.
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

#include "catch2/catch.hpp"
#include "logger.hpp"
#include "math/matrix.hpp"

TEST_CASE( "Matrix Transpose: check individual values" ) {
    math::Matrix m({0, 1,
                    2, 3,
                    4, 5}, 3, 2);
    REQUIRE( m(0,0) == 0 );
    REQUIRE( m(0,1) == 1 );
    REQUIRE( m(1,0) == 2 );
    REQUIRE( m(1,1) == 3 );
    REQUIRE( m(2,0) == 4 );
    REQUIRE( m(2,1) == 5 );

    REQUIRE( m.T()(0,0) == 0 );
    REQUIRE( m.T()(1,0) == 1 );
    REQUIRE( m.T()(0,1) == 2 );
    REQUIRE( m.T()(1,1) == 3 );
    REQUIRE( m.T()(0,2) == 4 );
    REQUIRE( m.T()(1,2) == 5 );
}
TEST_CASE( "Matrix Transpose: equalities" ) {
    math::Matrix m({0, 1,
                    2, 3,
                    4, 5}, 3, 2);
    math::Matrix m2 = m.T();
    math::Matrix m3 = m2.T();
    math::Matrix mT({0, 2, 4,
                     1, 3, 5},
                     2, 3);
    math::Matrix m4({1, 2,
                     3, 4},
                     2, 2); 
    math::Matrix m5({1, 2,
                     2, 4},
                     2, 2); 
    REQUIRE( m.T() == m2 );
    REQUIRE( m2 == m.T() );
    REQUIRE( m2.T() == m3 );
    REQUIRE( m3 == m2.T() );
    REQUIRE( m3 == m );
    REQUIRE( m == m3 );
    REQUIRE( m != m.T() );
    REQUIRE( m != m2 );
    REQUIRE( m2 != m3 );
    REQUIRE( m.T() != m3 );
    REQUIRE( m.T() == mT );
    REQUIRE( mT == m.T() );
    REQUIRE( m4 != m4.T() );
    REQUIRE( m5 == m5.T() );
}
TEST_CASE( "Matrix Transpose: operations" ) {
    math::Matrix m1({0, 1,
                     2, 3,
                     4, 5}, 3, 2);
    math::Matrix m2({0, 1, 2,
                     3, 4, 5}, 2, 3);
    math::Matrix m3({1, 2,
                     3, 4},
                     2, 2); 
    math::Matrix m4({1, 2, 3,
                     4, 5, 6,
                     7, 8, 9},
                     3, 3); 
    REQUIRE( m1 + m2.T() ==
            math::Matrix({
            0, 4,
            3, 7,
            6, 10,
            }, 3, 2));
    REQUIRE( m1.T() - m2 ==
            math::Matrix({
            0, 1, 2,
            -2, -1, 0,
            }, 2, 3));
    REQUIRE( 5*m1.T() + m2 ==
            math::Matrix({
            0, 11, 22,
            8, 19, 30,
            }, 2, 3));
    REQUIRE( 5*(m1.T() - m2) ==
            math::Matrix({
            0, 5, 10,
            -10, -5, 0,
            }, 2, 3));
    REQUIRE( m1*m2 ==
            math::Matrix({
            3, 4, 5,
            9, 14, 19,
            15, 24, 33,
            }, 3, 3));
    REQUIRE( m1*m1.T() ==
            math::Matrix({
            1, 3, 5,
            3, 13, 23,
            5, 23, 41,
            }, 3, 3));
    REQUIRE( m2*m2.T() ==
            math::Matrix({
            5, 14,
            14, 50,
            }, 2, 2));
    REQUIRE( m1.T()*m1 ==
            math::Matrix({
            20, 26,
            26, 35,
            }, 2, 2));
    REQUIRE( m2.T()*m2 ==
            math::Matrix({
            9, 12, 15,
            12, 17, 22,
            15, 22, 29,
            }, 3, 3));
    REQUIRE( m1*m3*m1.T() ==
            math::Matrix({
            4, 18, 32,
            16, 70, 124,
            28, 122, 216,
            }, 3, 3));
    REQUIRE( m1*m3.T()*m1.T() ==
            math::Matrix({
            4, 16, 28,
            18, 70, 122,
            32, 124, 216,
            }, 3, 3));
    REQUIRE( m2.T()*m3*m2 ==
            math::Matrix({
            36, 57, 78,
            54, 85, 116,
            72, 113, 154,
            }, 3, 3));
    REQUIRE( m2.T()*m3.T()*m2 ==
            math::Matrix({
            36, 54, 72,
            57, 85, 113,
            78, 116, 154,
            }, 3, 3));
    REQUIRE( m2*m4*m2.T() ==
            math::Matrix({
            69, 258,
            222, 816,
            }, 2, 2));
    REQUIRE( m2*m4.T()*m2.T() ==
            math::Matrix({
            69, 222,
            258, 816,
            }, 2, 2));
    REQUIRE( m1.T()*m4*m1 ==
            math::Matrix({
            276, 402,
            378, 549,
            }, 2, 2));
    REQUIRE( m1.T()*m4.T()*m1 ==
            math::Matrix({
            276, 378,
            402, 549,
            }, 2, 2));
}
