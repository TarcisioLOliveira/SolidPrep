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
#include "math/matrix.hpp"

inline bool equals(double a, double b, double eps = 1e-7){
    return std::abs(a - b) < eps;
}

TEST_CASE( "Vector: check individual values (initializer list)" ) {
    math::Vector v({0, 1, 2, 3, 4});
    REQUIRE( v[0] == 0 );
    REQUIRE( v[1] == 1 );
    REQUIRE( v[2] == 2 );
    REQUIRE( v[3] == 3 );
    REQUIRE( v[4] == 4 );
}
TEST_CASE( "Vector: check individual values (vector)" ) {
    std::vector<double> V{0, 1, 2, 3, 4};
    math::Vector v(V);
    REQUIRE( v[0] == 0 );
    REQUIRE( v[1] == 1 );
    REQUIRE( v[2] == 2 );
    REQUIRE( v[3] == 3 );
    REQUIRE( v[4] == 4 );
}
TEST_CASE( "Vector: check individual values (array)" ) {
    double V[]{0, 1, 2, 3, 4};
    math::Vector v(V, 5);
    REQUIRE( v[0] == 0 );
    REQUIRE( v[1] == 1 );
    REQUIRE( v[2] == 2 );
    REQUIRE( v[3] == 3 );
    REQUIRE( v[4] == 4 );
}

TEST_CASE( "Vector: check individual values (transpose)" ){
    math::Vector v({0, 1, 2, 3, 4});
    REQUIRE( v.T()[0] == 0 );
    REQUIRE( v.T()[1] == 1 );
    REQUIRE( v.T()[2] == 2 );
    REQUIRE( v.T()[3] == 3 );
    REQUIRE( v.T()[4] == 4 );

    auto vT = v.T();
    REQUIRE( vT[0] == 0 );
    REQUIRE( vT[1] == 1 );
    REQUIRE( vT[2] == 2 );
    REQUIRE( vT[3] == 3 );
    REQUIRE( vT[4] == 4 );

    auto vT2 = 1*v.T();
    REQUIRE( vT2[0] == 0 );
    REQUIRE( vT2[1] == 1 );
    REQUIRE( vT2[2] == 2 );
    REQUIRE( vT2[3] == 3 );
    REQUIRE( vT2[4] == 4 );

    auto vTT = vT.T();
    REQUIRE( vTT[0] == 0 );
    REQUIRE( vTT[1] == 1 );
    REQUIRE( vTT[2] == 2 );
    REQUIRE( vTT[3] == 3 );
    REQUIRE( vTT[4] == 4 );
}

TEST_CASE( "Vector: check individual values (fill)" ) {
    math::Vector v(5, 10);
    REQUIRE( v[0] == 10 );
    REQUIRE( v[1] == 10 );
    REQUIRE( v[2] == 10 );
    REQUIRE( v[3] == 10 );
    REQUIRE( v[4] == 10 );
}

TEST_CASE( "Vector: equalities" ) {
    math::Vector v1({1, 2, 3, 4});
    math::Vector v2({1, 2, 3, 4});
    math::Vector v3({1, 2, 3, 4, 5});
    math::Vector v4({5, 6, 7, 8});
    math::Vector v5({1, 2, 3, 5});

    REQUIRE(v1 == v2);
    REQUIRE(v2 == v1);
    REQUIRE(v1 == v2.T());
    REQUIRE(v2 == v1.T());
    REQUIRE(v1.T() == v2);
    REQUIRE(v2.T() == v1);
    REQUIRE(v1.T() == v2.T());
    REQUIRE(v2.T() == v1.T());

    REQUIRE(v1.is_equal(v2));
    REQUIRE(v2.is_equal(v1));
    REQUIRE(v1.is_equal(v2.T()));
    REQUIRE(v2.is_equal(v1.T()));
    REQUIRE(v1.T().is_equal(v2));
    REQUIRE(v2.T().is_equal(v1));
    REQUIRE(v1.T().is_equal(v2.T()));
    REQUIRE(v2.T().is_equal(v1.T()));

    REQUIRE(v1 != v3);
    REQUIRE(v2 != v3);
    REQUIRE(v3 != v1);
    REQUIRE(v3 != v2);
    REQUIRE(v1.T() != v3);
    REQUIRE(v2.T() != v3);
    REQUIRE(v3.T() != v1);
    REQUIRE(v3.T() != v2);
    REQUIRE(v1 != v3.T());
    REQUIRE(v2 != v3.T());
    REQUIRE(v3 != v1.T());
    REQUIRE(v3 != v2.T());
    REQUIRE(v1.T() != v3.T());
    REQUIRE(v2.T() != v3.T());
    REQUIRE(v3.T() != v1.T());
    REQUIRE(v3.T() != v2.T());

    REQUIRE(v1 != v4);
    REQUIRE(v2 != v4);
    REQUIRE(v4 != v1);
    REQUIRE(v4 != v2);
    REQUIRE(v1.T() != v4);
    REQUIRE(v2.T() != v4);
    REQUIRE(v4.T() != v1);
    REQUIRE(v4.T() != v2);
    REQUIRE(v1 != v4.T());
    REQUIRE(v2 != v4.T());
    REQUIRE(v4 != v1.T());
    REQUIRE(v4 != v2.T());
    REQUIRE(v1.T() != v4.T());
    REQUIRE(v2.T() != v4.T());
    REQUIRE(v4.T() != v1.T());
    REQUIRE(v4.T() != v2.T());

    REQUIRE(v1 != v5);
    REQUIRE(v2 != v5);
    REQUIRE(v5 != v1);
    REQUIRE(v5 != v2);
    REQUIRE(v1.T() != v5);
    REQUIRE(v2.T() != v5);
    REQUIRE(v5.T() != v1);
    REQUIRE(v5.T() != v2);
    REQUIRE(v1 != v5.T());
    REQUIRE(v2 != v5.T());
    REQUIRE(v5 != v1.T());
    REQUIRE(v5 != v2.T());
    REQUIRE(v1.T() != v5.T());
    REQUIRE(v2.T() != v5.T());
    REQUIRE(v5.T() != v1.T());
    REQUIRE(v5.T() != v2.T());
}

TEST_CASE( "Vector: operations" ) {
    auto m1 = math::Vector({1, 2, 3});
    auto m2 = math::Vector({10, 9, 5});
    auto m3 = math::Matrix({
        1, 2, 3,
        4, 5, 6,
        7, 8, 9,
        }, 3, 3);
    auto m4 = math::Vector({4, -6, 7, -9});
    auto m5 = math::Matrix({
        4, 6, 7, 1,
        6, 8, -1, 0,
        -3, 4, 5, -6,
        }, 3, 4);
    auto m6 = math::Vector({8, 4, 1});
    auto m7 = math::Vector({1, 5});
    auto m8 = math::Vector({8, 4});
    auto m9 = math::Vector({-1, 7});

    REQUIRE( m1 + m2 == math::Vector({11, 11, 8}) );
    REQUIRE( m1 - m2 == math::Vector({-9, -7, -2}) );
    REQUIRE( m1 + 5*m2 == math::Vector({51, 47, 28}) );
    REQUIRE( m1 - 10*m2 == math::Vector({-99, -88, -47}) );
    REQUIRE( 6*m1 + m2 == math::Vector({16, 21, 23}) );
    REQUIRE( 9*m1 - m2 == math::Vector({-1, 9, 22}) );
    REQUIRE( 8*(m1 + m2) == math::Vector({88, 88, 64}) );
    REQUIRE( 5*(m1 - m2) == math::Vector({-45, -35, -10}) );
    REQUIRE( m1 + m2/6 == math::Vector({2.666666666666667, 3.5, 3.8333333333333335}) );
    REQUIRE( m1 - m2/8 == math::Vector({-0.25, 0.875, 2.375}) );
    REQUIRE( m1/2 + m2 == math::Vector({10.5, 10.0, 6.5}) );
    REQUIRE( m1/7 - m2 == math::Vector({-9.857142857142858, -8.714285714285714, -4.571428571428571}) );
    REQUIRE( (m1 + m2)/8 == math::Vector({1.375, 1.375, 1.0}) );
    REQUIRE( (m1 - m2)/10 == math::Vector({-0.9, -0.7, -0.2}) );
    REQUIRE( m1.T() * m2 == 43 );
    REQUIRE( m1.T() * (5*m2) == 215 );
    REQUIRE( (6*m1).T() * m2 == 258 );
    REQUIRE( (6*m1.T()) * m2 == 258 );
    REQUIRE( m1 * m2.T() == math::Matrix({
    10, 9, 5,
    20, 18, 10,
    30, 27, 15,
    }, 3, 3));
    REQUIRE( m1 * (7*m2.T()) == math::Matrix({
    70, 63, 35,
    140, 126, 70,
    210, 189, 105,
    }, 3, 3));
    REQUIRE( m1/8 * m2.T() == math::Matrix({
    1.25, 1.125, 0.625,
    2.5, 2.25, 1.25,
    3.75, 3.375, 1.875,
    }, 3, 3));
    REQUIRE( (m1 * m2.T())*7 == math::Matrix({
    70, 63, 35,
    140, 126, 70,
    210, 189, 105,
    }, 3, 3));
    REQUIRE( m3 * m2 == math::Vector({43, 115, 187}) );
    REQUIRE( m3 * m1 == math::Vector({14, 32, 50}) );
    REQUIRE( m2.T() * m3 == math::Vector({81, 105, 129}) );
    REQUIRE( m1.T() * m3 == math::Vector({30, 36, 42}) );
    REQUIRE( m1.T() * m3 * m2 == 834 );
    REQUIRE( m3.T() * m2 == math::Vector({81, 105, 129}) );
    REQUIRE( m3.T() * m1 == math::Vector({30, 36, 42}) );
    REQUIRE( m2.T() * m3.T() == math::Vector({43, 115, 187}) );
    REQUIRE( m1.T() * m3.T() == math::Vector({14, 32, 50}) );
    REQUIRE( m1.T() * m3.T() * m2 == 678 );
    REQUIRE( m1.T() * m3 * m3.T() * (m1 + m2) == 15588 );
    REQUIRE( (m1 + m2).T() * m3 * m3.T() * m2 == 45855 );
    REQUIRE( m1.T() * m3 * m3.T() * (m3 * m1) == 64656 );
    REQUIRE( m1.T() * m3 * m3.T() * (m3.T() * m1) == 63504 );
    REQUIRE( m1.T() * m3 * m3.T() * (m2*7 + m3.T() * m1) == 144900 );
    REQUIRE( (m1 + 3*m2).T() * (5*m3) * m3.T() * (m2*7 + m3.T() * m1) == 7117695 );

    REQUIRE( m5 * m4 == math::Vector({20, -31, 53}) );
    REQUIRE( (m5 * m4).T() == math::Vector({20, -31, 53}) );
    REQUIRE( m4.T() * m5.T() == math::Vector({20, -31, 53}) );
    REQUIRE( m2.T() * m5 * m4 == 186 );
    REQUIRE( m4.T() * m5.T() * m2 == 186 );
    REQUIRE( m4.T() * m5.T() * m3 * m5 * m4 == 14364 );
    REQUIRE( m1.T() * m5 * m5.T() * m2 == 7781 );

    REQUIRE( m1.norm() == 3.7416573867739413 );
    REQUIRE( m2.norm() == 14.352700094407323 );
    REQUIRE( m4.norm() == 13.490737563232042 );
    REQUIRE( m6.norm() == 9.0 );
    REQUIRE( m7.norm() == 5.0990195135927845 );
    REQUIRE( m8.norm() == 8.94427190999916 );
    REQUIRE( m9.norm() == 7.0710678118654755 );

    REQUIRE( m1.normalized() == math::Vector({0.2672612419124244, 0.5345224838248488, 0.8017837257372732}) );
    REQUIRE( m2.normalized() == math::Vector({0.6967330142916177, 0.6270597128624559, 0.34836650714580886}) );
    REQUIRE( m4.normalized() == math::Vector({0.29649972666444047, -0.4447495899966607, 0.5188745216627708, -0.667124384994991}) );
    REQUIRE( m6.normalized() == math::Vector({0.8888888888888888, 0.4444444444444444, 0.1111111111111111}) );
    REQUIRE( m7.normalized() == math::Vector({0.19611613513818404, 0.9805806756909202}) );
    REQUIRE( m8.normalized() == math::Vector({0.8944271909999159, 0.4472135954999579}) );
    REQUIRE( m9.normalized() == math::Vector({-0.1414213562373095, 0.9899494936611665}) );

    REQUIRE( m1.cross(m2) == math::Vector({-17, 25, -11}) );
    REQUIRE( m2.cross(m1) == math::Vector({17, -25, 11}) );
    REQUIRE( m1.cross(m6) == math::Vector({-10, 23, -12}) );
    REQUIRE( m6.cross(m1) == math::Vector({10, -23, 12}) );
    REQUIRE( m2.cross(m6) == math::Vector({-11, 30, -32}) );
    REQUIRE( m6.cross(m2) == math::Vector({11, -30, 32}) );

    REQUIRE( m7.cross(m8) == math::Vector({0, 0, -36}) );
    REQUIRE( m8.cross(m7) == math::Vector({0, 0, 36}) );
    REQUIRE( m7.cross(m9) == math::Vector({0, 0, 12}) );
    REQUIRE( m9.cross(m7) == math::Vector({0, 0, -12}) );
    REQUIRE( m8.cross(m9) == math::Vector({0, 0, 60}) );
    REQUIRE( m9.cross(m8) == math::Vector({0, 0, -60}) );
}

TEST_CASE( "LU" ) {
    auto v1 = math::Vector({10, 9, 5});
    math::Matrix m1({2, 4, 8,
                     4, 7, 6,
                     1, 7, 8},
                     3, 3);
    math::Matrix m2(m1);
    math::Matrix m3({
        1, 2, 3,
        4, 5, 6,
        7, 8, 9,
        }, 3, 3);

    auto v2 = v1;
    auto m4 = m3;

    math::LU lu1(m1);
    REQUIRE( m2.determinant() != 0 );
    REQUIRE( equals(lu1.determinant(), m2.determinant()) );
    lu1.solve(v1);
    lu1.solve(m2);
    lu1.solve(m3);
    lu1.solve_transposed(v2);
    lu1.solve_transposed(m4);
    REQUIRE( v1.is_equal(math::Vector({2.13043478, -0.95652174, 1.19565217})) );
    REQUIRE( m2 == math::Matrix::identity(3) );
    REQUIRE( m3.is_equal(
        math::Matrix(
        {-1.23913043, -1.17391304, -1.10869565,
          1.58695652,  1.60869565,  1.63043478,
         -0.35869565, -0.26086957, -0.16304348}, 3, 3)) );
    REQUIRE( v2.is_equal(math::Vector({0.11956522, 2.84782609, -1.63043478})) );
    REQUIRE( m4.is_equal(
            math::Matrix(
        {0.61956522,  0.7173913 ,  0.81521739,
        -0.15217391,  0.08695652,  0.32608696,
         0.36956522,  0.2173913 ,  0.06521739}, 3, 3)) );

}

TEST_CASE( "LU copy" ) {
    auto v1 = math::Vector({10, 9, 5});
    math::Matrix m1({2, 4, 8,
                     4, 7, 6,
                     1, 7, 8},
                     3, 3);
    math::Matrix m2(m1);
    math::Matrix m3({
        1, 2, 3,
        4, 5, 6,
        7, 8, 9,
        }, 3, 3);

    auto v2 = v1;
    auto m4 = m3;

    math::LU lu1(m1, true);
    REQUIRE( m2.determinant() != 0 );
    REQUIRE( equals(lu1.determinant(), m2.determinant()) );
    lu1.solve(v1);
    lu1.solve(m2);
    lu1.solve(m3);
    lu1.solve_transposed(v2);
    lu1.solve_transposed(m4);
    REQUIRE( v1.is_equal(math::Vector({2.13043478, -0.95652174, 1.19565217})) );
    REQUIRE( m2 == math::Matrix::identity(3) );
    REQUIRE( m3.is_equal(
        math::Matrix(
        {-1.23913043, -1.17391304, -1.10869565,
          1.58695652,  1.60869565,  1.63043478,
         -0.35869565, -0.26086957, -0.16304348}, 3, 3)) );
    REQUIRE( v2.is_equal(math::Vector({0.11956522, 2.84782609, -1.63043478})) );
    REQUIRE( m4.is_equal(
            math::Matrix(
        {0.61956522,  0.7173913 ,  0.81521739,
        -0.15217391,  0.08695652,  0.32608696,
         0.36956522,  0.2173913 ,  0.06521739}, 3, 3)) );

}

TEST_CASE( "Cholesky" ) {
    auto v1 = math::Vector({10, 9, 5});
    math::Matrix m1({10, 3, 1,
                     3, 10, 3,
                     1, 3, 10},
                     3, 3);
    math::Matrix m2(m1);
    math::Matrix m3({
        5, -1, 0,
        -1, 10, -1,
        0, -1, 5,
        }, 3, 3);

    math::Cholesky ch1(m1);
    REQUIRE( m2.determinant() != 0 );
    REQUIRE( equals(ch1.determinant(), m2.determinant()) );
    ch1.solve(v1);
    ch1.solve(m2);
    ch1.solve(m3);
    REQUIRE( v1.is_equal(math::Vector({0.79951691, 0.58695652, 0.24396135})) );
    REQUIRE( m2 == math::Matrix::identity(3) );
    REQUIRE( m3.is_equal(
        math::Matrix(
        { 0.5821256 , -0.43478261,  0.02657005,
         -0.2826087 ,  1.26086957, -0.2826087 ,
          0.02657005, -0.43478261,  0.5821256 }, 3, 3)) );
}

TEST_CASE( "Cholesky copy" ) {
    auto v1 = math::Vector({10, 9, 5});
    math::Matrix m1({10, 3, 1,
                     3, 10, 3,
                     1, 3, 10},
                     3, 3);
    math::Matrix m2(m1);
    math::Matrix m3({
        5, -1, 0,
        -1, 10, -1,
        0, -1, 5,
        }, 3, 3);

    math::Cholesky ch1(m1, true);
    REQUIRE( m2.determinant() != 0 );
    REQUIRE( equals(ch1.determinant(), m2.determinant()) );
    ch1.solve(v1);
    ch1.solve(m2);
    ch1.solve(m3);
    REQUIRE( v1.is_equal(math::Vector({0.79951691, 0.58695652, 0.24396135})) );
    REQUIRE( m2 == math::Matrix::identity(3) );
    REQUIRE( m3.is_equal(
        math::Matrix(
        { 0.5821256 , -0.43478261,  0.02657005,
         -0.2826087 ,  1.26086957, -0.2826087 ,
          0.02657005, -0.43478261,  0.5821256 }, 3, 3)) );
}
