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

TEST_CASE( "Matrix: check individual values (initializer list)" ) {
    math::Matrix m({0, 1, 2,
                    3, 4, 5,
                    6, 7, 8},
                    3, 3);
    REQUIRE( m(0,0) == 0 );
    REQUIRE( m(0,1) == 1 );
    REQUIRE( m(0,2) == 2 );
    REQUIRE( m(1,0) == 3 );
    REQUIRE( m(1,1) == 4 );
    REQUIRE( m(1,2) == 5 );
    REQUIRE( m(2,0) == 6 );
    REQUIRE( m(2,1) == 7 );
    REQUIRE( m(2,2) == 8 );
}
TEST_CASE( "Matrix: check individual values (vector)" ) {
    std::vector<double> mv{0, 1, 2,
                           3, 4, 5,
                           6, 7, 8};
    math::Matrix m(mv,
                    3, 3);
    REQUIRE( m(0,0) == 0 );
    REQUIRE( m(0,1) == 1 );
    REQUIRE( m(0,2) == 2 );
    REQUIRE( m(1,0) == 3 );
    REQUIRE( m(1,1) == 4 );
    REQUIRE( m(1,2) == 5 );
    REQUIRE( m(2,0) == 6 );
    REQUIRE( m(2,1) == 7 );
    REQUIRE( m(2,2) == 8 );
}
TEST_CASE( "Matrix: check individual values (array)" ) {
    double mv[9] = {0, 1, 2,
                    3, 4, 5,
                    6, 7, 8};
    math::Matrix m(mv,
                    3, 3);
    REQUIRE( m(0,0) == 0 );
    REQUIRE( m(0,1) == 1 );
    REQUIRE( m(0,2) == 2 );
    REQUIRE( m(1,0) == 3 );
    REQUIRE( m(1,1) == 4 );
    REQUIRE( m(1,2) == 5 );
    REQUIRE( m(2,0) == 6 );
    REQUIRE( m(2,1) == 7 );
    REQUIRE( m(2,2) == 8 );
}
TEST_CASE( "Matrix: check individual values (fill)" ) {
    math::Matrix m(3, 3, 10);
    REQUIRE( m(0,0) == 10 );
    REQUIRE( m(0,1) == 10 );
    REQUIRE( m(0,2) == 10 );
    REQUIRE( m(1,0) == 10 );
    REQUIRE( m(1,1) == 10 );
    REQUIRE( m(1,2) == 10 );
    REQUIRE( m(2,0) == 10 );
    REQUIRE( m(2,1) == 10 );
    REQUIRE( m(2,2) == 10 );
}
TEST_CASE( "Matrix: equalities" ) {
    math::Matrix m1({0, 1, 2,
                     3, 4, 5,
                     6, 7, 8},
                     3, 3);
    math::Matrix m2({0, 1, 2,
                     3, 4, 5,
                     6, 7, 8},
                     3, 3);
    REQUIRE( m1 == m2 );
    m2 = m1;
    REQUIRE( m1 == m2 );
    m2(1,1) = 10;
    REQUIRE( m2(1,1) == 10 );
    REQUIRE( m1 != m2 );
    m1 = math::Matrix(3,4);
    REQUIRE( m1 != m2 );
    m2 = math::Matrix(3,4);
    REQUIRE( m1 == m2 );
    m2(1,1) = 1;
    REQUIRE( m1 != m2 );

    REQUIRE( m1.get_W() == m2.get_W() );
    REQUIRE( m1.get_H() == m2.get_H() );
    REQUIRE( m1.get_W() != m2.get_H() );
    REQUIRE( m1.get_H() != m2.get_W() );

    m1 = 
      math::Matrix({0, 1, 2,
                    3, 4, 5,
                    6, 7, 8},
                    3, 3);
    m2 = 
      math::Matrix({0, 1, 2,
                    3, 4, 5,
                    6, 7, 8},
                    3, 3);
    math::Matrix m3(m1);
    REQUIRE( m1 == m3 );
    {
        math::Matrix m4(m1);
    }
    REQUIRE( m1 == m2 );
    math::Matrix m4(
      math::Matrix({0, 1, 2,
                    3, 4, 5,
                    6, 7, 8},
                    3, 3)
    );
    REQUIRE( m1 == m4 );
    m4 = 
      math::Matrix({0, 1, 2,
                    3, 4, 5,
                    6, 7, 8},
                    3, 3);
    REQUIRE( m1 == m4 );
}
TEST_CASE( "Matrix: identity" ) {
    auto I = math::Matrix::identity(3);
    math::Matrix m({1, 0, 0,
                    0, 1, 0,
                    0, 0, 1},
                    3, 3);
    math::Matrix m1({0, 1, 2,
                     3, 4, 5,
                     6, 7, 8},
                     3, 3);
    REQUIRE( I == m );
    REQUIRE( I*m == I );
    REQUIRE( I*m1 == m1 );
}
TEST_CASE( "Matrix: determinant" ) {
    math::Matrix m1({0, 1, 2,
                     3, 4, 5,
                     6, 7, 8},
                     3, 3);
    math::Matrix m2({2, 4, 8,
                     4, 7, 6,
                     1, 7, 8},
                     3, 3);
    math::Matrix m3({2, 4, 8,
                     1, 7, 8,
                     4, 7, 6},
                     3, 3);
    REQUIRE( m1.determinant() ==  0 );
    REQUIRE( m2.determinant() ==  92.0 );
    REQUIRE( m3.determinant() == -92.0 );
}
TEST_CASE( "Matrix: inverse" ) {
    auto I = math::Matrix::identity(3);
    math::Matrix m2({2, 4, 8,
                     4, 7, 6,
                     1, 7, 8},
                     3, 3);
    auto m3 = m2.get_inverted_LU();
    auto m4 = m3*m2;
    math::Matrix minv({
             0.15217391304347826081,  0.26086956521739130436 , -0.34782608695652173909,
            -0.28260869565217391305,  0.086956521739130434773,  0.21739130434782608695,
             0.22826086956521739132, -0.10869565217391304348 , -0.021739130434782608698,
        }, 3, 3);
    REQUIRE( m4.is_equal(I) );
    REQUIRE( (m3*m2).is_equal(I) );
    REQUIRE( (m2*m3).is_equal(I) );
    REQUIRE( equals(m2.determinant(), 1.0/m3.determinant()) );
    REQUIRE( m3.is_equal(minv) );

    m2 = math::Matrix({10, 3, 1,
                       3, 10, 3,
                       1, 3, 10},
                       3, 3);
    m3 = m2.get_inverted_cholesky();
    m4 = m3*m2;
    minv = math::Matrix({
              0.10990338, -0.0326087 , -0.00120773,
             -0.0326087 ,  0.11956522, -0.0326087 ,
             -0.00120773, -0.0326087 ,  0.10990338
        }, 3, 3);
    REQUIRE( m4.is_equal(I) );
    REQUIRE( (m3*m2).is_equal(I) );
    REQUIRE( (m2*m3).is_equal(I) );
    REQUIRE( equals(m2.determinant(), 1.0/m3.determinant()) );
    REQUIRE( m3.is_equal(minv) );
}

TEST_CASE( "Matrix: operations" ) {
    math::Matrix m1({0, 1, 2,
                     3, 4, 5,
                     6, 7, 8},
                     3, 3);
    math::Matrix m2({2, 4, 8,
                     4, 7, 6,
                     1, 7, 8},
                     3, 3);

    REQUIRE(m1 + m2 ==
            math::Matrix({
                2, 5, 10,
                7, 11, 11,
                7, 14, 16,
                }, 3, 3)
           );
    REQUIRE(m1 - m2 ==
                math::Matrix({
                -2, -3, -6,
                -1, -3, -1,
                5, 0, 0,
                }, 3, 3)
           );
    REQUIRE(m1 * m2 ==
                math::Matrix({
                6, 21, 22,
                27, 75, 88,
                48, 129, 154,
                }, 3, 3)
            );
    REQUIRE(m2 * m1 ==
                math::Matrix({
                60, 74, 88,
                57, 74, 91,
                69, 85, 101,
                }, 3, 3)
            );
    REQUIRE(5 * m1 ==
                math::Matrix({
                0, 5, 10,
                15, 20, 25,
                30, 35, 40,
                }, 3, 3) 
            );
    REQUIRE(m1 * 5 ==
                math::Matrix({
                0, 5, 10,
                15, 20, 25,
                30, 35, 40,
                }, 3, 3) 
            );
    REQUIRE(m1 / 5 ==
                math::Matrix({
                0.0, 0.2, 0.4,
                0.6, 0.8, 1.0,
                1.2, 1.4, 1.6,
                }, 3, 3)
            );
    REQUIRE(+m2 ==
                math::Matrix({
                2, 4, 8,
                4, 7, 6,
                1, 7, 8,
                }, 3, 3)
            );
    REQUIRE(-m2 ==
                math::Matrix({
                -2, -4, -8,
                -4, -7, -6,
                -1, -7, -8,
                }, 3, 3)
            );
    REQUIRE(m2 * m1 * m2 ==
                math::Matrix({
                504, 1374, 1628,
                501, 1383, 1628,
                579, 1578, 1870,
                }, 3, 3)
            );
    REQUIRE(((m2 * m2) / 7 - m1 * m2 + m2 * 6 - m2 * m1 * m2).is_equal(
                math::Matrix({
                -494.0, -1357.857142857143, -1587.142857142857,
                -498.0, -1400.7142857142858, -1662.5714285714287,
                -615.5714285714286, -1649.4285714285713, -1959.7142857142858,
                }, 3, 3))
            );
}

TEST_CASE( "Matrix: eigen" ) {
    math::Matrix m1({10, 3, 1,
                     3, 10, 3,
                     1, 3, 10},
                     3, 3);
    math::Matrix m2({
        5, -1, 0,
        -1, 10, -1,
        0, -1, 5,
        }, 3, 3);

    math::Matrix P({
        0, 0, 1,
        0, 1, 0,
        1, 0, 0,
        }, 3, 3);

    math::Eigen e1(m1);
    math::Eigen e2(m2);

    REQUIRE( e1.eigenvalues.is_equal(
        P*math::Vector({
            14.77200187,  9, 6.22799813
        }))
    );

    REQUIRE( e1.eigenvectors.is_equal(
        P*math::Matrix({
            -5.28450837e-01, -7.07106781e-01,  4.69829451e-01,
            -6.64439182e-01, -9.99640728e-17, -7.47342340e-01,
            -5.28450837e-01,  7.07106781e-01,  4.69829451e-01
        }, 3, 3)*P)
    );

    REQUIRE( e2.eigenvalues.is_equal(
        P*math::Vector({
            10.37228132, 5, 4.62771868
        }))
    );

    REQUIRE( e2.eigenvectors.is_equal(
        P*math::Matrix({
            -1.80008139e-01, -7.07106781e-01,  6.83810697e-01,
             9.67054362e-01, -4.36107634e-16,  2.54569951e-01,
            -1.80008139e-01,  7.07106781e-01,  6.83810697e-01
        }, 3, 3)*P)
    );

    const auto e1sqrt = e1.square_root();
    const auto e2sqrt = e2.square_root();

    REQUIRE( m1.is_equal(e1sqrt*e1sqrt) );
    REQUIRE( m2.is_equal(e2sqrt*e2sqrt) );

    REQUIRE( e1sqrt.is_equal(
        math::Matrix({
            3.12419619, 0.47325838, 0.12419619,
            0.47325838, 3.09063958, 0.47325838,
            0.12419619, 0.47325838, 3.12419619
        }, 3, 3))
    );

    REQUIRE( e2sqrt.is_equal(
        math::Matrix({
             2.22829199, -0.18615679, -0.00777599,
            -0.18615679,  3.15129994, -0.18615679,
            -0.00777599, -0.18615679,  2.22829199
        }, 3, 3))
    );
        

}
