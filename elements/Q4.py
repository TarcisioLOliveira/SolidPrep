#!/usr/bin/python3

import numpy as np
import sympy
from sympy.core.function import *
import sympy.physics.vector as spv
import sys
import re

# Partially based on: 
# MATLAB Guide to Finite Elements: An Interactive Approach (Kattan, 2008)

xi = sympy.symbols('xi')
eta = sympy.symbols('eta')
N = [0]*4
t = sympy.symbols("t")
r = sympy.symbols("r")
D = []
B = sympy.Matrix(np.zeros((3, 8)))
J = 0
x, y = sympy.symbols('x, y')

def init_N_nat():
    """ 
        Initializes the interpolation matrix N with natural coordinates.
    """
    global N

    a = sympy.symbols('a:4')
    b = sympy.symbols('b:4')
    c = sympy.symbols('c:4')
    d = sympy.symbols('d:4')

    for i in range(4):
        N[i] = a[i] + b[i]*x + c[i]*y + d[i]*x*y

def init_DB_nat():
    """ 
        Initializes the linear elasticity matrix B and the constitutive matrix 
        D with natural coordinates.
    """
    global B
    global D
    global J

    init_N_nat()

    for i in range(4):
        dNidx = sympy.diff(N[i], x)
        dNidy = sympy.diff(N[i], y)
        B[0, 2*i] = dNidx
        B[1, 2*i + 1] = dNidy
        B[2, 2*i] = dNidy
        B[2, 2*i + 1] = dNidx

    d = sympy.symbols('D:9')
    D = sympy.Matrix([[d[0], d[1], d[2]],
                      [d[3], d[4], d[5]],
                      [d[6], d[7], d[8]]])

def init_N_norm():
    """ 
        Initializes the interpolation matrix N with normalized coordinates.
    """
    global N

    N[0] = 0.25*(1-xi)*(1-eta)
    N[1] = 0.25*(1+xi)*(1-eta)
    N[2] = 0.25*(1+xi)*(1+eta)
    N[3] = 0.25*(1-xi)*(1+eta)

def init_DB_norm():
    """ 
        Initializes the linear elasticity matrix B and the constitutive matrix 
        D with normalized coordinates.
    """
    global B
    global D
    global J

    init_N_norm()

    x = sympy.symbols('x:4')
    y = sympy.symbols('y:4')

    a = sympy.Rational(0.25)*(y[0]*(xi-1)+ y[1]*(-1-xi)+y[2]*(1+xi) +y[3]*(1-xi))
    b = sympy.Rational(0.25)*(y[0]*(eta-1)+y[1]*(1-eta)+y[2]*(1+eta)+y[3]*(-1-eta))
    c = sympy.Rational(0.25)*(x[0]*(eta-1)+x[1]*(1-eta)+x[2]*(1+eta)+x[3]*(-1-eta))
    d = sympy.Rational(0.25)*(x[0]*(xi-1)+ x[1]*(-1-xi)+x[2]*(1+xi) +x[3]*(1-xi))
    a = sympy.expand(a)
    b = sympy.expand(b)
    c = sympy.expand(c)
    d = sympy.expand(d)

    for i in range(4):
        dNidxi = sympy.diff(N[i], xi)
        dNideta = sympy.diff(N[i], eta)
        B[0, 2*i] = a*dNidxi - b*dNideta
        B[1, 2*i + 1] = c*dNideta - d*dNidxi
        B[2, 2*i] = c*dNideta - d*dNidxi
        B[2, 2*i + 1] = a*dNidxi - b*dNideta

    yy = sympy.Matrix([y[0], y[1], y[2], y[3]])
    xx = sympy.Matrix([x[0], x[1], x[2], x[3]])

    Jn = sympy.Matrix(
            [[0, 1-eta, eta-xi, xi-1],
             [eta-1, 0, xi+1, -xi-eta],
             [xi-eta, -xi-1, 0, eta+1],
             [1-xi, xi+eta, -eta-1, 0]])

    J = sympy.simplify(sympy.collect(sympy.expand((sympy.Rational(1, 8)*xx.T*(Jn*yy))[0]), [xi**2, eta**2, xi*eta, xi, eta]), rational=True)
    B = sympy.expand(B)

    d = sympy.symbols('D:9')
    D = sympy.Matrix([[d[0], d[1], d[2]],
                      [d[3], d[4], d[5]],
                      [d[6], d[7], d[8]]])

def make_B():
    """
        Creates the B matrix.
    """
    init_DB_nat()

    print("std::vector<double> B{")
    for i in range(len(B)):
        # Prepare for printing
        B[i] = sympy.simplify(sympy.nsimplify(sympy.collect(sympy.expand(B[i]), [x, y])), rational=True)
        # Format output for use with C++
        formatted = str(B[i])
        formatted = re.sub(r"([abcdD])(\d)", r"\1[\2]", formatted)

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def make_DB():
    """
        Creates the DB matrix.
    """
    init_DB_nat()

    DB = D*B

    print("std::vector<double> DB{")
    for i in range(len(DB)):
        # Prepare for printing
        DB[i] = sympy.simplify(sympy.nsimplify(sympy.collect(sympy.expand(DB[i]), [x, y])), rational=True)
        # Format output for use with C++
        formatted = str(DB[i])
        formatted = re.sub(r"([abcdD])(\d)", r"\1[\2]", formatted)

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def print_matrix(M, name):
    """
        Prints the matrix M, item by item, in C++ format, with the specified
        variable name.
    """

    print("std::vector<double> "+name+"{")
    for i in range(len(M)):
        M[i] = sympy.simplify(sympy.nsimplify(sympy.expand(M[i]), rational=True))

        # Format output for use with C++
        formatted = str(M[i])
        formatted = re.sub(r"(eta)\*\*3", r"\1*\1*\1", formatted)
        formatted = re.sub(r"(xi)\*\*3", r"\1*\1*\1", formatted)
        formatted = re.sub(r"(eta)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"(xi)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([xyr]\d?)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([xy])(\d)", r"\1[\2]", formatted)

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def main():
    # Backwards compatible switch statement
    args = {
        "-DB": make_DB,
        "-B": make_B,
    }
    for i in range(1, len(sys.argv)):
        args[sys.argv[i]]()


if __name__ == "__main__":
    main()
