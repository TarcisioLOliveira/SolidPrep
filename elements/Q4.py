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

def make_Nf():
    """
        Creates the Nf matrix.
    """
    init_N_nat()

    NN1 = [N[0],   0 , N[1],   0 , N[2],   0 , N[3],   0 ]
    NN2 = [  0 , N[0],   0 , N[1],   0 , N[2],   0 , N[3]] 
    NN = t*sympy.Matrix([NN1, NN2]).T

    # Set up variables for line integral
    x0, x1, y0, y1 = sympy.symbols("x0 x1 y0 y1")

    M = spv.ReferenceFrame("M")

    s = sympy.symbols("s")

    rx = ((1-s)*x0+s*x1)
    ry = ((1-s)*y0+s*y1)
    dr = (x1-x0)*M.x + (y1-y0)*M.y

    drnorm = sympy.sqrt((x1-x0)**2 + (y1-y0)**2)

    print("std::vector<double> Nf{")
    for i in range(len(NN)):
        # Apply line integration
        NN[i] = NN[i].subs({x:rx, y:ry})
        NN[i] = NN[i]*drnorm
        NN[i] = NN[i].integrate((s, 0, 1))


        # Prepare for printing
        NN[i] = sympy.simplify(sympy.nsimplify(NN[i], rational=True))

        # Format output for use with C++
        formatted = str(NN[i])
        formatted = re.sub(r"([abcdxy]\d)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcdxy])(\d)", r"\1[\2]", formatted)

        formatted = formatted.replace("sqrt", "std::sqrt")

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

def make_k():
    """
        Creates the k matrix (non-integrated, as numerical integration is
        necessary).
    """
    init_DB_norm()

    # Create the non-integrated matrix
    k = sympy.Matrix(t*B.T*D*B/J)

    print("std::vector<double> k{")
    for i in range(len(k)):
        # No integration step in this case, as the integral is too complex
        k[i] = sympy.simplify(sympy.collect(sympy.expand(k[i]), [xi**2, eta**2, xi*eta, xi, eta]), rational=True)

        # Format output for use with C++
        formatted = str(k[i])
        formatted = re.sub(r"([xyD]\d)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"(\w+)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([xyD])(\d)", r"\1[\2]", formatted)

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def init_J_matrix():
    """
        Initializes the global variable J as a matrix instead of as the
        determinant of J.
    """
    global J

    x = sympy.symbols('x:4')
    y = sympy.symbols('y:4')

    xx = sympy.Rational(1,4)*((1-xi)*(1-eta)*x[0]+(1+xi)*(1-eta)*x[1]+(1+xi)*(1+eta)*x[2]+(1-xi)*(1+eta)*x[3])
    yy = sympy.Rational(1,4)*((1-xi)*(1-eta)*y[0]+(1+xi)*(1-eta)*y[1]+(1+xi)*(1+eta)*y[2]+(1-xi)*(1+eta)*y[3])

    J = sympy.Matrix([[xx.diff(xi), yy.diff(xi)],[xx.diff(eta), yy.diff(eta)]])

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

def make_h():
    """
        Creates the elemental Helmholtz tensor.
    """
    init_N_norm()
    init_J_matrix()
    detJ = sympy.nsimplify(sympy.expand(J.det()), rational=True)
    invJ = J.inv()

    NN = sympy.Matrix([N[0], N[1], N[2], N[3]]).T
    dNN = sympy.Matrix([NN.diff(xi), NN.diff(eta)])
    for i in range(len(dNN)):
        dNN[i] = sympy.expand(sympy.nsimplify(dNN[i], rational=True))

    k = (r**2)*t*dNN.T*dNN
    k2 = t*NN.T*NN*detJ

    print_matrix(J, "J")
    print_matrix(sympy.Matrix([detJ]), "detJ")
    print_matrix(dNN, "dN")
    print_matrix(k2, "h2")

def main():
    # Backwards compatible switch statement
    args = {
        "-k":  make_k,
        "-DB": make_DB,
        "-Nf": make_Nf,
        "-h": make_h
    }
    for i in range(1, len(sys.argv)):
        args[sys.argv[i]]()


if __name__ == "__main__":
    main()
