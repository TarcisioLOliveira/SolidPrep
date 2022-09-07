#!/usr/bin/python3

import numpy as np
import sympy
from sympy.core.function import *
import sympy.physics.vector as spv
import sys
import re

xi = sympy.symbols('xi')
eta = sympy.symbols('eta')
N = [0]*4
t = sympy.symbols("t")
D = []
B = sympy.Matrix(np.zeros((3, 8)))
J = 0
a = sympy.symbols('a:4')
b = sympy.symbols('b:4')
c = sympy.symbols('c:4')
d = sympy.symbols('d:4')
x, y = sympy.symbols('x, y')

def init_N():
    """ 
        Initializes the interpolation matrix N.
    """
    global N

    for i in range(4):
        N[i] = a[i] + b[i]*x + c[i]*y + d[i]*x*y

def init_DB():
    """ 
        Initializes the linear elasticity matrix B and the constitutive matrix D.
    """
    global B
    global D
    global J

    init_N()

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

def make_Nf():
    """
        Creates the Nf matrix.
    """
    init_N()

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
    init_DB()

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
        Creates the k matrix.
    """
    init_DB()

    # Create the non-integrated matrix
    k = t*B.T*D*B

    print("std::vector<double> k{")
    for i in range(len(k)):
        # Prepare for integration
        k[i] = sympy.simplify(sympy.collect(sympy.expand(k[i]), [x, y]))


        # The integration step
        k[i] = k[i].integrate((x, -1, 1))
        k[i] = k[i].integrate((y, -1, 1))

        k[i] = sympy.simplify(sympy.nsimplify(sympy.expand(k[i]), rational=True))

        # Format output for use with C++
        formatted = str(k[i])
        formatted = re.sub(r"([abcdD]\d)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcdD])(\d)", r"\1[\2]", formatted)

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def main():
    # Backwards compatible switch statement
    args = {
        "-k":  make_k,
        "-DB": make_DB,
        "-Nf": make_Nf
    }
    for i in range(1, len(sys.argv)):
        args[sys.argv[i]]()


if __name__ == "__main__":
    main()
