#!/usr/bin/python3

import numpy as np
import sympy
from sympy.abc import x, y, delta
from sympy.core.function import *
import sympy.physics.vector as spv
import sys
import re

a = sympy.symbols('a:4')
b = sympy.symbols('b:4')
c = sympy.symbols('c:4')
d = sympy.symbols('d:4')
L = [0]*4
t = sympy.symbols("t")
r = sympy.symbols("r")
N = [0]*4
D = []
B = []

def init_L_function():
    """
        Initialize L as standard symbolic functions, for symbolic
        integration.
    """
    global L
    for i in range(4):
        L[i] = a[i] + b[i]*x + c[i]*y + d[i]*x*y

def init_L_unitary():
    """
        Initialize L as unit symbolic functions, for symbolic
        integration.
    """
    global L
    L[0] = 1 - x - y - 3*x*y
    L[1] = x - 3*x*y
    L[2] = y - 3*x*y
    L[3] = 9*x*y
    

def init_N():
    """ 
        Initialize the interpolation matrix N.
    """
    global N

    for i in range(4):
        N[i] = sympy.Matrix([[L[i],    0],
                             [   0, L[i]]])

def init_DB():
    """ 
        Initialize the linear elasticity matrix B and the constitutive matrix D.
    """
    global B
    global D

    B0 = [0]*8
    B1 = [0]*8
    B2 = [0]*8

    for i in range(4):
        B0[2*i+0] = sympy.simplify(sympy.expand(sympy.diff(N[i][0,0], x)))
        B1[2*i+1] = sympy.simplify(sympy.expand(sympy.diff(N[i][1,1], y)))
        B2[2*i+0] = sympy.simplify(sympy.expand(sympy.diff(N[i][0,0], y)))
        B2[2*i+1] = sympy.simplify(sympy.expand(sympy.diff(N[i][1,1], x)))

    B = sympy.Matrix([B0, B1, B2])

    dd = sympy.symbols('D:9')
    D = sympy.Matrix([[dd[0], dd[1], dd[2]],
                      [dd[3], dd[4], dd[5]],
                      [dd[6], dd[7], dd[8]]])

def make_Nf():
    """
        Creates the Nf matrix.
    """
    init_L_function()
    init_N()

    # NN is the full interpolation matrix
    NN1 = [N[0][0,0], N[0][0,1], N[1][0,0], N[1][0,1], N[2][0,0], N[2][0,1], N[3][0,0], N[3][0,1]]
    NN2 = [N[0][1,0], N[0][1,1], N[1][1,0], N[1][1,1], N[2][1,0], N[2][1,1], N[3][1,0], N[3][1,1]]
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
    init_L_function()
    init_N()
    init_DB()

    DB = D*B

    print("std::vector<double> DB{")
    for i in range(len(DB)):
        # Prepare for printing
        DB[i] = sympy.simplify(sympy.nsimplify(sympy.collect(sympy.expand(DB[i]), L)), rational=True)
        # Format output for use with C++
        formatted = str(DB[i]).replace("(x, y)","")
        formatted = re.sub(r"([abcdD])(\d)", r"\1[\2]", formatted)

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def make_k():
    """
        Creates the k matrix.
    """
    init_L_unitary()
    init_N()
    init_DB()

    # Create the non-integrated matrix
    k = t*B.T*D*B

    # Integration part 1
    k = 2*delta*k

    print("std::vector<double> k{")
    for i in range(len(k)):

        k[i] = sympy.integrate(k[i], (y, 0, 1-x), (x, 0, 1))
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
        "-Nf": make_Nf,
        #"-h": make_h,
        #"-phir": make_phi_radial,
        #"-phig": make_phi_grad,
        #"-phiu": make_phi_unid,
    }
    for i in range(1, len(sys.argv)):
        args[sys.argv[i]]()


if __name__ == "__main__":
    main()
