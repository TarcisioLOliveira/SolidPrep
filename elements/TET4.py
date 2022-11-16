#!/usr/bin/python3

import numpy as np
import sympy
from sympy.abc import x, y, z
from sympy.core.function import *
import sympy.physics.vector as spv
import sys
import re


a = sympy.symbols('a:4')
b = sympy.symbols('b:4')
c = sympy.symbols('c:4')
d = sympy.symbols('d:4')
L = [0]*4
V = sympy.symbols("V")
r = sympy.symbols("r")
N = [0]*4
D = []
B = []

def init_L_symbolic():
    """ 
        Initializes the interpolation functions L as symbols, with standard
        derivation but 'replacement' integration (see `make_k()`).
    """
    global L

    # Define the classes for L and their behavior for derivation
    class L0(AppliedUndef):
        def fdiff(self, argindex=1):
            if argindex == 1:
                return b[0]/(6*V)
            elif argindex == 2:
                return c[0]/(6*V)
            elif argindex == 3:
                return d[0]/(6*V)
            else:
                return 0

    class L1(AppliedUndef):
        def fdiff(self, argindex=1):
            if argindex == 1:
                return b[1]/(6*V)
            elif argindex == 2:
                return c[1]/(6*V)
            elif argindex == 3:
                return d[1]/(6*V)
            else:
                return 0

    class L2(AppliedUndef):
        def fdiff(self, argindex=1):
            if argindex == 1:
                return b[2]/(6*V)
            elif argindex == 2:
                return c[2]/(6*V)
            elif argindex == 3:
                return d[2]/(6*V)
            else:
                return 0

    class L3(AppliedUndef):
        def fdiff(self, argindex=1):
            if argindex == 1:
                return b[3]/(6*V)
            elif argindex == 2:
                return c[3]/(6*V)
            elif argindex == 3:
                return d[3]/(6*V)
            else:
                return 0

    # Assign global variable
    L[0] = L0(x,y,z)
    L[1] = L1(x,y,z)
    L[2] = L2(x,y,z)
    L[3] = L3(x,y,z)

def init_L_function():
    """
        Initialize L as standard symbolic functions, for symbolic
        integration.
    """
    global L
    for i in range(4):
        L[i] = (a[i] + b[i]*x + c[i]*y + d[i]*z)/(6*V)

def init_N():
    """ 
        Initialize the interpolation matrix N.
    """
    global N

    for i in range(4):
        j = (i+1) % 4
        k = (i+2) % 4

        N[i] = sympy.Matrix([[L[i],    0,   0],
                             [   0, L[i],   0],
                             [   0,    0, L[i]]])

def init_DB():
    """ 
        Initialize the linear elasticity matrix B and the constitutive matrix D.
    """
    global B
    global D

    B0 = [0]*12
    B1 = [0]*12
    B2 = [0]*12
    B3 = [0]*12
    B4 = [0]*12
    B5 = [0]*12

    for i in range(4):
        B0[3*i+0] = sympy.simplify(sympy.expand(sympy.diff(N[i][0,0], x)))

    for i in range(4):
        B1[3*i+1] = sympy.simplify(sympy.expand(sympy.diff(N[i][1,1], y)))

    for i in range(4):
        B2[3*i+2] = sympy.simplify(sympy.expand(sympy.diff(N[i][2,2], z)))

    for i in range(4):
        B3[3*i+0] = sympy.simplify(sympy.expand(sympy.diff(N[i][0,0], y)))
        B3[3*i+1] = sympy.simplify(sympy.expand(sympy.diff(N[i][1,1], x)))

    for i in range(4):
        B4[3*i+1] = sympy.simplify(sympy.expand(sympy.diff(N[i][1,1], z)))
        B4[3*i+2] = sympy.simplify(sympy.expand(sympy.diff(N[i][2,2], y)))

    for i in range(4):
        B5[3*i+0] = sympy.simplify(sympy.expand(sympy.diff(N[i][0,0], z)))
        B5[3*i+2] = sympy.simplify(sympy.expand(sympy.diff(N[i][2,2], x)))

    B = sympy.Matrix([B0, B1, B2, B3, B4, B5])

    dd = sympy.symbols('D:36')
    D = sympy.Matrix(np.asarray(dd).reshape(6,6))

def make_Nf():
    """
        Creates the Nf matrix.
    """
    init_L_function()
    init_N()

    # NN is the full interpolation matrix
    NN1 = []
    NN2 = []
    NN3 = []
    for i in range(4):
        for j in range(3):
            NN1.append(N[i][0,j])
            NN2.append(N[i][1,j])
            NN3.append(N[i][2,j])

    NN = sympy.Matrix([NN1, NN2, NN3]).T

    # Set up variables for line integral
    x0, x1, x2, y0, y1, y2, z0, z1, z2 = sympy.symbols("x0 x1 x2 y0 y1 y2 z0 z1 z2")

    M = spv.ReferenceFrame("M")

    u, v = sympy.symbols("u v")

    rx = x0 + (x1-x0)*u + (x2-x0)*v
    ry = y0 + (y1-y0)*u + (y2-y0)*v
    rz = z0 + (z1-z0)*u + (z2-z0)*v
    v1 = sympy.Matrix([(x1-x0), (y1-y0), (z1-z0)])
    v2 = sympy.Matrix([(x2-x0), (y2-y0), (z2-z0)])

    vv = v1.cross(v2)
    drnorm = sympy.sqrt(vv.dot(vv))

    print("std::vector<double> Nf{")
    for i in range(len(NN)):
        # Apply line integration
        NN[i] = NN[i].subs({x:rx, y:ry, z:rz})
        NN[i] = NN[i]*drnorm
        NN[i] = NN[i].integrate((u, 0, 1),(v, 0, 1))


        # Prepare for printing
        NN[i] = sympy.simplify(sympy.nsimplify(NN[i], rational=True))

        # Format output for use with C++
        formatted = str(NN[i])
        formatted = re.sub(r"([abcdxyz]\d)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcdxyz])(\d)", r"\1[\2]", formatted)

        formatted = formatted.replace("sqrt", "std::sqrt")

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def make_DB():
    """
        Creates the DB matrix.
    """
    init_L_symbolic()
    init_N()
    init_DB()

    DB = D*B

    print("std::vector<double> DB{")
    for i in range(len(DB)):
        # Prepare for printing
        DB[i] = sympy.simplify(sympy.nsimplify(sympy.collect(sympy.expand(DB[i]), L)), rational=True)
        # Format output for use with C++
        formatted = str(DB[i]).replace("(x, y, z)","")
        formatted = re.sub(r"([abcdD])(\d+)", r"\1[\2]", formatted)

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def make_k():
    """
        Creates the k matrix.
    """
    init_L_symbolic()
    init_N()
    init_DB()

    # Create the non-integrated matrix
    k = B.T*D*B

    # Integration
    k = V*k

    print("std::vector<double> k{")
    for i in range(len(k)):

        k[i] = sympy.simplify(sympy.nsimplify(sympy.expand(k[i]), rational=True))

        # Format output for use with C++
        formatted = str(k[i])
        formatted = re.sub(r"([abcdD]\d+)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcdD])(\d+)", r"\1[\2]", formatted)

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def make_h():
    """
        Creates the elemental Helmholtz tensor.
    """
    init_L_symbolic()
    NN = sympy.Matrix([L[0], L[1], L[2], L[3]]).T
    dNN = sympy.Matrix([NN.diff(x), NN.diff(y)])
    k = r**2*dNN.T*dNN*V
    k2 = V*NN.T*NN

    # Prepare for "integration" by defining the variables that should be
    # collected
    LL = []
    for l1 in L:
        for l2 in L:
            LL.append(l1*l2)

    LL.extend(L)

    print("std::vector<double> h{")
    for i in range(len(k)):
        # Prepare for integration
        k2[i] = sympy.simplify(sympy.collect(sympy.expand(k2[i]), LL))

        # The "integration" step
        # Replace the combinations of Li*Lj with the respective constants
        # obtained by integration.
        # The formula is: area integral of (L1^a)*(L2^b)*(L3^c)
        #                 = (a!b!c!/(a+b+c+2)!)*2delta
        # (HUEBNER, 2001, p. 156)
        for l in L:
            k2[i] = k2[i].subs(l*l, 2*2/(4*3*2))

        for l in L:
            k2[i] = k2[i].subs(l, 2/(3*2))

        k[i] = sympy.simplify(sympy.nsimplify(sympy.expand(k[i]+k2[i]), rational=True))

        # Format output for use with C++
        formatted = str(k[i])
        formatted = re.sub(r"([abcrV]\d?)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abc])(\d)", r"\1[\2]", formatted)

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
        "-h": make_h
    }
    for i in range(1, len(sys.argv)):
        args[sys.argv[i]]()


if __name__ == "__main__":
    main()


