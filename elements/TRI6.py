#!/usr/bin/python3

import numpy as np
import sympy
from sympy.abc import x, y, delta
from sympy.core.function import *
import sympy.physics.vector as spv
import sys
import re

a = sympy.symbols('a:3')
b = sympy.symbols('b:3')
c = sympy.symbols('c:3')
L = [0]*6
t = sympy.symbols("t")
r = sympy.symbols("r")
N = [0]*6
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
                return b[0]/(2*delta)
            elif argindex == 2:
                return c[0]/(2*delta)
            else:
                return 0

    class L1(AppliedUndef):
        def fdiff(self, argindex=1):
            if argindex == 1:
                return b[1]/(2*delta)
            elif argindex == 2:
                return c[1]/(2*delta)
            else:
                return 0

    class L2(AppliedUndef):
        def fdiff(self, argindex=1):
            if argindex == 1:
                return b[2]/(2*delta)
            elif argindex == 2:
                return c[2]/(2*delta)
            else:
                return 0

    # Assign global variable
    L[0] = L0(x,y)
    L[1] = L1(x,y)
    L[2] = L2(x,y)

def init_L_function():
    """
        Initialize L as standard symbolic functions, for symbolic
        integration.
    """
    global L
    for i in range(3):
        L[i] = (a[i] + b[i]*x + c[i]*y)/(2*delta)

def init_N():
    """ 
        Initialize the interpolation matrix N.
    """
    global N

    for i in range(3):
        N[i] = L[i]*(2*L[i] - 1)

    for i in range(3):
        j = (i+1) % 3

        N[3 + i] = 4*L[i]*L[j]
        
def init_DB():
    """ 
        Initialize the linear elasticity matrix B and the constitutive matrix D.
    """
    global B
    global D

    B0 = [0]*12
    B1 = [0]*12
    B2 = [0]*12

    for i in range(6):
        B0[2*i+0] = sympy.simplify(sympy.expand(sympy.diff(N[i], x)))

        B1[2*i+1] = sympy.simplify(sympy.expand(sympy.diff(N[i], y)))

        B2[2*i+0] = sympy.simplify(sympy.expand(sympy.diff(N[i], y)))
        B2[2*i+1] = sympy.simplify(sympy.expand(sympy.diff(N[i], x)))

    B = sympy.Matrix([B0, B1, B2])

    d = sympy.symbols('D:9')
    D = sympy.Matrix([[d[0], d[1], d[2]],
                      [d[3], d[4], d[5]],
                      [d[6], d[7], d[8]]])

def make_Nf():
    """
        Creates the Nf matrix.
    """
    init_L_function()
    init_N()

    # NN is the full interpolation matrix
    NN1 = [0]*12
    NN2 = [0]*12
    for i in range(6):
        NN1[2*i] = N[i]
        NN2[2*i+1] = N[i]
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
        formatted = re.sub(r"(delta)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcxy]\d)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcxy])(\d)", r"\1[\2]", formatted)

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
        formatted = str(DB[i]).replace("(x, y)","")
        formatted = re.sub(r"([abcD])(\d)", r"\1[\2]", formatted)

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
    k = t*B.T*D*B

    # First integration step
    k = delta*k

    # Prepare for "integration" by defining the variables that should be
    # collected
    LL = []
    for l1 in L:
        for l2 in L:
            LL.append(l1*l2)

    for l1 in L:
        for l2 in L:
            for l3 in L:
                LL.append(l1*l2*l3)

    LL.extend(L)

    print("std::vector<double> k{")
    for i in range(len(k)):
        # Prepare for integration
        k[i] = sympy.simplify(sympy.collect(sympy.expand(k[i]), LL, exact=True))

        # The "integration" step
        # Replace the combinations of Li*Lj with the respective constants
        # obtained by integration.
        # The formula is: area integral of (L1^a)*(L2^b)*(L3^c)
        #                 = (a!b!c!/(a+b+c+2)!)*2delta
        # (HUEBNER, 2001, p. 156)
        for l in L:
            k[i] = k[i].subs(l*l*l, 3*2*2/(4*3*2))

        for l1 in L:
            for l2 in L:
                for l3 in L:
                    if l1 != l2 and l1 != l3 and l2 != l3:
                        k[i] = k[i].subs(l1*l2*l3, 2/(5*4*3*2))

        for l in L:
            k[i] = k[i].subs(l*l, 2*2/(4*3*2))

        for l1 in L:
            for l2 in L:
                if l1 != l2:
                    k[i] = k[i].subs(l1*l2, 2/(4*3*2))

        for l in L:
            k[i] = k[i].subs(l, 2/(3*2))

        k[i] = sympy.simplify(sympy.nsimplify(sympy.expand(k[i]), rational=True))
        if(k[i] == sympy.Rational(1, 3)):
            k[i] = 0 # workaround for a bug

        # Format output for use with C++
        formatted = str(k[i])
        formatted = re.sub(r"([abcD]\d)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcD])(\d)", r"\1[\2]", formatted)

        print(formatted)
        if i < len(k)-1:
            print(",")
    print("};")

def main():
    # Backwards compatible switch statement
    args = {
        "-k":  make_k,
        "-DB": make_DB,
        "-Nf": make_Nf,
        # "-h": make_h,
        # "-phir": make_phi_radial,
        # "-phig": make_phi_grad,
        # "-phiu": make_phi_unid,
        # "-phio": make_phi_omni,
    }
    for i in range(1, len(sys.argv)):
        args[sys.argv[i]]()


if __name__ == "__main__":
    main()
