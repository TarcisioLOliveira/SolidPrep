#!/usr/bin/python3

import numpy as np
import sympy
from sympy.abc import x, y, delta
from sympy.core.function import *
import sympy.physics.vector as spv
import sys


a = sympy.symbols('a:3')
b = sympy.symbols('b:3')
c = sympy.symbols('c:3')
L = [0]*3
t = sympy.symbols("t")
N = [0]*3
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

    Nut = [0]*3
    Nvt = [0]*3

    for i in range(3):
        j = (i+1) % 3
        k = (i+2) % 3

        Nut[i] = 0.5*L[i]*(b[k]*L[j]-b[j]*L[k])
        Nvt[i] = 0.5*L[i]*(c[k]*L[j]-c[j]*L[k])
        N[i] = sympy.Matrix([[L[i],    0, Nut[i]],
                             [   0, L[i], Nvt[i]]])

def init_DB():
    """ 
        Initialize the linear elasticity matrix B and the constitutive matrix D.
    """
    global B
    global D

    B0 = [0]*9
    B1 = [0]*9
    B2 = [0]*9

    for i in range(3):
        B0[3*i+0] = sympy.simplify(sympy.expand(sympy.diff(N[i][0,0], x)))
        B0[3*i+2] = sympy.simplify(sympy.expand(sympy.diff(N[i][0,2], x)))

    for i in range(3):
        B1[3*i+1] = sympy.simplify(sympy.expand(sympy.diff(N[i][1,1], y)))
        B1[3*i+2] = sympy.simplify(sympy.expand(sympy.diff(N[i][1,2], y)))

    for i in range(3):
        B2[3*i+0] = sympy.simplify(sympy.expand(sympy.diff(N[i][0,0], y)))
        B2[3*i+1] = sympy.simplify(sympy.expand(sympy.diff(N[i][1,1], x)))
        B2[3*i+2] = sympy.simplify(sympy.expand(sympy.diff(N[i][0,2], y) + sympy.diff(N[i][1,2], x)))

    B = sympy.Matrix([B0, B1, B2])

    d = sympy.symbols('d:9')
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
    NN1 = [N[0][0,0], N[0][0,1], N[0][0,2], N[1][0,0], N[1][0,1], N[1][0,2], N[2][0,0], N[2][0,1], N[2][0,2]]
    NN2 = [N[0][1,0], N[0][1,1], N[0][1,2], N[1][1,0], N[1][1,1], N[1][1,2], N[2][1,0], N[2][1,1], N[2][1,2]]
    NN = t*sympy.Matrix([NN1, NN2]).T

    # Set up variables for line integral
    x1, x2, y1, y2 = sympy.symbols("x1 x2 y1 y2")

    M = spv.ReferenceFrame("M")

    s = sympy.symbols("s")

    rx = ((1-s)*x1+s*x2)
    ry = ((1-s)*y1+s*y2)
    dr = (x2-x1)*M.x + (y2-y1)*M.y

    drnorm = sympy.sqrt((x2-x1)**2 + (y2-y1)**2)

    print("std::vector<double> Nf{")
    for i in range(len(NN)):
        # Apply line integration
        NN[i] = NN[i].subs({x:rx, y:ry})
        NN[i] = NN[i]*drnorm
        NN[i] = NN[i].integrate((s, 0, 1))


        # Prepare for printing
        NN[i] = sympy.simplify(NN[i])

        # Format output for use with C++
        formatted = str(NN[i])
        for j in range(3):
            ai = "a"+str(j)
            bi = "b"+str(j)
            ci = "c"+str(j)
            formatted = formatted.replace(ai+"**2", ai+"*"+ai)
            formatted = formatted.replace(bi+"**2", bi+"*"+bi)
            formatted = formatted.replace(ci+"**2", ci+"*"+ci)

        for j in range(1,3):
            xi = "x"+str(j)
            yi = "y"+str(j)
            formatted = formatted.replace(xi+"**2", xi+"*"+xi)
            formatted = formatted.replace(yi+"**2", yi+"*"+yi)

        formatted = formatted.replace("delta**2", "(delta*delta)")
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
        DB[i] = sympy.simplify(sympy.collect(sympy.expand(DB[i]), L))
        # Format output for use with C++
        formatted = str(DB[i]).replace("(x, y)","")

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

    # First integration step
    k = delta*k

    print("std::vector<double> k{")
    for i in range(len(k)):
        # Prepare for integration
        k[i] = sympy.simplify(sympy.collect(sympy.expand(k[i]), LL))

        # The "integration" step
        # Replace the combinations of Li*Lj with the respective constants
        # obtained by integration.
        # The formula is: area integral of (L1^a)*(L2^b)*(L3^c)
        #                 = (a!b!c!/(a+b+c+2)!)*2delta
        # (HUEBNER, 2001, p. 156)
        for l in L:
            k[i] = k[i].subs(l*l, 2*2/(4*3*2))

        for l1 in L:
            for l2 in L:
                if l1 != l2:
                    k[i] = k[i].subs(l1*l2, 2/(4*3*2))

        for l in L:
            k[i] = k[i].subs(l, 2/(3*2))

        k[i] = sympy.simplify(k[i])

        # Format output for use with C++
        formatted = str(k[i])
        for j in range(3):
            ai = "a"+str(j)
            bi = "b"+str(j)
            ci = "c"+str(j)
            formatted = formatted.replace(ai+"**2", ai+"*"+ai)
            formatted = formatted.replace(bi+"**2", bi+"*"+bi)
            formatted = formatted.replace(ci+"**2", ci+"*"+ci)

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

