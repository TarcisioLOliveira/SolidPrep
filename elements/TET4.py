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

    formatted = str(sympy.simplify(sympy.nsimplify(sympy.expand(drnorm/8), rational=True)))
    formatted = re.sub(r"([abcdxyz]\d)\*\*2", r"\1*\1", formatted)
    formatted = re.sub(r"([abcdxyz])(\d)", r"\1[\2]", formatted)
    formatted = re.sub(r"sqrt", r"std::sqrt", formatted)

    print("const double A3 = ", formatted)
    print("")

    print("std::vector<double> Nf{")
    print("    A3,  0,  0,")
    print("     0, A3,  0,")
    print("     0,  0, A3,")
    print("    A3,  0,  0,")
    print("     0, A3,  0,")
    print("     0,  0, A3,")
    print("    A3,  0,  0,")
    print("     0, A3,  0,")
    print("     0,  0, A3,")
    print("    A3,  0,  0,")
    print("     0, A3,  0,")
    print("     0,  0, A3")
    print("};")

def make_B():
    """
        Creates the B matrix.
    """
    init_L_symbolic()
    init_N()
    init_DB()

    print("std::vector<double> B{")
    for i in range(len(B)):
        # Prepare for printing
        B[i] = sympy.simplify(sympy.nsimplify(sympy.collect(sympy.expand(B[i]), L)), rational=True)
        # Format output for use with C++
        formatted = str(B[i]).replace("(x, y)","")
        formatted = re.sub(r"([abcdD])(\d)", r"\1[\2]", formatted)

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

def make_diff():
    """
        Creates the 1 degree-of-freedom diffusion matrix.
    """

    init_L_symbolic()
    NN = sympy.Matrix([L[0], L[1], L[2], L[3]]).T
    dNN = sympy.Matrix([NN.diff(x), NN.diff(y), NN.diff(z)])
    A = sympy.symbols("A:9")
    Am = sympy.Matrix([[A[0], A[1], A[2]],
                       [A[3], A[4], A[5]],
                       [A[6], A[7], A[8]]])

    LEN = 4

    D = V*dNN.T*Am*dNN

    print("Eigen::MatrixXd M{{")
    for i in range(LEN):
        for j in range(LEN):

            D.row(i)[j] = sympy.simplify(sympy.expand(D.row(i)[j]), rational=True)

            # Format output for use with C++
            formatted = str(D.row(i)[j])
            formatted = re.sub(r"([abcdA]\d)\*\*2", r"\1*\1", formatted)
            formatted = re.sub(r"([abcdA])(\d)", r"\1[\2]", formatted)

            if j > 0:
                print(",")
            print(formatted)
        if i < LEN - 1:
            print("},{")
        else:
            print("}")
    print("};")

def make_adv():
    """
        Creates the 1 degree-of-freedom advection matrix.
    """

    init_L_symbolic()
    NN = sympy.Matrix([L[0], L[1], L[2], L[3]]).T
    dNN = sympy.Matrix([NN.diff(x), NN.diff(y), NN.diff(z)])
    v = sympy.symbols("v:3")
    vv = sympy.Matrix([v[0], v[1], v[2]])

    LEN = 4

    D = V*dNN.T*(vv*NN)

    # Prepare for "integration" by defining the variables that should be
    # collected
    LL = []
    for l1 in L:
        for l2 in L:
            LL.append(l1*l2)

    LL.extend(L)

    print("Eigen::MatrixXd M{{")
    for i in range(LEN):
        for j in range(LEN):
            H = sympy.collect(sympy.expand(D.row(i)[j]), LL, exact=True)

            # The "integration" step
            # Replace the combinations of Li*Lj with the respective constants
            # obtained by integration.
            # The formula is: area integral of (L1^a)*(L2^b)*(L3^c)*(L4^d)
            #                 = (a!b!c!d!/(a+b+c+d+3)!)*6V
            # (HUEBNER, 2001, p. 156)
            for l in L:
                H = H.subs(l*l, 2*6/(5*4*3*2))

            for l1 in L:
                for l2 in L:
                    H = H.subs(l1*l2, 1*6/(5*4*3*2))

            for l in L:
                H = H.subs(l, 6/(4*3*2))

            H = sympy.simplify(sympy.expand(H), rational=True)

            # Format output for use with C++
            formatted = str(H)
            formatted = re.sub(r"([abcdv]\d)\*\*2", r"\1*\1", formatted)
            formatted = re.sub(r"([abcdv])(\d)", r"\1[\2]", formatted)

            if j > 0:
                print(",")
            print(formatted)
        if i < LEN - 1:
            print("},{")
        else:
            print("}")
    print("};")

def make_abs():
    """
        Creates the 1 degree-of-freedom absorption matrix.
    """

    init_L_symbolic()
    NN = sympy.Matrix([L[0], L[1], L[2], L[3]]).T
    dNN = sympy.Matrix([NN.diff(x), NN.diff(y), NN.diff(z)])

    LEN = 4

    D = V*NN.T*NN

    # Prepare for "integration" by defining the variables that should be
    # collected
    LL = []
    for l1 in L:
        for l2 in L:
            LL.append(l1*l2)

    LL.extend(L)

    print("Eigen::MatrixXd M{{")
    for i in range(LEN):
        for j in range(LEN):
            H = sympy.collect(sympy.expand(D.row(i)[j]), LL, exact=True)

            # The "integration" step
            # Replace the combinations of Li*Lj with the respective constants
            # obtained by integration.
            # The formula is: area integral of (L1^a)*(L2^b)*(L3^c)*(L4^d)
            #                 = (a!b!c!d!/(a+b+c+d+3)!)*6V
            # (HUEBNER, 2001, p. 156)
            for l in L:
                H = H.subs(l*l, 2*6/(5*4*3*2))

            for l1 in L:
                for l2 in L:
                    H = H.subs(l1*l2, 1*6/(5*4*3*2))

            for l in L:
                H = H.subs(l, 6/(4*3*2))

            H = sympy.simplify(sympy.expand(H), rational=True)

            # Format output for use with C++
            formatted = str(H)
            formatted = re.sub(r"([abcdv]\d)\*\*2", r"\1*\1", formatted)
            formatted = re.sub(r"([abcdv])(\d)", r"\1[\2]", formatted)

            if j > 0:
                print(",")
            print(formatted)
        if i < LEN - 1:
            print("},{")
        else:
            print("}")
    print("};")


def make_src():
    """
        Creates the 1 degree-of-freedom source vector.
    """

    init_L_symbolic()
    NN = sympy.Matrix([L[0], L[1], L[2], L[3]]).T
    dNN = sympy.Matrix([NN.diff(x), NN.diff(y), NN.diff(z)])

    LEN = 4

    D = V*NN

    # Prepare for "integration" by defining the variables that should be
    # collected
    LL = []
    for l1 in L:
        for l2 in L:
            LL.append(l1*l2)

    LL.extend(L)

    print("Eigen::VectorXd M{")
    for i in range(LEN):
        H = sympy.collect(sympy.expand(D[i]), LL, exact=True)

        # The "integration" step
        # Replace the combinations of Li*Lj with the respective constants
        # obtained by integration.
        # The formula is: area integral of (L1^a)*(L2^b)*(L3^c)*(L4^d)
        #                 = (a!b!c!d!/(a+b+c+d+3)!)*6V
        # (HUEBNER, 2001, p. 156)
        for l in L:
            H = H.subs(l*l, 2*6/(5*4*3*2))

        for l1 in L:
            for l2 in L:
                H = H.subs(l1*l2, 1*6/(5*4*3*2))

        for l in L:
            H = H.subs(l, 6/(4*3*2))

        H = sympy.simplify(sympy.expand(H), rational=True)

        # Format output for use with C++
        formatted = str(H)
        formatted = re.sub(r"([abcv]\d)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcv])(\d)", r"\1[\2]", formatted)

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def main():
    # Backwards compatible switch statement
    args = {
        "-k":  make_k,
        "-DB": make_DB,
        "-B": make_B,
        "-Nf": make_Nf,
        "-diff": make_diff,
        "-adv": make_adv,
        "-abs": make_abs,
        "-src": make_src,
    }
    for i in range(1, len(sys.argv)):
        args[sys.argv[i]]()


if __name__ == "__main__":
    main()



