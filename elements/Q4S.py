#!/usr/bin/python3

import numpy as np
import sympy
from sympy.core.function import *
import sympy.physics.vector as spv
import sys
import re

# Partially based on: 
# https://github.com/williamhunter/topy/blob/master/topy/data/Q4_K.py

xi = sympy.symbols('xi')
eta = sympy.symbols('eta')
a, b = sympy.symbols('a b')
N = [0]*4
t = sympy.symbols("t")
r = sympy.symbols("r")
D = []
B = sympy.Matrix(np.zeros((3, 8)))

def init_N():
    """ 
        Initializes the interpolation matrix N with normalized coordinates.
    """
    global N

    N[0] = (a-xi)*(b-eta)/(4*a*b)
    N[1] = (a+xi)*(b-eta)/(4*a*b)
    N[2] = (a+xi)*(b+eta)/(4*a*b)
    N[3] = (a-xi)*(b+eta)/(4*a*b)

def init_DB():
    """ 
        Initializes the linear elasticity matrix B and the constitutive matrix 
        D with normalized coordinates.
    """

    global B
    global D

    init_N()

    for i in range(4):
        dNidx = sympy.diff(N[i], xi)
        dNidy = sympy.diff(N[i], eta)
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
        NN[i] = NN[i].subs({xi:rx, eta:ry})
        NN[i] = NN[i]*drnorm
        NN[i] = NN[i].integrate((s, 0, 1))


        # Prepare for printing
        NN[i] = sympy.simplify(sympy.nsimplify(sympy.expand(NN[i]), rational=True))

        # Format output for use with C++
        formatted = str(NN[i])
        formatted = re.sub(r"([abcdxy]\d)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcdxy])(\d)", r"\1[\2]", formatted)

        formatted = formatted.replace("sqrt", "std::sqrt")

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def make_B():
    """
        Creates the B matrix.
    """
    init_DB()

    print("std::vector<double> B{")
    for i in range(len(B)):
        # Prepare for printing
        B[i] = sympy.simplify(sympy.expand(B[i]), rational=True)
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
    init_DB()

    DB = D*B

    print("std::vector<double> DB{")
    for i in range(len(DB)):
        # Prepare for printing
        DB[i] = sympy.simplify(sympy.expand(DB[i]), rational=True)
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
    k = sympy.Matrix(t*B.T*D*B)

    print("std::vector<double> k{")
    for i in range(len(k)):
        # Integration step
        k[i] = sympy.integrate(k[i], (xi, -a, a), (eta, -b, b))
        k[i] = sympy.ratsimp(sympy.simplify(sympy.expand(k[i]), rational=True))

        # Format output for use with C++
        formatted = str(k[i])
        formatted = re.sub(r"([ab])\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([D]\d)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([D])(\d)", r"\1[\2]", formatted)

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def make_R():
    """
        Creates the R matrix.
    """
    init_N()

    # NN is the full interpolation matrix
    NN1 = [N[0],   0 , N[1],   0 , N[2],   0 , N[3],   0 ]
    NN2 = [  0 , N[0],   0 , N[1],   0 , N[2],   0 , N[3]] 
    NN = sympy.Matrix([NN1, NN2]).T
    kk = sympy.symbols('K:4')
    K = sympy.Matrix(np.asarray(kk).reshape(2,2))
    R = t*NN*K*NN.T

    # Set up variables for line integral
    x0, x1, y0, y1 = sympy.symbols("x0 x1 y0 y1")

    M = spv.ReferenceFrame("M")

    s = sympy.symbols("s")

    rx = ((1-s)*x0+s*x1)
    ry = ((1-s)*y0+s*y1)
    dr = (x1-x0)*M.x + (y1-y0)*M.y

    drnorm = sympy.sqrt((x1-x0)**2 + (y1-y0)**2)

    print("std::vector<double> R{")
    for i in range(len(R)):
        # Apply line integration
        R[i] = R[i].subs({xi:rx, eta:ry})
        R[i] = R[i]*drnorm
        R[i] = R[i].integrate((s, 0, 1))


        # Prepare for printing
        R[i] = sympy.simplify(sympy.nsimplify(R[i], rational=True))

        # Format output for use with C++
        formatted = str(R[i])
        formatted = re.sub(r"([Kabxy]\d?)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([Kabxy])(\d)", r"\1[\2]", formatted)

        formatted = formatted.replace("sqrt", "std::sqrt")

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def make_Rf():
    """
        Creates the Rf matrix.
    """
    init_N()

    NN1 = [N[0],   0 , N[1],   0 , N[2],   0 , N[3],   0 ]
    NN2 = [  0 , N[0],   0 , N[1],   0 , N[2],   0 , N[3]] 
    NN = t*sympy.Matrix([NN1, NN2]).T

    # Set up variables for line integral
    x0, x1, y0, y1 = sympy.symbols("x0 x1 y0 y1")

    M = spv.ReferenceFrame("M")

    s = sympy.symbols("s")
    cx, cy = sympy.symbols("cx cy")

    xyz = sympy.Matrix([xi - cx, eta - cy])

    Sl = sympy.symbols("S:4")
    S = sympy.Matrix([[Sl[0], Sl[1]],
                      [Sl[2], Sl[3]]])
    Fl = sympy.symbols("F:2")
    F = sympy.Matrix([Fl[0],Fl[1]])

    NN = NN*(S*xyz + F)

    rx = ((1-s)*x0+s*x1)
    ry = ((1-s)*y0+s*y1)
    dr = (x1-x0)*M.x + (y1-y0)*M.y

    drnorm = sympy.sqrt((x1-x0)**2 + (y1-y0)**2)

    print("std::vector<double> Rf{")
    for i in range(len(NN)):
        # Apply line integration
        NN[i] = NN[i].subs({xi:rx, eta:ry})
        NN[i] = NN[i]*drnorm
        NN[i] = NN[i].integrate((s, 0, 1))


        # Prepare for printing
        NN[i] = sympy.simplify(sympy.nsimplify(sympy.expand(NN[i]), rational=True))

        # Format output for use with C++
        formatted = str(NN[i])
        formatted = re.sub(r"([abcdxySF]\d)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcdxySF])(\d)", r"\1[\2]", formatted)

        formatted = formatted.replace("sqrt", "std::sqrt")

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def make_diff():
    """
        Creates the 1 degree-of-freedom diffusion matrix.
    """

    init_N()
    NN = sympy.Matrix([N[0], N[1], N[2], N[3]]).T
    dNN = sympy.Matrix([NN.diff(xi), NN.diff(eta)])
    A = sympy.symbols("A:9")
    Am = sympy.Matrix([[A[0], A[1]],
                       [A[3], A[4]]])

    LEN = 4

    D = t*dNN.T*Am*dNN

    print("Eigen::MatrixXd M{{")
    for i in range(LEN):
        for j in range(LEN):

            H = sympy.integrate(D.row(i)[j], (xi, -a, a), (eta, -b, b))

            H = sympy.simplify(sympy.expand(H), rational=True)

            # Format output for use with C++
            formatted = str(H)
            formatted = re.sub(r"([abcA]\d)\*\*2", r"\1*\1", formatted)
            formatted = re.sub(r"([abcA])(\d)", r"\1[\2]", formatted)

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

    init_N()
    NN = sympy.Matrix([N[0], N[1], N[2], N[3]]).T
    dNN = sympy.Matrix([NN.diff(xi), NN.diff(eta)])
    v = sympy.symbols("v:3")
    vv = sympy.Matrix([v[0], v[1]])

    LEN = 4

    D = t*dNN.T*(vv*NN)

    print("Eigen::MatrixXd M{{")
    for i in range(LEN):
        for j in range(LEN):
            H = sympy.integrate(D.row(i)[j], (xi, -a, a), (eta, -b, b))

            H = sympy.simplify(sympy.expand(H), rational=True)

            # Format output for use with C++
            formatted = str(H)
            formatted = re.sub(r"([abcv]\d)\*\*2", r"\1*\1", formatted)
            formatted = re.sub(r"([abcv])(\d)", r"\1[\2]", formatted)

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

    init_N()
    NN = sympy.Matrix([N[0], N[1], N[2], N[3]]).T
    dNN = sympy.Matrix([NN.diff(xi), NN.diff(eta)])

    LEN = 4

    D = t*NN.T*NN

    print("Eigen::MatrixXd M{{")
    for i in range(LEN):
        for j in range(LEN):
            H = sympy.integrate(D.row(i)[j], (xi, -a, a), (eta, -b, b))

            H = sympy.simplify(sympy.expand(H), rational=True)

            # Format output for use with C++
            formatted = str(H)
            formatted = re.sub(r"([abcv]\d)\*\*2", r"\1*\1", formatted)
            formatted = re.sub(r"([abcv])(\d)", r"\1[\2]", formatted)

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

    init_N()
    NN = sympy.Matrix([N[0], N[1], N[2], N[3]]).T
    dNN = sympy.Matrix([NN.diff(xi), NN.diff(eta)])

    LEN = 4

    D = t*NN

    print("Eigen::VectorXd M{")
    for i in range(LEN):
        H = sympy.integrate(D[i], (xi, -a, a), (eta, -b, b))

        H = sympy.simplify(sympy.expand(H), rational=True)

        # Format output for use with C++
        formatted = str(H)
        formatted = re.sub(r"([abcv]\d)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcv])(\d)", r"\1[\2]", formatted)

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def make_flow():
    """
        Creates the 1 degree-of-freedom flow vector.
    """
    init_N()

    NN = t*sympy.Matrix([N[0], N[1], N[2], N[3]]).T

    # Set up variables for line integral
    x0, x1, y0, y1 = sympy.symbols("x0 x1 y0 y1")

    M = spv.ReferenceFrame("M")

    s = sympy.symbols("s")

    rx = ((1-s)*x0+s*x1)
    ry = ((1-s)*y0+s*y1)
    dr = (x1-x0)*M.x + (y1-y0)*M.y

    drnorm = sympy.sqrt((x1-x0)**2 + (y1-y0)**2)

    print("Eigen::VectorXd M{")
    for i in range(len(NN)):
        # Apply line integration
        NN[i] = NN[i].subs({xi:rx, eta:ry})
        NN[i] = NN[i]*drnorm
        NN[i] = NN[i].integrate((s, 0, 1))


        # Prepare for printing
        NN[i] = sympy.simplify(sympy.nsimplify(sympy.expand(NN[i]), rational=True))

        # Format output for use with C++
        formatted = str(NN[i])
        formatted = re.sub(r"([abcdxy]\d)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcdxy])(\d)", r"\1[\2]", formatted)

        formatted = formatted.replace("sqrt", "std::sqrt")

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
        "-R": make_R,
        "-Rf": make_Rf,
        "-diff": make_diff,
        "-adv": make_adv,
        "-abs": make_abs,
        "-src": make_src,
        "-flow": make_flow,
    }
    for i in range(1, len(sys.argv)):
        args[sys.argv[i]]()


if __name__ == "__main__":
    main()



