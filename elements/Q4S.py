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

def make_h():
    """
        Creates the elemental Helmholtz tensor.
    """
    init_N()
    NN = sympy.Matrix([N[0], N[1], N[2], N[3]]).T
    dNN = sympy.Matrix([NN.diff(xi), NN.diff(eta)])
    k = r**2*t*dNN.T*dNN
    k2 = t*NN.T*NN

    print("std::vector<double> h{")
    for i in range(len(k)):
        # Prepare for integration
        k[i] = sympy.simplify(sympy.nsimplify(sympy.collect(sympy.expand(k[i]+k2[i]), (xi, eta, xi**2, eta**2)), rational=True))
        k[i] = sympy.integrate(k[i], (xi, -a, a), (eta, -b, b))
        k[i] = sympy.simplify(sympy.nsimplify(sympy.expand(k[i]), rational=True))

        # Format output for use with C++
        formatted = str(k[i])
        formatted = re.sub(r"(delta)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcr]\d?)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abc])(\d)", r"\1[\2]", formatted)

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def make_phi_radial():
    """
        Creates the elemental Helmholtz tensor.
    """
    init_N()
    # l = sympy.symbols("l")
    beta = sympy.symbols("beta")
    rho = sympy.symbols("rho")
    dv = sympy.symbols("dv")
    vn = sympy.symbols("vn")
    vp = sympy.symbols("vp")
    ax, ay = sympy.symbols("ax ay")
    NN = sympy.Matrix([N[0], N[1], N[2], N[3]]).T
    A = sympy.Matrix([ax, ay]).T
    x0, y0 = sympy.symbols("x0 y0")
    cx, cy = sympy.symbols("cx cy")
    P = sympy.Matrix([x0 + a + xi,
                      y0 + b + eta]).T
    C = sympy.Matrix([cx, cy]).T
    PC = C - P
    v = (PC - PC.dot(A)*A).T
    dv = v[0].diff(xi) + v[1].diff(eta)
    dNN = sympy.Matrix([NN.diff(xi), NN.diff(eta)])
    # k = t*(beta*rho*NN.T*NN + l*l*dNN.T*dNN + l*NN.T*(v.T*dNN) + dv*NN.T*NN)
    # k = t*(beta*rho*NN.T*NN + vn*dNN.T*A*dNN + vp*NN.T*(v.T*dNN) + vp*dv*NN.T*NN)
    k = t*(beta*rho*NN.T*NN + v.dot(v)*dNN.T*dNN + vp*NN.T*(v.T*dNN) + vp*dv*NN.T*NN)

    print("std::vector<double> phi{")
    for i in range(len(k)):
        # Prepare for integration
        k[i] = sympy.simplify(sympy.nsimplify(sympy.collect(sympy.expand(k[i]), (xi, eta, xi**2, eta**2)), rational=True))
        k[i] = sympy.integrate(k[i], (xi, -a, a), (eta, -b, b))
        k[i] = sympy.simplify(sympy.nsimplify(sympy.expand(k[i]), rational=True))

        # Format output for use with C++
        formatted = str(k[i])
        formatted = re.sub(r"(delta)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcdlxy]+\d?)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcdlxy]+\d?)\*\*3", r"\1*\1*\1", formatted)
        formatted = re.sub(r"([abcdlxy]+\d?)\*\*4", r"\1*\1*\1*\1", formatted)
        formatted = re.sub(r"([abcdv])(\d)", r"\1[\2]", formatted)

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def make_phi_grad():
    """
        Creates the elemental Helmholtz tensor.
    """
    init_N()
    beta = sympy.symbols("beta")
    NN = sympy.Matrix([N[0], N[1], N[2], N[3]]).T
    dNN = sympy.Matrix([NN.diff(xi), NN.diff(eta)])
    k = t*(beta*NN.T*NN)

    print("std::vector<double> phi{")
    for i in range(len(k)):
        # Prepare for integration
        k[i] = sympy.simplify(sympy.nsimplify(sympy.collect(sympy.expand(k[i]), (xi, eta, xi**2, eta**2)), rational=True))
        k[i] = sympy.integrate(k[i], (xi, -a, a), (eta, -b, b))
        k[i] = sympy.simplify(sympy.nsimplify(sympy.expand(k[i]), rational=True))

        # Format output for use with C++
        formatted = str(k[i])
        formatted = re.sub(r"(delta)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcdl]\d?)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcdv])(\d)", r"\1[\2]", formatted)

        if i > 0:
            print(",")
        print(formatted)
    print("};")

def make_phi_unid():
    """
        Creates the elemental Helmholtz tensor.
    """
    init_N()
    l = sympy.symbols("l")
    beta = sympy.symbols("beta")
    vv = sympy.symbols("v:2")
    v = sympy.Matrix([vv[0], vv[1]])
    vn = sympy.symbols("vn")
    NN = sympy.Matrix([N[0], N[1], N[2], N[3]]).T
    dNN = sympy.Matrix([NN.diff(xi), NN.diff(eta)])
    k = t*(-beta*NN.T*NN - l*l*dNN.T*dNN + l*vn*NN.T*(v.T*dNN))

    print("std::vector<double> phi{")
    for i in range(len(k)):
        # Prepare for integration
        k[i] = sympy.simplify(sympy.nsimplify(sympy.collect(sympy.expand(k[i]), (xi, eta, xi**2, eta**2)), rational=True))
        k[i] = sympy.integrate(k[i], (xi, -a, a), (eta, -b, b))
        k[i] = sympy.simplify(sympy.nsimplify(sympy.expand(k[i]), rational=True))

        # Format output for use with C++
        formatted = str(k[i])
        formatted = re.sub(r"(delta)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcdl]\d?)\*\*2", r"\1*\1", formatted)
        formatted = re.sub(r"([abcdv])(\d)", r"\1[\2]", formatted)

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
        "-h": make_h,
        "-phir": make_phi_radial,
        "-phig": make_phi_grad,
        "-phiu": make_phi_unid
    }
    for i in range(1, len(sys.argv)):
        args[sys.argv[i]]()


if __name__ == "__main__":
    main()



