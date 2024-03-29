#!/usr/bin/python

import numpy as np
import sympy
from sympy.core.function import *
import sys
import re

# Anisotropic Elasticity: Theory and Applications
# (Ting, 1996)

def format_and_print(var, val):
    val = sympy.simplify(sympy.expand(val), rational=True);
    string = str(val)
    string = re.sub(r"sqrt", r"std::sqrt", string)
    string = re.sub(r"([des]\d+)\*\*2", r"\1*\1", string)
    string = re.sub(r"([cs][ab])\*\*2", r"\1*\1", string)
    string = re.sub(r"([cs][ab])\*\*3", r"\1*\1*\1", string)
    string = re.sub(r"([cs][ab])\*\*4", r"\1*\1*\1*\1", string)
    string = re.sub(r"([des]\d+)\*\*0\.5", r"std::sqrt(\1)", string)
    string = re.sub(r"([des])(\d+)", r"\1[\2]", string)
    print("const double " + var + " = " + string + ";")

d = sympy.symbols("d:36")
e = sympy.symbols("e:36")
s = sympy.symbols("s:36")

ca, sa, cb, sb = sympy.symbols("ca sa cb sb")

D2D = sympy.Matrix([[d[0], d[1],    0],
                    [d[1], d[4],    0],
                    [   0,   0, d[8]]]);

D3D = sympy.Matrix([[d[ 0], d[ 1], d[ 2],     0,     0,     0],
                    [d[ 1], d[ 7], d[ 8],     0,     0,     0],
                    [d[ 2], d[ 8], d[14],     0,     0,     0],
                    [    0,     0,     0, d[21],     0,     0],
                    [    0,     0,     0,     0, d[28],     0],
                    [    0,     0,     0,     0,     0, d[35]]]);

S2D = sympy.Matrix([[s[0], s[1],    0],
                    [s[1], s[4],    0],
                    [   0,   0, s[8]]]);
 
S3D = sympy.Matrix([[s[ 0], s[ 1], s[ 2],     0,     0,     0],
                    [s[ 1], s[ 7], s[ 8],     0,     0,     0],
                    [s[ 2], s[ 8], s[14],     0,     0,     0],
                    [    0,     0,     0, s[21],     0,     0],
                    [    0,     0,     0,     0, s[28],     0],
                    [    0,     0,     0,     0,     0, s[35]]]);

print("2D:")

T2D = sympy.Matrix([[ca**2, sa**2, ca*sa],
                    [sa**2, ca**2, -ca*sa],
                    [-ca*sa, ca*sa, ca**2 - sa**2]])

D2DT = T2D.T*D2D*T2D
S2DT = T2D*S2D*T2D.T

format_and_print("E", D2DT[0])
format_and_print("G", D2DT[8])
format_and_print("S12", S2DT[1])

print("")
print("3D:")

T3DA = sympy.Matrix([[1, 0, 0, 0, 0, 0],
                     [0, ca**2, sa**2, 0, 0, ca*sa],
                     [0, sa**2, ca**2, 0, 0, -ca*sa],
                     [0, 0, 0, ca, -sa, 0],
                     [0, 0, 0, sa, ca, 0],
                     [0, -ca*sa, ca*sa, 0, 0, ca**2 - sa**2]])
T3DAinv = sympy.Matrix([[1, 0, 0, 0, 0, 0],
                        [0, ca**2, sa**2, 0, 0, -ca*sa],
                        [0, sa**2, ca**2, 0, 0, ca*sa],
                        [0, 0, 0, ca, sa, 0],
                        [0, 0, 0, -sa, ca, 0],
                        [0, ca*sa, -ca*sa, 0, 0, ca**2 - sa**2]])

T3DB = sympy.Matrix([[cb**2, sb**2, 0, cb*sb, 0, 0],
                     [sb**2, cb**2, 0, -cb*sb, 0, 0],
                     [0, 0, 1, 0, 0, 0],
                     [-cb*sb, cb*sb, 0, cb**2 - sb**2, 0, 0],
                     [0, 0, 0, 0, cb, sb],
                     [0, 0, 0, 0, -sb, cb]])
T3DBinv = sympy.Matrix([[cb**2, sb**2, 0, -cb*sb, 0, 0],
                        [sb**2, cb**2, 0, cb*sb, 0, 0],
                        [0, 0, 1, 0, 0, 0],
                        [cb*sb, -cb*sb, 0, cb**2 - sb**2, 0, 0],
                        [0, 0, 0, 0, cb, sb],
                        [0, 0, 0, 0, -sb, cb]])

D3DT = T3DB*T3DA*D3D*T3DA.T*T3DB.T
S3DT = T3DBinv.T*T3DAinv.T*S3D*T3DAinv*T3DBinv

format_and_print("E", D3DT[0])
format_and_print("G1", D3DT[21])
format_and_print("G2", D3DT[28])
format_and_print("G3", D3DT[35])

format_and_print("S12", S3DT[1])
format_and_print("S13", S3DT[2])

