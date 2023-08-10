#!/usr/bin/python3

import numpy as np
import sympy
from sympy.core.function import *
import sys
import re

def format_and_print(mat):
    print("std::vector<double> D = {")
    for i in range(len(mat)):
        mat[i] = sympy.simplify(mat[i], rational=True);
        string = str(mat[i])
        string = re.sub(r"sqrt", r"std::sqrt", string)
        string = re.sub(r"([de]\d+)\*\*2", r"\1*\1", string)
        string = re.sub(r"([de]\d+)\*\*0\.5", r"std::sqrt(\1)", string)
        string = re.sub(r"([de])(\d+)", r"\1[\2]", string)
        print(string)
        if i < len(mat) - 1:
            print(",")
    print("};")

d = sympy.symbols("d:36")
e = sympy.symbols("e:36")

D2D = sympy.Matrix([[d[0], d[1],    0],
                    [d[3], d[4],    0],
                    [   0,   0, d[8]]]);

D3D = sympy.Matrix([[d[ 0], d[ 1], d[ 2],     0,     0,     0],
                    [d[ 6], d[ 7], d[ 8],     0,     0,     0],
                    [d[12], d[13], d[14],     0,     0,     0],
                    [    0,     0,     0, d[21],     0,     0],
                    [    0,     0,     0,     0, d[28],     0],
                    [    0,     0,     0,     0,     0, d[35]]]);

E2D = sympy.Matrix([[e[0], e[1],    0],
                    [e[3], e[4],    0],
                    [   0,   0, e[8]]]);

E3D = sympy.Matrix([[e[ 0], e[ 1], e[ 2],     0,     0,     0],
                    [e[ 6], e[ 7], e[ 8],     0,     0,     0],
                    [e[12], e[13], e[14],     0,     0,     0],
                    [    0,     0,     0, e[21],     0,     0],
                    [    0,     0,     0,     0, e[28],     0],
                    [    0,     0,     0,     0,     0, e[35]]]);

format_and_print(D2D.inv())
print("")
format_and_print(D3D.inv())
print("")
format_and_print(D2D*E2D)
print("")
format_and_print(D3D*E3D)
print("")
format_and_print(D2D*D2D)
print("")
format_and_print(D3D*D3D)
print("")
format_and_print(D2D**(1/2))
print("")

exit(0)
E3D[:3,:3] = D3D[:3,:3]**(1/2)
E3D[3,3] = sympy.sqrt(D3D[3,3])
E3D[4,4] = sympy.sqrt(D3D[4,4])
E3D[5,5] = sympy.sqrt(D3D[5,5])

format_and_print(E3D)
print("")
