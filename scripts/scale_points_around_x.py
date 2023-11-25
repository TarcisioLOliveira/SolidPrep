#!/bin/python

import sys
import os
import math

if len(sys.argv) < 3:
    print("Missing arguments!")
    exit(-1)

fpath = sys.argv[1]
r = float(sys.argv[2])

file = open(fpath, "r")
lines = file.readlines()

points = []
center = [0, 0, 0]

y = [1e100, -1e100]
z = [1e100, -1e100]

for l in lines:
    p = l.split(" ")
    p = [float(p[0]), float(p[1]), float(p[2])]
    points.append(p)
    y[0] = min(y[0], p[1])
    y[1] = max(y[1], p[1])
    z[0] = min(z[0], p[2])
    z[1] = max(z[1], p[2])

file.close()

size = len(points)
center[1] = (y[1] + y[0])/2
center[2] = (z[1] + z[0])/2

sy = 1 - r/(y[1] - center[1])
sz = 1 - r/(z[1] - center[2])

fname = os.path.splitext(fpath)[0]
file2 = open(fname+"_inner.txt", "w")

for p in points:
    v = [0, sy*(p[1] - center[1]), sz*(p[2] - center[2])]
    p = [p[0], center[1]+v[1], center[2]+v[2]]
    file2.write(str(p[0]) + " " + str(p[1]) + " " + str(p[2]) + "\n")

file2.close()
