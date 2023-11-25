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

for l in lines:
    p = l.split(" ")
    p = [float(p[0]), float(p[1]), float(p[2])]
    points.append(p)
    center[0] += p[0]
    center[1] += p[1]
    center[2] += p[2]

size = len(points)
center[0] /= size
center[1] /= size
center[2] /= size

file.close()

fname = os.path.splitext(fpath)[0]
file2 = open(fname+"_inner.txt", "w")

for p in points:
    v = [center[0] - p[0], center[1] - p[1], center[2] - p[2]]
    vnorm = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    v = [r*v[0]/vnorm, r*v[1]/vnorm, r*v[2]/vnorm]
    p = [p[0]+v[0], p[1]+v[1], p[2]+v[2]]
    file2.write(str(p[0]) + " " + str(p[1]) + " " + str(p[2]) + "\n")

file2.close()
