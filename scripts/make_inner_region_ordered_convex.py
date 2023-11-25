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

for l in lines:
    p = l.split(" ")
    p = [float(p[0]), float(p[1]), float(p[2])]
    points.append(p)

size = len(points)

file.close()

fname = os.path.splitext(fpath)[0]
file2 = open(fname+"_inner.txt", "w")

normal = []
for i in range(size):
    im = (i + size - 1) % size
    ip = (i + 1) % size

    pm = points[im]
    pi = points[i]
    pp = points[ip]

    v1 = [pm[0] - pi[0],
          pm[1] - pi[1],
          pm[2] - pi[2]]
    v1n = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    v2 = [pp[0] - pi[0],
          pp[1] - pi[1],
          pp[2] - pi[2]]
    v2n = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
    v = [(v1[0]/v1n + v2[0]/v2n)/2,
         (v1[1]/v1n + v2[1]/v2n)/2,
         (v1[2]/v1n + v2[2]/v2n)/2]
    vnorm = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    if vnorm > 1e-7:
        normal = [
            v1[1]*v2[2] - v1[2]*v2[1],
          -(v1[0]*v2[2] - v1[2]*v2[0]),
            v1[0]*v2[1] - v1[1]*v2[0]
        ]

for i in range(size):
    im = (i + size - 1) % size
    ip = (i + 1) % size

    pm = points[im]
    pi = points[i]
    pp = points[ip]

    v1 = [pm[0] - pi[0],
          pm[1] - pi[1],
          pm[2] - pi[2]]
    v1n = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    v2 = [pp[0] - pi[0],
          pp[1] - pi[1],
          pp[2] - pi[2]]
    v2n = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
    v = [(v1[0]/v1n + v2[0]/v2n)/2,
         (v1[1]/v1n + v2[1]/v2n)/2,
         (v1[2]/v1n + v2[2]/v2n)/2]
    vnorm = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    if vnorm < 1e-7:
        v =  [
            normal[1]*v1[2] - normal[2]*v1[1],
          -(normal[0]*v1[2] - normal[2]*v1[0]),
            normal[0]*v1[1] - normal[1]*v1[0]
        ]
        vnorm = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)


    vv = [r*v[0]/vnorm, r*v[1]/vnorm, r*v[2]/vnorm]
    p = [pi[0]+vv[0], pi[1]+vv[1], pi[2]+vv[2]]
    file2.write(str(p[0]) + " " + str(p[1]) + " " + str(p[2]) + "\n")

file2.close()
