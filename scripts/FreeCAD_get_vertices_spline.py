#!/bin/python

import FreeCAD

N = 100

e = obj.getSubObject("Edge1")
f = open("points.txt", "w")

begin = e.FirstParameter
end = e.LastParameter

step = (end - begin)/N

i = begin
while i < end:
    v = e.Curve.getD1(i)
    n = v[1]
    n.normalize()
    vv = v[0] - 2*n
    string = str(vv.x) + " " + str(vv.y) + " " + str(vv.z) + "\n"
    f.write(string)
    i += step

f.close()
