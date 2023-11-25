#!/bin/python

import FreeCAD

sub = obj.getSubObject("Face1")
f = open("points.txt", "w")

for v in sub.Vertexes:
	string = str(v.X) + " " + str(v.Y) + " " + str(v.Z) + "\n"
	f.write(string)

f.close()
