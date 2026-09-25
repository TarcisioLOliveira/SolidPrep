#!/bin/python

import numpy as np
import scipy
from scipy.spatial.transform import Rotation


#Rf = np.asarray([[0.468451786249991,	-0.883489062727543,	0],
#                 [0.883489062727543,	 0.468451786249991,	0],
#                 [ 0                ,   0                ,  1]])

Rf = np.asarray([[0.16443142319296,  0.98638851730276, 0],
                 [-0.98638851730276, 0.16443142319296, 0],
                 [0, 0, 1]])

Rz = np.asarray([
    [1, 0, 0],
    [0,-1, 0],
    [0, 0,-1]
    ])
R0 = Rotation.from_rotvec(-85.95*np.array([1, 0, 0]), degrees=True)
R1 = Rotation.from_rotvec( +7.80*np.array([0, 1, 0]), degrees=True)
R2 = Rotation.from_rotvec( +0.70*np.array([0, 0, 1]), degrees=True)

#R = Rf @ R2.as_matrix() @ R1.as_matrix() @ R0.as_matrix() @ Rz
#R = Rf @ R2.as_matrix() @ R1.as_matrix() @ R0.as_matrix() @ Rz

# Inverted process relative to what I did in FreeCAD, might be related to how
# they store successive rotations
R = Rf @ R0.as_matrix() @ R1.as_matrix() @ R2.as_matrix() # @ Rz
print(R)
print(scipy.linalg.det(R))

O = np.asarray([12.7, 4.5125, -68.5])
P = np.asarray([-60.99918094865, -26.31951913311, -23.673946289])

S = np.asarray([
[1/0.3125, 0, 0],
[0, 1/0.3125, 0],
[0, 0, 1/1.25]
])

print(R.T @ P)
print(R.T @ P - O)
print(S @ (R.T @ P - O))
print(S @ R.T @ P - O)



# test = [[1, 417.382, -14.7079, 1006.41],
#         [1, 417.336, -14.9332, 1006.64],
#         [1, 417.359, -15.0605, 1006.26],
#         [1, 417.615, -14.9718, 1006.57]]
# 
# print(scipy.linalg.inv(test))
