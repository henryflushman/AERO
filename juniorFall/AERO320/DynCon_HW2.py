import DynCon_Functions as dc

import os


os.system('cls')


""" Problem 1.2 """
import sympy as sp
# Constants
Ax = sp.Matrix([['0','-az','ay'],
              ['az','0','-ax'],
              ['-ay','ax','0']])

#Analysis
Ax_cube = Ax ** 3
print(f"Ax Cubed is equal to:\n{sp.simplify(Ax_cube)}\n\nNegative Ax is equal to:\n {-Ax}\n")
print(f"\nThe two end up being the same as you can take the -(ax^2, ay^2, az^2) term out of the matrix for Ax cubed. This term is equivalent to a*a.T which is equal to 1. ")


""" Problem 2 """
import numpy as np
# Constants
rvec = np.array([6783, 3391, 1953])
vvec = np.array([-3.5, 4.39, 4.44])
zlvlh = -rvec/np.linalg.norm(rvec)
ylvlh = -np.cross(rvec,vvec)/np.linalg.norm(np.cross(rvec,vvec))
xlvlh = np.cross(ylvlh, zlvlh)
Flvlh = np.array([[xlvlh], [ylvlh], [zlvlh]])
Feci = np.array([[1, 0, 0],
                 [0, 1, 0],
                 [0, 0, 1]])
C_lvlh2eci = Flvlh @ Feci.T
print(C_lvlh2eci)