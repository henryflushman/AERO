""" Dynamics and Controls Homework 1 """

# Package imports
import os
os.system('cls')
import numpy as np
import DynCon_Functions as dcf

""" Problem 1.1 """
print("Problem 1.1\n")

# Constants
roll = 2
pitch = -1
yaw = 40
r_b = np.array([[.1], [-.2], [.95]])

# Rotation matrices
C_lb = dcf.rotation_matrix(roll, pitch, yaw) # LVLH to body rotation matrix
print(f"Rotation matrix from LVLH to body, C_lvlh-b:\n{C_lb}\n")
print("We can relate the rotation matrix C_lvlh-b to r_lvlh through the equation:\n")
print("r_b = C_lvlh-b * r_lvlh\n")
print("Multiplying both sides by (C_lvlh-b) Transpose, we find the equation:\n")
print("(C_lvlh-b)^T * r_b = r_lvlh\n")
r_l = np.dot(C_lb.T, r_b)
print(f"Position vector in LVLH frame, r_lvlh:\n{r_l}\n")

""" Problem 1.2 """
print("Problem 1.2\n")

""" Problem 1.3 """
print("Problem 1.3\n")


""" Problem 2 """
print("Problem 2\n")

x1 = np.array([[1], [0], [0]])
y1 = np.array([[0], [1], [0]])
z1 = np.array([[0], [0], [1]])
b1 = [x1, y1, z1]

print("A.\n")
print(f"The component vector of each basis vector x1, y1, z1 in F1 components is:\nx1:{x1}\ny1:{y1}\nz1:{z1}\n")
print("B.\n")

C_12 = dcf.rotation_matrix(0, 0, 25)
b1_in_2 = np.dot(C_12.T, b1)
print(f"The component vector of each basis vector x1, y1, z1 in F2 components is:\nx1:{b1_in_2[0]}\ny1:{b1_in_2[1]}\nz1:{b1_in_2[2]}\n")

print("C.\n")
print(f"The rotation matrix, Cz: \n{C_12}\n")
print("This rotation matrix is related to the basis matrix, b1, in F2 such that:\n")
print("b1_in_2 = Cz^T * b1\n")

""" Problem 3 """
print("Problem 3\n")