import DynCon_Functions as dc
import numpy as np
import sympy as sp
import os

os.system('cls')

""" Problem 1 """

F_ac = np.array([[1, 0, 0],
               [0, 1, 0],
               [0, 0, 1]])
F_enu = np.array([[1, 0, 0],
               [0, -1, 0],
               [0, 0, -1]])

C_enu_ac = F_ac.T @ F_enu

""" Question 1"""
print(f"Rotation matrix from ENU to AC, C_enu2ac:\n{C_enu_ac}\n")

""" Question 2 """
C_angle = dc.rotation_matrix(10, -5, 20)
C_ac0_acf = C_enu_ac @ C_angle

print(f"Rotation matrix from AC0 to ACf, C_ac02acf:\n{C_ac0_acf}\n")

# 2.1
yaw, pitch, roll = np.rad2deg(dc.rotm_to_euler321(C_enu_ac.T))

print(f"Yaw: {yaw}\nPitch: {pitch}\nRoll: {roll}\n")

# 2.2

r_enu = np.array([[0], [1], [0]])

C_acf_enu = C_enu_ac @ C_ac0_acf

print(f"Rotation matrix from ACf to ENU, C_acf2enu:\n{C_acf_enu}\n")

r_acf = C_enu_ac @ C_ac0_acf @ r_enu

print(f"North vector in ACf frame, r_acf:\n{r_acf}\n")

# 2.3
eta, eps = dc.quaternion_from_rotation_matrix(C_acf_enu)
eta = float(eta)
print(f"The associated quaternion of C_acf2enu has eta = {eta} and eps = {eps}")