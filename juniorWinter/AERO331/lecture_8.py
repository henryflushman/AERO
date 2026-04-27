import numpy as np
import matplotlib.pyplot as plt
import libFRC_v1_1 as frc

ply = frc.Ply(
    E1 = 1.4e11,
    E2 = 1e10,
    G12 = 7e9,
    nu12 = 0.3,
    theta = 60,
    zL = -5e-4,
    zH = 5e-4
)

nu12 = ply.nu21
Q = ply.Q

def transformation_matrix_sigma(theta):
    T_mat = np.array([
        [np.cos(theta)**2, np.sin(theta)**2, 2*np.sin(theta)*np.cos(theta)],
        [np.sin(theta)**2, np.cos(theta)**2, -2*np.sin(theta)*np.cos(theta)],
        [-np.sin(theta)*np.cos(theta), np.sin(theta)*np.cos(theta), np.cos(theta)**2 - np.sin(theta)**2]
        ])
    return T_mat

def transformation_matrix_eps(theta):
    T_mat = np.array([
        [np.cos(theta)**2, np.sin(theta)**2, np.sin(theta)*np.cos(theta)],
        [np.sin(theta)**2, np.cos(theta)**2, -np.sin(theta)*np.cos(theta)],
        [-2*np.sin(theta)*np.cos(theta), 2*np.sin(theta)*np.cos(theta), np.cos(theta)**2 - np.sin(theta)**2]
        ])
    return T_mat

T_matrix_sigma = transformation_matrix_sigma(np.deg2rad(30))
T_matrix_eps = transformation_matrix_eps(np.deg2rad(30))


Q_bar = np.linalg.inv(T_matrix_sigma) @ Q @ T_matrix_eps

print(ply.Qbar)
print(np.linalg.inv(transformation_matrix_sigma(np.deg2rad(60))) @ Q @ transformation_matrix_eps(np.deg2rad(60)))